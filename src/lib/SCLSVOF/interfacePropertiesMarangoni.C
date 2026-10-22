/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2026 Suguru Shiratori
License
    This file is part of sharpInterFoam, a derivative of OpenFOAM.

    sharpInterFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    sharpInterFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with sharpInterFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Temperature-dependent surface tension
        sigma(T) = sigma + sigmaTemp*(T - TRef)       (TRef default 0)
    for the CSF model:
    - sigmaK(): sigma(Tmean) K with the mean interface temperature Tmean
      (used in surfaceTensionForce(), balanced with the pressure gradient)
    - marangoniForce(): the deviation sigma(T) - sigma(Tmean), normal and tangential
      (thermocapillary), as the divergence of a surface stress, evaluated
      as a face-normal force density so that it can be added to
      surfaceTensionForce() in the momentum and pressure equations.

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesLS.H"
#include "fvCFD.H"
#include "symmTransformField.H"
#include "linear.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::interfacePropertiesLS::interfaceMeanT
(
    const volScalarField& T
) const
{
    // Mean temperature on the interface(s), weighted by the delta function
    const volScalarField& H =
        params_.densityScaled ? scaledHeaviside_ : Heaviside_;
    const scalarField w(mag(fvc::grad(H))().primitiveField()*psi_.mesh().V());
    const scalar sw = gSum(w);
    return sw > VSMALL ? gSum(w*T.primitiveField())/sw : 0;
}


Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesLS::sigmaK() const
{
    // Surface tension at the mean interface temperature, sigma(Tmean),
    // balanced with the pressure gradient. The deviation
    // sigma(T) - sigma(Tmean) is in marangoniForce().
    const volScalarField* TPtr =
        psi_.mesh().findObject<volScalarField>("T");

    if (!TPtr || mag(sigmaTemp_.value()) < VSMALL)
    {
        return sigmaPtr_->sigma()*K_;
    }

    const dimensionedScalar Tmean(dimTemperature, interfaceMeanT(*TPtr));

    return (sigmaPtr_->sigma() + sigmaTemp_*(Tmean - TRef_))*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesLS::marangoniForce() const
{
    const fvMesh& mesh = psi_.mesh();

    auto tmaf = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "marangoniForce",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), Zero)
    );
    surfaceScalarField& maf = tmaf.ref();
    maf.setOriented();

    // Isothermal problem: no temperature field or zero coefficient
    const volScalarField* TPtr = mesh.findObject<volScalarField>("T");

    if (!TPtr || mag(sigmaTemp_.value()) < VSMALL)
    {
        return tmaf;
    }

    const volScalarField& T = *TPtr;

    // Temperature-dependent part of the surface force in conservative form
    //     f = div( (sigma(T) - sigma0) (I - n n) delta )
    //       = (sigma(T) - sigma0) K n delta + grad_s(sigma(T)) delta
    // with sigma0 = sigma(Tmean) (Tmean: mean interface temperature), so
    // that sigma(T) - sigma0 = sigmaTemp*(T - Tmean). The constant part
    // sigma0 K n delta is in surfaceTensionForce() (balanced with the
    // pressure gradient); choosing sigma0 at the interface keeps the part
    // that is not balanced small. Written as the divergence
    // of a surface stress that vanishes away from the interface, the
    // discrete net force over a closed interface is zero (as it is for the
    // exact solution); a separate normal and tangential discretisation
    // left a spurious net force that drives drops in a temperature
    // gradient. delta = |grad(H)| (|grad(scaled H)| for the density-scaled
    // model), consistent with surfaceTensionForce().
    const volScalarField& H =
        params_.densityScaled ? scaledHeaviside_ : Heaviside_;

    const volScalarField deltaC(mag(fvc::grad(H)));

    const volSymmTensorField stress
    (
        IOobject
        (
            "marangoniStress",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        sigmaTemp_*(T - dimensionedScalar(dimTemperature, interfaceMeanT(T)))
       *deltaC*(symmTensor::I - sqr(nHatv_))
    );

    // Gauss divergence with linear face interpolation (independent of the
    // schemes in fvSchemes)
    volVectorField f
    (
        IOobject
        (
            "marangoniForceDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        fvc::surfaceIntegrate(mesh.Sf() & linearInterpolate(stress)),
        extrapolatedCalculatedFvPatchVectorField::typeName
    );
    f.correctBoundaryConditions();

    if (!params_.densityScaled)
    {
        maf = linearInterpolate(f) & (mesh.Sf()/mesh.magSf());
        return tmaf;
    }

    // Density-scaled model: the force per unit delta, interpolated between
    // the two cells weighted with their delta, is applied with the
    // density-scaled face delta (as in surfaceTensionForce()). A linear
    // interpolation of the cell force would put force on faces of the
    // light phase, where it drives large spurious velocities
    // (thermocapillary cavity, gas density 1e-3: max|U| 1.4 instead of 0.3).
    const scalarField& dI = deltaC.primitiveField();
    const scalar small = 1e-12*gMax(dI) + VSMALL;
    const vectorField g(f.primitiveField()/max(dI, small));
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& magSf = mesh.magSf().primitiveField();
    scalarField& mI = maf.primitiveFieldRef();

    forAll(own, facei)
    {
        const label o = own[facei];
        const label n = nei[facei];
        const scalar w = dI[o] + dI[n];
        if (w > small)
        {
            const vector gf((dI[o]*g[o] + dI[n]*g[n])/w);
            mI[facei] = (gf & Sf[facei])/magSf[facei]*scaledDeltaFace_[facei];
        }
    }

    // Coupled patches: same weighting with the neighbour-cell values
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& pp = mesh.boundary()[patchi];
        if (!pp.coupled())
        {
            continue;
        }
        const labelUList& fc = pp.faceCells();
        const scalarField dN(deltaC.boundaryField()[patchi].patchNeighbourField());
        const vectorField fN(f.boundaryField()[patchi].patchNeighbourField());
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];
        const scalarField& pDelta = scaledDeltaFace_.boundaryField()[patchi];
        scalarField& pm = maf.boundaryFieldRef()[patchi];

        forAll(fc, i)
        {
            const scalar dO = dI[fc[i]];
            const scalar w = dO + dN[i];
            if (w > small)
            {
                const vector gf((dO*g[fc[i]] + fN[i])/w);
                pm[i] = (gf & pSf[i])/pMagSf[i]*pDelta[i];
            }
        }
    }

    return tmaf;
}


// ************************************************************************* //
