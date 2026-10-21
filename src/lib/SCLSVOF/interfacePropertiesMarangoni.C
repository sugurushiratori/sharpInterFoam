/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2026 Suguru Shiratori, Tokyo City University
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
    - sigmaK(): normal force sigma(T) K n delta (used in
      surfaceTensionForce())
    - marangoniForce(): tangential (thermocapillary) force
        F_M = sigmaTemp * (I - n n) & grad(T) * delta
      evaluated as a face-normal force density so that it can be added to
      surfaceTensionForce() in the momentum and pressure equations.

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesLS.H"
#include "fvCFD.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesLS::sigmaK() const
{
    // Normal surface-tension force with sigma(T). Together with the
    // tangential Marangoni force (marangoniForce) this is the full
    // surface force div_s(sigma (I - nn)) = sigma K n + grad_s(sigma);
    // with a constant sigma in the normal part the integral of the
    // Marangoni force over a closed interface is not balanced and a drop
    // in a temperature gradient is driven by a spurious net force.
    const volScalarField* TPtr =
        psi_.mesh().findObject<volScalarField>("T");

    if (!TPtr || mag(sigmaTemp_.value()) < VSMALL)
    {
        return sigmaPtr_->sigma()*K_;
    }

    return (sigmaPtr_->sigma() + sigmaTemp_*(*TPtr - TRef_))*K_;
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

    // Tangential projection of grad(T) at cell centres
    const volVectorField gradT(fvc::grad(T));
    const volVectorField projectedGradT((I - sqr(nHatv_)) & gradT);

    // Face-normal component
    const surfaceScalarField pgtf
    (
        fvc::interpolate(projectedGradT) & (mesh.Sf()/mesh.magSf())
    );

    if (params_.densityScaled)
    {
        maf = sigmaTemp_*pgtf*scaledDeltaFace_;
    }
    else
    {
        maf = sigmaTemp_*pgtf*DeltaFace_;
    }

    return tmaf;
}


// ************************************************************************* //
