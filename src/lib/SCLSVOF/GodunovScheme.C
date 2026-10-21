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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gaussConvectionScheme.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "upwind.H"
#include "downwind.H"
#include "interfacePropertiesLS.H"



void Foam::interfacePropertiesLS::godunovDifferences
(
    const volScalarField& psi
) const
{
    // One-sided (upwind and downwind) differences of psi in the x, y and z
    // directions, evaluated as Gauss sums over the faces:
    //     GU_d = (1/V) sum_f Sf_d psi_f(upwind w.r.t. e_d)
    //     GD_d = (1/V) sum_f Sf_d psi_f(downwind w.r.t. e_d)
    // This is the same as fvc::div(e_d & Sf * upwind/downwind interpolate)
    // but uses one loop over the faces without temporary fields.
    // Values are valid in the band cells C3 only.

    const fvMesh& mesh = psi.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& V = mesh.V();
    const scalarField& psiI = psi.primitiveField();

    // Only the band cells C3 are needed: faces without a band cell on
    // either side are skipped
    const labelList& band = narrowBand_.C3();
    godunovMark_.setSize(mesh.nCells());
    godunovMark_ = false;
    for (const label ci : band)
    {
        godunovMark_[ci] = true;
    }

    for (direction d = 0; d < 6; ++d)
    {
        godunovG_[d].setSize(mesh.nCells());
        for (const label ci : band)
        {
            godunovG_[d][ci] = 0;
        }
    }

    forAll(own, facei)
    {
        const label o = own[facei];
        const label n = nei[facei];
        if (!godunovMark_[o] && !godunovMark_[n])
        {
            continue;
        }
        const scalar po = psiI[o];
        const scalar pn = psiI[n];

        for (direction d = 0; d < 3; ++d)
        {
            const scalar s = Sf[facei][d];
            // upwind w.r.t. e_d: owner value for s >= 0
            const scalar fu = s*(s >= 0 ? po : pn);
            const scalar fd = s*(s >= 0 ? pn : po);
            godunovG_[d][o] += fu;
            godunovG_[d][n] -= fu;
            godunovG_[d + 3][o] += fd;
            godunovG_[d + 3][n] -= fd;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatchScalarField& pp = psi.boundaryField()[patchi];
        if (!pp.size())
        {
            continue;
        }
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        if (pp.coupled())
        {
            const scalarField pI(pp.patchInternalField());
            const scalarField pN(pp.patchNeighbourField());
            forAll(fc, i)
            {
                if (!godunovMark_[fc[i]]) continue;
                for (direction d = 0; d < 3; ++d)
                {
                    const scalar s = pSf[i][d];
                    godunovG_[d][fc[i]] += s*(s >= 0 ? pI[i] : pN[i]);
                    godunovG_[d + 3][fc[i]] += s*(s >= 0 ? pN[i] : pI[i]);
                }
            }
        }
        else
        {
            forAll(fc, i)
            {
                if (!godunovMark_[fc[i]]) continue;
                for (direction d = 0; d < 3; ++d)
                {
                    const scalar f = pSf[i][d]*pp[i];
                    godunovG_[d][fc[i]] += f;
                    godunovG_[d + 3][fc[i]] += f;
                }
            }
        }
    }

    for (direction d = 0; d < 6; ++d)
    {
        for (const label ci : band)
        {
            godunovG_[d][ci] /= V[ci];
        }
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesLS::calMagGradPsi
(
    const volScalarField& signPsi0,
    const volScalarField& psi
) const
{
    // |grad(psi)| with the Godunov Hamiltonian, in the narrow band C3

    const fvMesh& mesh = psi_.mesh();

    auto tmagGradPsi = tmp<volScalarField>::New
    (
        IOobject
        (
            "magGradPsi",
            psi.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(psi.dimensions()/dimLength, Zero),
        zeroGradientFvPatchField<scalar>::typeName
    );
    volScalarField& magGradPsi = tmagGradPsi.ref();

    godunovDifferences(psi);

    const scalarField& GxU = godunovG_[0];
    const scalarField& GyU = godunovG_[1];
    const scalarField& GzU = godunovG_[2];
    const scalarField& GxD = godunovG_[3];
    const scalarField& GyD = godunovG_[4];
    const scalarField& GzD = godunovG_[5];

    for (const label ci : narrowBand_.C3())
    {
        if (signPsi0[ci] > 0.0)
        {
            const scalar vx = max(sqr(max(scalar(0), GxU[ci])), sqr(min(scalar(0), GxD[ci])));
            const scalar vy = max(sqr(max(scalar(0), GyU[ci])), sqr(min(scalar(0), GyD[ci])));
            const scalar vz = max(sqr(max(scalar(0), GzU[ci])), sqr(min(scalar(0), GzD[ci])));

            magGradPsi[ci] = sqrt(vx + vy + vz);
        }
        else
        {
            const scalar vx = max(sqr(min(scalar(0), GxU[ci])), sqr(max(scalar(0), GxD[ci])));
            const scalar vy = max(sqr(min(scalar(0), GyU[ci])), sqr(max(scalar(0), GyD[ci])));
            const scalar vz = max(sqr(min(scalar(0), GzU[ci])), sqr(max(scalar(0), GzD[ci])));

            magGradPsi[ci] = sqrt(vx + vy + vz);
        }
    }

    magGradPsi.correctBoundaryConditions();

    return tmagGradPsi;
}


// ************************************************************************* //
