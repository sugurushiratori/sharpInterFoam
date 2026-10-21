/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Extension of the interface curvature off the interface for the
    density-scaled CSF model (paper Eq. 46):

        dK/dtau + S(psi) n & grad(K) = 0,   S(psi) = psi/sqrt(psi^2 + h^2)

    The characteristic velocity S(psi) n points away from the interface on
    both sides, so the curvature evaluated at the interface is transported
    outward along the normal. The equation is discretised with first-order
    upwinding in non-conservative form, using a local pseudo time step
    (Jacobi iteration):

        K_i <- K_i + c * sum_in( |F_f| (K_nb - K_i) ) / sum_in |F_f|,

    where F_f = S_f (n_f & S_f) and the sums run over the inflow faces of
    cell i. Every update is a convex combination (bounded, monotone) and the
    result is independent of the length scale of the problem.
    Coupled (processor/cyclic) patches are handled through the patch
    neighbour values; other boundaries contribute no inflow.

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesLS.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfacePropertiesLS::curvatureInterpolation()
{
    const fvMesh& mesh = psi_.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const scalar h = meshSize().value();

    // Pseudo-time CFL number (<= 1 for a convex update)
    const scalar cfl = 0.5;

    // Number of iterations: information has to travel over the band where
    // the (scaled) delta function is non-zero, |psi| < eps, i.e.
    // eps/h = widthFactor cells. Two-fold margin, cfl cells per iteration.
    // (Computed from widthFactor, not eps/h, so that round-off in the
    // globally reduced mesh size cannot change the count between serial
    // and parallel runs.)
    const label nIter =
        label(std::ceil(2.0*params_.widthFactor/cfl - SMALL)) + 1;

    // Face flux of the characteristic velocity S(psi) n
    const surfaceScalarField psif(fvc::interpolate(psi_));
    const surfaceScalarField F
    (
        nHatf_*psif/sqrt(sqr(psif) + dimensionedScalar(sqr(dimLength), sqr(h)))
    );
    const scalarField& Fi = F.primitiveField();

    // Cells to update: band where the scaled delta function is non-zero
    // (with margin)
    boolList active(mesh.nCells(), false);
    for (const label celli : narrowBand_.C3())
    {
        active[celli] = true;
    }

    // Sum of inflow face fluxes per cell (constant during the iteration)
    scalarField sumIn(mesh.nCells(), Zero);

    forAll(own, facei)
    {
        const scalar f = Fi[facei];
        if (f < 0) sumIn[own[facei]] -= f;     // enters owner
        else       sumIn[nei[facei]] += f;     // enters neighbour
    }

    forAll(mesh.boundary(), patchi)
    {
        if (K_.boundaryField()[patchi].coupled())
        {
            const fvsPatchScalarField& pF = F.boundaryField()[patchi];
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(pF, facei)
            {
                if (pF[facei] < 0) sumIn[fc[facei]] -= pF[facei];
            }
        }
    }

    scalarField dK(mesh.nCells());

    for (label iter = 0; iter < nIter; ++iter)
    {
        K_.correctBoundaryConditions();
        const scalarField& K = K_.primitiveField();

        dK = Zero;

        forAll(own, facei)
        {
            const scalar f = Fi[facei];
            const label o = own[facei];
            const label n = nei[facei];

            if (f < 0)
            {
                dK[o] -= f*(K[n] - K[o]);
            }
            else
            {
                dK[n] += f*(K[o] - K[n]);
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatchScalarField& pK = K_.boundaryField()[patchi];

            if (pK.coupled())
            {
                const fvsPatchScalarField& pF = F.boundaryField()[patchi];
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const scalarField Knb(pK.patchNeighbourField());

                forAll(pF, facei)
                {
                    if (pF[facei] < 0)
                    {
                        const label c = fc[facei];
                        dK[c] -= pF[facei]*(Knb[facei] - K[c]);
                    }
                }
            }
        }

        scalarField& Kref = K_.primitiveFieldRef();

        forAll(Kref, celli)
        {
            if (active[celli] && sumIn[celli] > VSMALL)
            {
                Kref[celli] += cfl*dK[celli]/sumIn[celli];
            }
        }
    }

    K_.correctBoundaryConditions();
}


// ************************************************************************* //
