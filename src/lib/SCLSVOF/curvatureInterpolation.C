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
    Curvature for the density-scaled CSF model (paper Sec. 3.4.5): the
    curvature must be the one of the nearest interface point, because the
    density-scaled delta function shifts the surface force off the
    interface.

    1. The raw curvature K = -div(n) is smoothed in the narrow band by two
       passes of a weighted average over the cell and
       its face neighbours (weights 1 + cos(pi psi/eps)), which suppresses
       the cell-to-cell noise of the curvature of the cut cells.
    No projection to the nearest interface point (K - psi n.grad(K)) is
    applied: grad(K) is unreliable in cells next to walls and band edges,
    and in the Marangoni cavity it made the interface relax back to flat.

    The evaluation is local (no iteration), independent of the length scale
    and of the domain decomposition (coupled patches are included).

    Note: an earlier pseudo-time extension (upwind solution of
    dK/dtau + S(psi) n.grad(K) = 0) gave a good static Laplace pressure but
    destroyed the oscillation of a deforming droplet, because it propagates
    the noisy value of the cut cell across the band.

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesLS.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfacePropertiesLS::curvatureInterpolation()
{
    const fvMesh& mesh = psi_.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalar eps = widthEpsilon().value();
    const scalar pi = constant::mathematical::pi;

    // 1. Weighted smoothing within the band
    volScalarField w
    (
        IOobject
        (
            "curvatureWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    for (const label celli : narrowBand_.C3())
    {
        const scalar r = psi_[celli]/eps;
        if (mag(r) < 1)
        {
            w[celli] = 1 + Foam::cos(pi*r);
        }
    }
    w.correctBoundaryConditions();

    const label nSmooth = 2;

    for (label pass = 0; pass < nSmooth; ++pass)
    {
        K_.correctBoundaryConditions();
        const scalarField K0(K_.primitiveField());
        const scalarField& wi = w.primitiveField();

        scalarField num(wi*K0);
        scalarField den(wi);

        forAll(own, facei)
        {
            const label o = own[facei];
            const label n = nei[facei];
            num[o] += wi[n]*K0[n];  den[o] += wi[n];
            num[n] += wi[o]*K0[o];  den[n] += wi[o];
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatchScalarField& pK = K_.boundaryField()[patchi];

            if (pK.coupled())
            {
                const scalarField Knb(pK.patchNeighbourField());
                const scalarField wnb
                (
                    w.boundaryField()[patchi].patchNeighbourField()
                );
                const labelUList& fc = mesh.boundary()[patchi].faceCells();

                forAll(fc, i)
                {
                    num[fc[i]] += wnb[i]*Knb[i];
                    den[fc[i]] += wnb[i];
                }
            }
        }

        scalarField& Kref = K_.primitiveFieldRef();

        for (const label celli : narrowBand_.C3())
        {
            // includes band cells with |psi| >= eps next to weighted cells
            if (den[celli] > SMALL)
            {
                Kref[celli] = num[celli]/den[celli];
            }
        }
    }

    K_.correctBoundaryConditions();
}


// ************************************************************************* //
