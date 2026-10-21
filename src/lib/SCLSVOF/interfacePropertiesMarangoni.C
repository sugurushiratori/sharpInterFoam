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
    Marangoni (thermocapillary) force for the CSF model:
        F_M = sigmaTemp * (I - n n) & grad(T) * delta
    evaluated as a face-normal force density so that it can be added to
    surfaceTensionForce() in the momentum and pressure equations.

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesLS.H"
#include "fvCFD.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
