/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd
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

#include "NewtonianThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(NewtonianThermo, 0);
    addToRunTimeSelectionTable(viscosityModel, NewtonianThermo, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::NewtonianThermo::NewtonianThermo
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    Newtonian(name, viscosityProperties, U, phi),
    kappa0_("kappa", dimViscosity, viscosityProperties_),
    kappa_
    (
        IOobject
        (
            name + "_kappa",    // must differ from the name of nu_
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        U_.mesh(),
        kappa0_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::NewtonianThermo::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);
    Newtonian::read(viscosityProperties);

    viscosityProperties_.readEntry("kappa", kappa0_);
    kappa_ = kappa0_;

    return true;
}


// ************************************************************************* //
