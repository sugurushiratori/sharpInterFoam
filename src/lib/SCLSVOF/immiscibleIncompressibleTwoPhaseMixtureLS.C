/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "immiscibleIncompressibleTwoPhaseMixtureLS.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immiscibleIncompressibleTwoPhaseMixtureLS::
immiscibleIncompressibleTwoPhaseMixtureLS
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixtureLS(U, phi),            // densityFunc is alpha at this point
    interfacePropertiesLS(alpha1(), U, *this)
{
    correct();
    correctSurfaceTensionForce();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::immiscibleIncompressibleTwoPhaseMixtureLS::read()
{
    return
        incompressibleTwoPhaseMixtureLS::read()
     && interfacePropertiesLS::read();
}

// ************************************************************************* //
