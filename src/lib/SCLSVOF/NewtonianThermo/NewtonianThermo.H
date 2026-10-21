/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

Class
    Foam::viscosityModels::Newtonian

Description
    Newtonian viscosity model that additionally reads the thermal diffusivity
    "kappa" [m^2/s]. Kept for backward compatibility of case files; the
    mixture reads "kappa" directly from the phase dictionaries.

SourceFiles
    NewtonianThermo.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonianThermo_H
#define NewtonianThermo_H

#include "Newtonian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Newtonian Declaration
\*---------------------------------------------------------------------------*/

class NewtonianThermo
:
    public Newtonian
{
    // Private data

        dimensionedScalar kappa0_;


public:

    // Protected data

        // thermal diffusivity
        volScalarField kappa_;


    //- Runtime type information
    TypeName("NewtonianThermo");


    // Constructors

        //- Construct from components
        NewtonianThermo
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~NewtonianThermo() = default;


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> kappa() const
        {
            return kappa_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> kappa(const label patchi) const
        {
            return kappa_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity (not appropriate, viscosity constant)
        virtual void correct()
        {}

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
