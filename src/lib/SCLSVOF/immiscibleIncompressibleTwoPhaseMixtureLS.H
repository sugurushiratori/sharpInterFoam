/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

Class
    Foam::immiscibleIncompressibleTwoPhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    immiscibleIncompressibleTwoPhaseMixtureLS.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleIncompressibleTwoPhaseMixtureLS_H
#define immiscibleIncompressibleTwoPhaseMixtureLS_H

#include "incompressibleTwoPhaseMixtureLS.H"
#include "interfacePropertiesLS.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class immiscibleIncompressibleTwoPhaseMixtureLS Declaration
\*---------------------------------------------------------------------------*/

class immiscibleIncompressibleTwoPhaseMixtureLS
:
    public incompressibleTwoPhaseMixtureLS,
    public interfacePropertiesLS
{
public:

    // Constructors

        //- Construct from components
        immiscibleIncompressibleTwoPhaseMixtureLS
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~immiscibleIncompressibleTwoPhaseMixtureLS() = default;


    // Member Functions
        
        //- Correct the transport and interface properties
        void correct()
        {
            interfacePropertiesLS::correctPsi();
            
            if(interfacePropertiesLS::densityFunctionHeaviside())
            {
                densityFunc_     == interfacePropertiesLS::Heaviside();
                densityFuncFace_ == interfacePropertiesLS::HeavisideFace();
            }
            else
            {
                densityFunc_     == max(scalar(0), min(scalar(1), alpha1() ));
                densityFuncFace_ == max(scalar(0), min(scalar(1), fvc::interpolate(alpha1()) ));
            }
            
            incompressibleTwoPhaseMixtureLS::correct();
        }
        
        void correctSurfaceTensionForce()
        {
            interfacePropertiesLS::correctCurvature();
        }
        
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
