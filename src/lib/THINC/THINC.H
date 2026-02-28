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

Class
    Foam::THINC

Description
    THINC (tangent of hyperbora for interface capturing) scheme
    implemented as surfaceInterpolationScheme.

Author
    Suguru SHIRATORI

SourceFiles
    THINC.C

\*---------------------------------------------------------------------------*/

#ifndef THINC_H
#define THINC_H

#include "fvMesh.H"
#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class THINC
:
     public refCount
{

    // Private data

        const fvMesh& mesh_;
        
        //- faceFlux field
        const surfaceScalarField& phi_;

        //- level-set function
        const volScalarField& psi_;
        
        //- unit vector normal to interface
        const volVectorField& nHatv_;
        
        const dictionary& thincDict_;
        
        // parameters
        //- parameter for interface width
              bool   enabled_;
        const scalar beta0_;
              scalar epsAlpha_;
              scalar epsBeta_;

    // Private Member Functions
    
        //- No copy construct
        THINC(const THINC&) = delete;
        
        //- No copy assignment
        void operator=(const THINC&) = delete;

        //- THINC core function (calculate flux)
        scalar flux1D
        (
            const scalar& fP,
            const scalar& fN,
            const scalar& normal,
            const scalar& beta,
            const scalar& phif,
            const scalar& dt,
            const scalar& dx
        ) const;
        
public:

    // Constructors

        //- Construct from faceFlux
        THINC
        (
            const fvMesh& mesh, 
            const surfaceScalarField& faceFlux,
            const volScalarField& psi,
            const volVectorField& nHatv,
            const dictionary& dict,
            const scalar beta0
        );
        

    // Destructor
        virtual ~THINC() = default;

    // Member Functions

        tmp<surfaceScalarField>
        getFlux(const volScalarField& vf) const;
        
        
        bool enabled() const {return enabled_;}
        
        bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
