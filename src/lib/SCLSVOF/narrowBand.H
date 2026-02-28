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
    Foam::narrowBand

Description
    calculate narrow-bands for Level-Set calculation

SourceFiles
    narrowBand.C

\*---------------------------------------------------------------------------*/

#ifndef narrowBand_H
#define narrowBand_H

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class narrowBand
{
    // Private data
    
        const volScalarField& psi_;

        labelList band1_;
        labelList band2_;
        labelList band3_;
        
    // Private Member Functions

        
public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        narrowBand
        (
            const volScalarField& psi
        );


    // Member Functions

    void update(const scalar widthEpsilon);
    
    void showSizes() const;
    
    const labelList& C1() const { return band1_; }
    const labelList& C2() const { return band2_; }
    const labelList& C3() const { return band3_; }

    
    label findNearestCell(const point& p) const;
    label findCell(const point& p) const;
    
    void addDummy();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
