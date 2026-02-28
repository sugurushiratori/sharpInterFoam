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

\*---------------------------------------------------------------------------*/

#include "narrowBand.H"
#include "fvc.H"

Foam::narrowBand::narrowBand
(
     const volScalarField& psi
)
 :
 psi_(psi),
 band1_(0),
 band2_(0),
 band3_(0)
{
}


void Foam::narrowBand::update(const scalar widthEpsilon)
{
    band1_.clear();
    band2_.clear();
    band3_.clear();
    
    DynamicList<label> list1, list2, list3;
    
    const scalar eps1 = widthEpsilon * 1.0;
    const scalar eps2 = widthEpsilon * 1.5;
    const scalar eps3 = widthEpsilon * 2.0;
    
    forAll(psi_, ci)
    {
        if( mag(psi_[ci]) <= eps1 )  list1.append(ci);
        if( mag(psi_[ci]) <= eps2 )  list2.append(ci);
        if( mag(psi_[ci]) <= eps3 )  list3.append(ci);
    }
    
    band1_.transfer(list1);
    band2_.transfer(list2);
    band3_.transfer(list3);

//    band1_ = list1.xfer();
//    band2_ = list2.xfer();
//    band3_ = list3.xfer();
}

void Foam::narrowBand::showSizes() const
{
     label num1 = band1_.size();
     label num2 = band2_.size();
     label num3 = band3_.size();
     reduce(num1, sumOp<label>());
     reduce(num2, sumOp<label>());
     reduce(num3, sumOp<label>());
          
     label nCells = psi_.mesh().nCells();
     reduce(nCells, sumOp<label>());
     
     Info << " (narrow-band: C1, C2, C3, All = " << num1 << ", "<< num2 << ", " << num3 <<", " << nCells  <<")" << endl;
}

Foam::label Foam::narrowBand::findNearestCell(const point& location) const
{
    const fvMesh& mesh = psi_.mesh();
    const vectorField& centres = mesh.cellCentres();
    
    label  nearestCelli = 0;
    scalar minProximity = magSqr(centres[band2_.first()] - location);
    
    forAllConstIter(labelList, band2_, iter)
    {
         const  label ci  = *iter;
         scalar proximity = magSqr(centres[ci] - location);
         if(proximity < minProximity)
         {
              nearestCelli = ci;
              minProximity = proximity;
         }
    }
    return nearestCelli;
}

Foam::label Foam::narrowBand::findCell(const point& location) const
{
    const fvMesh& mesh = psi_.mesh();
    
    const label num = band2_.size();
    if(num == 0) return -1;
    
    label  celli = findNearestCell(location);
    
    if(mesh.pointInCell(location, celli))
    {
        return celli;
    }
    else
    {
        bool cellFound = false;
        forAllConstIter(labelList, band2_, iter)
        {
            const label ci  = *iter;
            if(mesh.pointInCell(location, ci) )
            {
                cellFound = true;
                celli = ci;
                break;
            }
        }
        
        if (cellFound)
        {
            return celli;
        }
        else
        {
            return -1;
        }
    }
}

void Foam::narrowBand::addDummy()
{
    // CPUによって補間対象セルの数が大きく異なる場合にエラーになるので
    // 数を合わせるための処置 (もっと良い方法ありそう・・)
    
    const label nCells = psi_.mesh().nCells();
    if (nCells == 0) return;
    
    label num = band2_.size();
    label max = num;
    reduce(max, maxOp<label>());
     
    DynamicList<label> list;
    
    forAllConstIter(labelList, band2_, iter)
    {
        list.append(*iter); 
    }
    
    label cid = 0;
    while(num < max) {list.append(cid % nCells); num++; cid++;}

    band2_.clear();
    band2_.transfer(list);
    
//    band2_ = list.xfer();

}

// ************************************************************************* //
