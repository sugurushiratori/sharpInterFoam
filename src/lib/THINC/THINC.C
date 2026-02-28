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

\*---------------------------------------------------------------------------*/

#include "THINC.H"
#include "fvCFD.H"
#include "coupledFvPatchField.H"
#include "fvc.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//---------------------------------------------
//- Constructor
//---------------------------------------------
Foam::THINC::THINC
(
    const fvMesh& mesh, 
    const surfaceScalarField& faceFlux,
    const volScalarField& psi,
    const volVectorField& nHatv,
    const dictionary& dict,
    const scalar beta0
) 
:
    mesh_(mesh),
    phi_(faceFlux),
    psi_(psi),
    nHatv_(nHatv),
    thincDict_(dict),
    
    enabled_(thincDict_.getOrDefault<bool>("enable", true) ),
    beta0_(beta0),
    epsAlpha_(thincDict_.getOrDefault<scalar>("epsAlpha", 1e-10) ),
    epsBeta_(thincDict_.getOrDefault<scalar>("epsBeta", 1e-4) )
{
    read();
}



bool Foam::THINC::read()
{
    enabled_  = thincDict_.lookupOrDefault<Switch>("enable", true);
    epsAlpha_ = thincDict_.lookupOrDefault<scalar>("epsAlpha", 1e-10);
    epsBeta_  = thincDict_.lookupOrDefault<scalar>("epsBeta", 1e-4);

    return true;
}


//---------------------------------------------
//- THINC core function
//  fP:      owner側のalpha
//  fN:      neighbour側のalpha
//  normal:  法線ベクトルのフェイス射影
//  phif:    速度のフェイス射影
//  dt:      時間刻み
//  dx:      格子幅
//---------------------------------------------
Foam::scalar
Foam::THINC::flux1D
(
     const scalar& fP,
     const scalar& fN,
     const scalar& normal,
     const scalar& betaArg,
     const scalar& phif,
     const scalar& dt,
     const scalar& dx
) const
{
      const scalar fup    = (phif >= 0.0 ? fP        : fN         );  // upwind value of alpha
      const scalar xface  = (phif >= 0.0 ? scalar(1) : scalar(0)  );  // face position in the upwind cell
      const scalar gamma  = ( normal >= 0.0 ? scalar(1) : scalar(-1) );  // direction of interface

      if( fup < epsAlpha_ || scalar(1) - epsAlpha_ < fup)
      {
         return fup * phif * dt;
      }
      else
      {
//         if(betaArg < epsBeta_) return fup * phif * dt;
         const scalar beta  = (betaArg < epsBeta_ ? epsBeta_ : betaArg);
         
         const scalar a1    = ::exp(beta/gamma * (scalar(2) * fup - scalar(1)) );
         const scalar a3    = ::exp(beta);
         const scalar dd    = scalar(1)/(scalar(2)*beta) * ::log( (a3*a3 -a1*a3) / (a1*a3 -scalar(1)) );
         const scalar a4    = ::cosh(beta*(xface -(phif * dt/dx) -dd));
         const scalar a5    = ::cosh(beta*(xface - dd));
         const scalar flux  = scalar(1)/scalar(2)*( phif * dt -gamma*dx/beta * ::log(a4/a5));
         return flux;
      }
}

// ************************************************************************* //
