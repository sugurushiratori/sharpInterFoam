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


#include "interfacePropertiesLS.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"
#include "CMULES.H"
#include "gaussConvectionScheme.H"
#include "surfaceInterpolate.H"
#include "upwind.H"
#include "downwind.H"
#include "unitConversion.H"

void Foam::interfacePropertiesLS::K_timeAdvanceRK3(volScalarField& K, const dimensionedScalar& deltaTau)
{
     volScalarField K1 = K + deltaTau * K_timeStepFunc(K);
     volScalarField K2 = 3.0/4.0 * K + 1.0/4.0 * K1 + 1.0/4.0*deltaTau * K_timeStepFunc(K1);
     volScalarField K3 = 1.0/3.0 * K + 2.0/3.0 * K2 + 2.0/3.0*deltaTau * K_timeStepFunc(K2);
     
     K = K3;
}

Foam::tmp<volScalarField>
Foam::interfacePropertiesLS::K_timeStepFunc(const volScalarField& K)
{
     const fvMesh& mesh  = K.mesh();
     tmp<volScalarField> tKDt
     (
        new volScalarField
        (
          IOobject
          (
              "KDt",
              K.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedScalar("KDt", dimless / dimLength, scalar(0)),
          zeroGradientFvPatchField<scalar>::typeName
        )
     );     
     volScalarField& KDt = tKDt.ref();
     
     surfaceScalarField nPsi = nHatf_ * fvc::interpolate(psi_);

     surfaceScalarField Kf = fvc::interpolate(K);

     KDt = fvc::div(nPsi * Kf);
     /*
     volScalarField signFunc = calSignFunc(psi_);
     volVectorField gradK = fvc::grad(K);
     volScalarField ngK = nHatv_ & gradK;

     forAll(KDt, ci)
     {
        KDt[ci] = signFunc[ci] * ngK[ci];
     }
     */
     KDt.correctBoundaryConditions();
     
     return tKDt;
}

void Foam::interfacePropertiesLS::curvatureInterpolation()
{
    Info << "interpolate the curvature for density-scaled CSF...";
    const fvMesh& mesh  = psi_.mesh();

    const scalar nsub = params_.denomDeltaTau;
    dimensionedScalar deltaTau( "deltaTau", dimless, scalar( meshSize().value() / nsub) );
    const label numLoop = label(widthFactor() * nsub * params_.factorNumLoop);

    label corr = 0;
    for(corr = 0; corr < numLoop; corr++)
    {         
        K_timeAdvanceRK3(K_, deltaTau);
    }

    K_.correctBoundaryConditions();
    
    Info << " finished." << nl;
}


// ************************************************************************* //
