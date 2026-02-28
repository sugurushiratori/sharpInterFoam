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


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfacePropertiesLS::initializePsi()
{
     Info << "initializing Level-Set function..";
     
     // Level-Set値の初期値設定
     if(params_.initializeAtanh)
     {
         Info << " (atanh) ";
         const scalar small = scalar(1e-6);
         const volScalarField limitedAlpha = max(small, min(scalar(1.0-small), alpha1_) );
         const scalar pi = Foam::constant::mathematical::pi;

         psi0_ == widthEpsilon() / pi * atanh(scalar(2)*limitedAlpha - scalar(1));
     }
     else 
     {
         Info << " (linear) ";
         psi0_ == (scalar(2) * alpha1_ - scalar(1)) * widthEpsilon() * scalar(1.5);
     }
     psi0_.correctBoundaryConditions();
     
     psi_ == psi0_;
     Info << " finished." << endl;
}

Foam::tmp<volScalarField>
Foam::interfacePropertiesLS::calSignFunc(const volScalarField& psi)
{
     const fvMesh& mesh  = psi.mesh();
     tmp<volScalarField> tsignPsi
     (
        new volScalarField
        (
          IOobject
          (
              "signPsi",
              psi.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedScalar("signPsi", dimVelocity, scalar(0)),
          zeroGradientFvPatchField<scalar>::typeName
        )
     );
     volScalarField& signPsi = tsignPsi.ref();

     const scalar eps = widthEpsilon().value();
     const scalar pi  = Foam::constant::mathematical::pi;
     
     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label cid = *iter;
         signPsi[cid] = psi0_[cid] / sqrt( sqr(psi0_[cid]) + sqr(eps) );
     }
     
/*
     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label ci = *iter;
         if     (psi0_[ci] < -eps) signPsi[ci] = scalar(-1);
         else if(psi0_[ci] >  eps) signPsi[ci] = scalar( 1);
         else                      signPsi[ci] = psi0_[ci] / eps + sin(pi * psi0_[ci] / eps) / pi;
     }

     const volScalarField gg = mag(fvc::grad(psi));
     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label ci = *iter;
         const scalar alp = max(scalar(1), gg[ci]);

         if     (psi[ci] < -alp*eps) signPsi[ci] = scalar(-1);
         else if(psi[ci] >  alp*eps) signPsi[ci] = scalar( 1);
         else                        signPsi[ci] = psi[ci] / (alp*eps) + sin(pi * psi[ci] / (alp*eps)) / pi;
     }
*/

     signPsi.correctBoundaryConditions();
     
     return tsignPsi;
}

Foam::tmp<volScalarField>
Foam::interfacePropertiesLS::forcingTerm(const volScalarField& psi)
{
     const fvMesh& mesh  = psi.mesh();
     tmp<volScalarField> tFF
     (
        new volScalarField
        (
          IOobject
          (
              "forcingTerm",
              psi.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedScalar("F", dimVelocity, scalar(0)),
          zeroGradientFvPatchField<scalar>::typeName
        )
     );
     volScalarField& FF = tFF.ref();
     
     volScalarField rr = psi0_;
     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label cid = *iter;
         scalar sum = 0.0;
         label count = 0;
         forAll(mesh.cellCells()[cid], ccid)
         {
              label ncid = mesh.cellCells()[cid][ccid];
              if(psi0_[ncid] * psi0_[cid] < 0.0)
              {
                  sum += psi0_[ncid];
                  count ++;
              }
         }
         if(count == 0)
             rr[cid] = 0;
         else 
             rr[cid] /= sum + VSMALL;
     }
     
     const dimensionedScalar eps = widthEpsilon();
     
     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label cid = *iter;
         scalar sum = 0.0;
         label  count = 0;
         forAll(mesh.cellCells()[cid], ccid)
         {
              label ncid = mesh.cellCells()[cid][ccid];
              if(psi[ncid] * psi[cid] < 0.0)
              {
                  sum += psi[ncid];
                  count ++;
              }
         }
         if(count == 0)
             FF[cid] = 0.0;
         else 
             FF[cid] = (rr[cid] * sum - psi[cid]) / eps.value();
     }
     

     FF.correctBoundaryConditions();
     
     return tFF;
}

void Foam::interfacePropertiesLS::psi_timeAdvanceRK3(volScalarField& psi, const dimensionedScalar& deltaTau)
{
     volScalarField psi1 = psi + deltaTau * psi_timeStepFunc(psi);
     volScalarField psi2 = 3.0/4.0 * psi + 1.0/4.0 * psi1 + 1.0/4.0*deltaTau * psi_timeStepFunc(psi1);
     volScalarField psi3 = 1.0/3.0 * psi + 2.0/3.0 * psi2 + 2.0/3.0*deltaTau * psi_timeStepFunc(psi2);
     
     psi = psi3;
}

void Foam::interfacePropertiesLS::psi_timeAdvanceRK3B(volScalarField& psi, const dimensionedScalar& deltaTau)
{
     volScalarField psi1 = psi + deltaTau * psi_timeStepFunc2(psi);
     volScalarField psi2 = 3.0/4.0 * psi + 1.0/4.0 * psi1 + 1.0/4.0*deltaTau * psi_timeStepFunc2(psi1);
     volScalarField psi3 = 1.0/3.0 * psi + 2.0/3.0 * psi2 + 2.0/3.0*deltaTau * psi_timeStepFunc2(psi2);
     
     psi = psi3;
}

Foam::tmp<volScalarField>
Foam::interfacePropertiesLS::psi_timeStepFunc(const volScalarField& psi)
{
     const fvMesh& mesh  = psi.mesh();
     tmp<volScalarField> tpsiDt
     (
        new volScalarField
        (
          IOobject
          (
              "psiDt",
              psi.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedScalar("psiDt", dimVelocity, scalar(0)),
          zeroGradientFvPatchField<scalar>::typeName
        )
     );     
     volScalarField& psiDt = tpsiDt.ref();
     
     
     volScalarField signFunc = calSignFunc(psi);
     volScalarField magGradPsi = calMagGradPsi(signFunc, psi);

     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
          const label cid = *iter;
          psiDt[cid] =  signFunc[cid] * ( scalar(1) - magGradPsi[cid]);
     }
     psiDt.correctBoundaryConditions();
     
     return tpsiDt;
}

Foam::tmp<volScalarField>
Foam::interfacePropertiesLS::psi_timeStepFunc2(const volScalarField& psi)
{
     const fvMesh& mesh  = psi.mesh();
     tmp<volScalarField> tpsiDt
     (
        new volScalarField
        (
          IOobject
          (
              "psiDt",
              psi.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedScalar("psiDt", dimVelocity, scalar(0)),
          zeroGradientFvPatchField<scalar>::typeName
        )
     );     
     volScalarField& psiDt = tpsiDt.ref();
     
     volScalarField FF = forcingTerm(psi);

     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
          const label cid = *iter;
          psiDt[cid] =  0.5 * FF[cid];
     }
     psiDt.correctBoundaryConditions();
     
     return tpsiDt;
}

void Foam::interfacePropertiesLS::solvePsi()
{
     Info << "solve the re-initialization equation of Level-Set function...";

     scalar resid0 = 1.0;
     scalar resid = 1.0;
     const scalar nsub = params_.denomDeltaTau;
     dimensionedScalar deltaTau( "deltaTau", dimTime, scalar( meshSize().value() / nsub) );
     const label numLoop = label(widthFactor() * nsub * params_.factorNumLoop);

     label corr = 0;
     for(corr = 0; corr < numLoop; corr++)
     {
          tmp<volScalarField> signPsi0 = calSignFunc(psi_);          
          tmp<volScalarField> magGradPsi = calMagGradPsi(signPsi0, psi_);
          
          resid = calResidualMagGradPsi(magGradPsi);
          if(corr == 0) resid0 = resid;
          
          if(resid / resid0 <  params_.tolMagGradPsi )
          {
               Info << " -- converged " ;
               break;
          }
          
//          Info << " (iter = " << corr << ", resid = " << resid << ")" << endl;

          psi_timeAdvanceRK3(psi_, deltaTau);
          psi_timeAdvanceRK3B(psi_, deltaTau);
          
          narrowBand_.update(widthEpsilon().value());
     }

     Info << " (iter = " << corr << ", resid / resid0 = " << resid << " / " << resid0 << ")";

     Info << " finished. (Max, Min) = (" << max(psi_).value() << ", " << min(psi_).value() << ")" << nl;
}

//----------------------------------------------
//- calculate 2-norm of  ( 1 - |grad(psi))| )
//----------------------------------------------
Foam::scalar Foam::interfacePropertiesLS::calResidualMagGradPsi
(
   const volScalarField& magGradPsi
) const
{
    scalar resid = 0;
    forAllConstIter(labelList, narrowBand_.C1(), iter)
    {
        const label ci = *iter;
        resid += pow(scalar(1) - magGradPsi[ci], 2);
    }
    
    reduce(resid, sumOp<scalar>());

    return sqrt(resid);
}


// ************************************************************************* //
