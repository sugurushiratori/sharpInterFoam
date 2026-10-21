/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
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
     
     // Initial value of the level-set function
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
Foam::interfacePropertiesLS::calSignFunc() const
{
     // Smoothed sign of psi0: S = psi0/sqrt(psi0^2 + eps^2).
     // Evaluated in all cells (only band cells are used), so that it can be
     // computed once per re-initialisation although the band changes.
     const scalar eps = widthEpsilon().value();

     auto tsignPsi = tmp<volScalarField>::New
     (
         IOobject
         (
             "signPsi",
             psi0_.time().timeName(),
             psi0_.mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             IOobject::NO_REGISTER
         ),
         psi0_.mesh(),
         dimensionedScalar(dimVelocity, Zero),
         zeroGradientFvPatchField<scalar>::typeName
     );
     volScalarField& signPsi = tsignPsi.ref();

     const scalarField& p0 = psi0_.primitiveField();
     signPsi.primitiveFieldRef() = p0/sqrt(sqr(p0) + sqr(eps));
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

     // Sum and number of face neighbours with the opposite sign of f,
     // including neighbours across coupled (processor/cyclic) patches so
     // that the result does not depend on the domain decomposition.
     const labelUList& own = mesh.owner();
     const labelUList& nei = mesh.neighbour();

     auto oppositeSignSum = [&]
     (
         const volScalarField& f,
         scalarField& sum,
         labelList& cnt
     )
     {
         sum = Zero;
         cnt = Zero;
         const scalarField& fi = f.primitiveField();

         forAll(own, facei)
         {
             const label o = own[facei];
             const label n = nei[facei];
             if (fi[o]*fi[n] < 0)
             {
                 sum[o] += fi[n]; ++cnt[o];
                 sum[n] += fi[o]; ++cnt[n];
             }
         }

         forAll(f.boundaryField(), patchi)
         {
             const fvPatchScalarField& pf = f.boundaryField()[patchi];
             if (pf.coupled())
             {
                 const scalarField fn(pf.patchNeighbourField());
                 const labelUList& fc = pf.patch().faceCells();
                 forAll(fc, i)
                 {
                     const label c = fc[i];
                     if (fi[c]*fn[i] < 0)
                     {
                         sum[c] += fn[i]; ++cnt[c];
                     }
                 }
             }
         }
     };

     scalarField sum(mesh.nCells());
     labelList cnt(mesh.nCells());

     // Ratio of psi0 to the sum of its opposite-sign neighbours: constant
     // during the re-initialisation, computed once per solvePsi()
     // (forcingRatio_ is emptied there)
     if (forcingRatio_.size() != mesh.nCells())
     {
         oppositeSignSum(psi0_, sum, cnt);
         forcingRatio_ = psi0_.primitiveField();
         forAll(forcingRatio_, cid)
         {
             if (cnt[cid] == 0)
                 forcingRatio_[cid] = 0;
             else
                 forcingRatio_[cid] /= sum[cid] + VSMALL;
         }
     }
     const scalarField& rr = forcingRatio_;

     const dimensionedScalar eps = widthEpsilon();

     oppositeSignSum(psi, sum, cnt);
     for (const label cid : narrowBand_.C3())
     {
         if (cnt[cid] == 0)
             FF[cid] = 0.0;
         else
             FF[cid] = (rr[cid]*sum[cid] - psi[cid])/eps.value();
     }

     FF.correctBoundaryConditions();
     
     return tFF;
}

void Foam::interfacePropertiesLS::psi_timeAdvanceRK3
(
    volScalarField& psi,
    const dimensionedScalar& deltaTau,
    const volScalarField& signPsi0,
    const volScalarField& magGradPsi0
)
{
     volScalarField psi1 = psi + deltaTau * psi_timeStepFunc(signPsi0, magGradPsi0);
     volScalarField psi2 = 3.0/4.0 * psi + 1.0/4.0 * psi1 + 1.0/4.0*deltaTau * psi_timeStepFunc(signPsi0, calMagGradPsi(signPsi0, psi1));
     volScalarField psi3 = 1.0/3.0 * psi + 2.0/3.0 * psi2 + 2.0/3.0*deltaTau * psi_timeStepFunc(signPsi0, calMagGradPsi(signPsi0, psi2));

     psi = psi3;
}

void Foam::interfacePropertiesLS::psi_timeAdvanceRK3B
(
    volScalarField& psi,
    const dimensionedScalar& deltaTau
)
{
     volScalarField psi1 = psi + deltaTau * psi_timeStepFunc2(psi);
     volScalarField psi2 = 3.0/4.0 * psi + 1.0/4.0 * psi1 + 1.0/4.0*deltaTau * psi_timeStepFunc2(psi1);
     volScalarField psi3 = 1.0/3.0 * psi + 2.0/3.0 * psi2 + 2.0/3.0*deltaTau * psi_timeStepFunc2(psi2);
     
     psi = psi3;
}

Foam::tmp<volScalarField>
Foam::interfacePropertiesLS::psi_timeStepFunc
(
    const volScalarField& signPsi0,
    const volScalarField& magGradPsi
) const
{
     auto tpsiDt = tmp<volScalarField>::New
     (
         IOobject
         (
             "psiDt",
             psi_.time().timeName(),
             psi_.mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             IOobject::NO_REGISTER
         ),
         psi_.mesh(),
         dimensionedScalar(dimVelocity, Zero),
         zeroGradientFvPatchField<scalar>::typeName
     );
     volScalarField& psiDt = tpsiDt.ref();

     for (const label cid : narrowBand_.C3())
     {
          psiDt[cid] = signPsi0[cid]*(scalar(1) - magGradPsi[cid]);
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

     // Constant during the re-initialisation
     const tmp<volScalarField> tsignPsi0 = calSignFunc();
     const volScalarField& signPsi0 = tsignPsi0();
     forcingRatio_.clear();

     label corr = 0;
     for(corr = 0; corr < numLoop; corr++)
     {
          // |grad(psi)| of the current psi: residual and first RK stage
          const tmp<volScalarField> magGradPsi = calMagGradPsi(signPsi0, psi_);

          resid = calResidualMagGradPsi(magGradPsi);
          if(corr == 0) resid0 = resid;

          if(resid / resid0 <  params_.tolMagGradPsi )
          {
               Info << " -- converged " ;
               break;
          }

          psi_timeAdvanceRK3(psi_, deltaTau, signPsi0, magGradPsi());
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
