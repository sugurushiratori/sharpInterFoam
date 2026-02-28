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

#include "fvCFD.H"
#include "gaussConvectionScheme.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "upwind.H"
#include "downwind.H"
#include "interfacePropertiesLS.H"



Foam::tmp<Foam::volVectorField>
Foam::interfacePropertiesLS::calGradPsi(const volScalarField& signPsi0, const volScalarField& psi) const
{
    // grad(psi)の計算: Godunov's scheme

     const fvMesh& mesh  = psi_.mesh();
     
     tmp<volVectorField> tgradPsi
     (
        new volVectorField
        (
          IOobject
          (
              "gradPsi",
              psi.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedVector("gradPsi", psi.dimensions() / dimLength, vector::zero),
          zeroGradientFvPatchField<vector>::typeName
        )
     );     
     volVectorField& gradPsi = tgradPsi.ref();
     
     const surfaceScalarField exf = dimensionedVector("ex", dimless, vector(1, 0, 0)) & mesh.Sf();
     const surfaceScalarField eyf = dimensionedVector("ey", dimless, vector(0, 1, 0)) & mesh.Sf();
     const surfaceScalarField ezf = dimensionedVector("ez", dimless, vector(0, 0, 1)) & mesh.Sf();
     
     upwind<scalar>   intpUx(mesh, exf),  intpUy(mesh, eyf),  intpUz(mesh, ezf);
     downwind<scalar> intpDx(mesh, exf),  intpDy(mesh, eyf),  intpDz(mesh, ezf);
     
     surfaceScalarField psifUx = intpUx.interpolate(psi);
     surfaceScalarField psifUy = intpUy.interpolate(psi);
     surfaceScalarField psifUz = intpUz.interpolate(psi);
     surfaceScalarField psifDx = intpDx.interpolate(psi);
     surfaceScalarField psifDy = intpDy.interpolate(psi);
     surfaceScalarField psifDz = intpDz.interpolate(psi);

     const volScalarField GxU = fvc::div(exf*psifUx);
     const volScalarField GyU = fvc::div(eyf*psifUy);
     const volScalarField GzU = fvc::div(ezf*psifUz);
     const volScalarField GxD = fvc::div(exf*psifDx);
     const volScalarField GyD = fvc::div(eyf*psifDy);
     const volScalarField GzD = fvc::div(ezf*psifDz);

     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label ci = *iter;
         const scalar posGxU = max(scalar(0), GxU[ci]);
         const scalar posGyU = max(scalar(0), GyU[ci]);
         const scalar posGzU = max(scalar(0), GzU[ci]);
         const scalar negGxU = min(scalar(0), GxU[ci]);
         const scalar negGyU = min(scalar(0), GyU[ci]);
         const scalar negGzU = min(scalar(0), GzU[ci]);
         
         const scalar posGxD = max(scalar(0), GxD[ci]);
         const scalar posGyD = max(scalar(0), GyD[ci]);
         const scalar posGzD = max(scalar(0), GzD[ci]);
         const scalar negGxD = min(scalar(0), GxD[ci]);
         const scalar negGyD = min(scalar(0), GyD[ci]);
         const scalar negGzD = min(scalar(0), GzD[ci]);

         if(signPsi0[ci] > 0.0)
         {
              gradPsi[ci][0] = (mag(posGxU) > mag(negGxD) ? posGxU : negGxD );
              gradPsi[ci][1] = (mag(posGyU) > mag(negGyD) ? posGyU : negGyD );
              gradPsi[ci][2] = (mag(posGzU) > mag(negGzD) ? posGzU : negGzD );
         }
         else
         {
              gradPsi[ci][0] = (mag(posGxD) > mag(negGxU) ? posGxD : negGxU );
              gradPsi[ci][1] = (mag(posGyD) > mag(negGyU) ? posGyD : negGyU );
              gradPsi[ci][2] = (mag(posGzD) > mag(negGzU) ? posGzD : negGzU );
         }
     }

     gradPsi.correctBoundaryConditions();

     return tgradPsi;
}

Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesLS::calMagGradPsi(const volScalarField& signPsi0, const volScalarField& psi) const
{
    // |grad(psi)|の計算: Godunov's scheme

     const fvMesh& mesh  = psi_.mesh();
     
     tmp<volScalarField> tmagGradPsi
     (
        new volScalarField
        (
          IOobject
          (
              "magGradPsi",
              psi.time().timeName(),
              mesh
          ),
          mesh,
          dimensionedScalar("magGradPsi", psi.dimensions() / dimLength, scalar(0)),
          zeroGradientFvPatchField<scalar>::typeName
        )
     );     
     volScalarField& magGradPsi = tmagGradPsi.ref();
     
     const surfaceScalarField exf = dimensionedVector("ex", dimless, vector(1, 0, 0)) & mesh.Sf();
     const surfaceScalarField eyf = dimensionedVector("ey", dimless, vector(0, 1, 0)) & mesh.Sf();
     const surfaceScalarField ezf = dimensionedVector("ez", dimless, vector(0, 0, 1)) & mesh.Sf();
     
     upwind<scalar>   intpUx(mesh, exf),  intpUy(mesh, eyf),  intpUz(mesh, ezf);
     downwind<scalar> intpDx(mesh, exf),  intpDy(mesh, eyf),  intpDz(mesh, ezf);
     
     surfaceScalarField psifUx = intpUx.interpolate(psi);
     surfaceScalarField psifUy = intpUy.interpolate(psi);
     surfaceScalarField psifUz = intpUz.interpolate(psi);
     surfaceScalarField psifDx = intpDx.interpolate(psi);
     surfaceScalarField psifDy = intpDy.interpolate(psi);
     surfaceScalarField psifDz = intpDz.interpolate(psi);

     const volScalarField GxU = fvc::div(exf*psifUx);
     const volScalarField GyU = fvc::div(eyf*psifUy);
     const volScalarField GzU = fvc::div(ezf*psifUz);
     const volScalarField GxD = fvc::div(exf*psifDx);
     const volScalarField GyD = fvc::div(eyf*psifDy);
     const volScalarField GzD = fvc::div(ezf*psifDz);

     forAllConstIter(labelList, narrowBand_.C3(), iter)
     {
         const label ci = *iter;
         const scalar posGxU = max(scalar(0), GxU[ci]);
         const scalar posGyU = max(scalar(0), GyU[ci]);
         const scalar posGzU = max(scalar(0), GzU[ci]);
         const scalar negGxU = min(scalar(0), GxU[ci]);
         const scalar negGyU = min(scalar(0), GyU[ci]);
         const scalar negGzU = min(scalar(0), GzU[ci]);
         
         const scalar posGxD = max(scalar(0), GxD[ci]);
         const scalar posGyD = max(scalar(0), GyD[ci]);
         const scalar posGzD = max(scalar(0), GzD[ci]);
         const scalar negGxD = min(scalar(0), GxD[ci]);
         const scalar negGyD = min(scalar(0), GyD[ci]);
         const scalar negGzD = min(scalar(0), GzD[ci]);

         if(signPsi0[ci] > 0.0)
         {
              const scalar vx = max(sqr(posGxU), sqr(negGxD));
              const scalar vy = max(sqr(posGyU), sqr(negGyD));
              const scalar vz = max(sqr(posGzU), sqr(negGzD));
              
              magGradPsi[ci] = sqrt(vx + vy + vz);
         }
         else
         {
              const scalar vx = max(sqr(negGxU), sqr(posGxD));
              const scalar vy = max(sqr(negGyU), sqr(posGyD));
              const scalar vz = max(sqr(negGzU), sqr(posGzD));
              
              magGradPsi[ci] = sqrt(vx + vy + vz);
         }
     }

     magGradPsi.correctBoundaryConditions();

     return tmagGradPsi;
}



// ************************************************************************* //
