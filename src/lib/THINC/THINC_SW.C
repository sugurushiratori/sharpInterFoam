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

//---------------------------------------------------------
//- Return surfaceScalarField interpolated by THINC
//---------------------------------------------------------
Foam::tmp<Foam::surfaceScalarField>
Foam::THINC::getFlux
(
   const volScalarField& vf
) const
{
    const fvMesh& mesh = mesh_;

    tmp<surfaceScalarField> tsf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "THINC::SW(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("", dimLength, scalar(0))
        )
    );
    surfaceScalarField& sf = tsf.ref();
    
    // Two cells adjacent to the target face
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    
    const surfaceVectorField SfHat = mesh.Sf() / mesh.magSf();
    const surfaceScalarField Uf  = phi_ / mesh.magSf();
    const scalar deltaT = vf.time().deltaT().value();
    
//    const volVectorField& CC  = mesh.C();
//    const surfaceVectorField& Cf  = mesh.Cf();
    
    volScalarField meshV(IOobject("V", mesh.time().timeName(), mesh), mesh, dimensionedScalar("V", dimVolume, scalar(1)));   
    meshV.ref() = mesh.V();
    meshV.correctBoundaryConditions();

    forAll(P, facei)
    {
        const scalar phif = Uf[facei];
        const label  ciup = (phif >= 0.0 ? P[facei] : N[facei] );
        const scalar fP   = vf[P[facei]];
        const scalar fN   = vf[N[facei]];
        const scalar nf   = nHatv_[ciup] & SfHat[facei];
        
        const scalar weight = ::min(scalar(1), mag(nf));    // weight for interface thickness
        const scalar beta = beta0_ * weight;

//        const scalar dx = 2.0 * mag(Cf[facei] - CC[ciup]);

        const scalar dx = meshV[ciup] / mesh.magSf()[facei];

        sf[facei] = flux1D(fP, fN, nf, beta, phif, deltaT, dx);
    }
    
    //- ‹«ŠEˆ—
    forAll(mesh.boundary(), patchi)
    {
       fvsPatchScalarField& psf      = sf.boundaryFieldRef()[patchi];
       const vectorField&   pSfHat   = SfHat.boundaryField()[patchi];
       
       if(vf.boundaryField()[patchi].coupled())
       {
           const fvsPatchScalarField& pUf  = Uf.boundaryField()[patchi];
           const scalarField vfP  = vf.boundaryField()[patchi].patchInternalField();
           const scalarField vfN  = vf.boundaryField()[patchi].patchNeighbourField();
           const vectorField nvP  = nHatv_.boundaryField()[patchi].patchInternalField();
           const vectorField nvN  = nHatv_.boundaryField()[patchi].patchNeighbourField();
           
//           const vectorField cvP  = CC.boundaryField()[patchi].patchInternalField();
//           const vectorField cvN  = CC.boundaryField()[patchi].patchNeighbourField();
           
           const scalarField VP  = meshV.boundaryField()[patchi].patchInternalField();
           const scalarField VN  = meshV.boundaryField()[patchi].patchNeighbourField();

           forAll(mesh.boundary()[patchi], facei)
           {
               const scalar phif  = pUf[facei];
               const vector nv    = (phif >= 0.0 ? nvP[facei]  : nvN[facei]);
//               const vector cV    = (phif >= 0.0 ? cvP[facei]  : cvN[facei]);
               const scalar nf    = nv & pSfHat[facei];
//               const vector ciF   = Cf.boundaryField()[patchi][facei];

               const scalar VV    = (phif >= 0.0 ? VP[facei] : VN[facei]);
               
               const scalar weight = ::min(scalar(1), mag(nf));    // weight for interface thickness
               const scalar beta = beta0_ * weight;
               
//             const scalar dx = 2.0 * mag(cV - ciF);
               const scalar dx = VV / mesh.magSf().boundaryField()[patchi][facei];

               psf[facei] = flux1D(vfP[facei], vfN[facei], nf, beta, phif, deltaT, dx);
           }
       }
       else
       {
           forAll(mesh.boundary()[patchi], facei)
           {
               psf[facei] = vf.boundaryField()[patchi][facei] * Uf.boundaryField()[patchi][facei] * deltaT;
           }
       }
    }
    
    return tsf;
}

// ************************************************************************* //
