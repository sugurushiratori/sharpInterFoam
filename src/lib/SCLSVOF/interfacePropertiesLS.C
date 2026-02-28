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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfacePropertiesLS::correctContactAngle
(
          surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = psi_.mesh();
    const volScalarField::Boundary& abf = psi_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}

//--------------------------------------------
//-  calculate normal vector & curvature
//--------------------------------------------
void Foam::interfacePropertiesLS::calculateK()
{
     const fvMesh& mesh = psi_.mesh();
     const surfaceVectorField& Sf = mesh.Sf();

    volVectorField gradPsi = fvc::grad(psi_, "nHat");;
    if (nAlphaSmoothCurvature_ >= 1)
    {
        // Smooth interface curvature to reduce spurious currents
        auto tpsiL = tmp<volScalarField>::New(psi_);
        auto& psiL = tpsiL.ref();

        for (int i = 0; i < nAlphaSmoothCurvature_; ++i)
        {
            psiL = fvc::average(fvc::interpolate(psiL));
        }

        gradPsi = fvc::grad(tpsiL, "nHat");
    }
    
     surfaceVectorField gradPsif(fvc::interpolate(gradPsi));
     
     // Interpolated face-gradient of psi (dimLess = dimLength / dimLength)
     surfaceVectorField nHatfv(gradPsif / (mag(gradPsif) + deltaN_ ));
     correctContactAngle(nHatfv.boundaryFieldRef(), gradPsif.boundaryField());
          
     // Face unit interface normal flux
     nHatf_ = nHatfv & Sf;

     nHatv_ == gradPsi / (mag(gradPsi) + deltaN_);
     
     // Simple expression for curvature
     K_ = -fvc::div(nHatf_);

}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::interfacePropertiesLS::meshSize() const
{
    // determine length scale from average cell volume
    return dimensionedScalar("meshSize", cbrt(average(alpha1_.mesh().V())) );
}

void Foam::interfacePropertiesLS::calDeltaFunc()
{
    const fvMesh& mesh = psi_.mesh();
    const scalar eps = widthEpsilon().value();
    const scalar pi  = Foam::constant::mathematical::pi;
    const surfaceScalarField psif = fvc::interpolate(psi_);
    
    surfaceScalarField& sf = DeltaFace_;
    scalarField& sif = sf.ref();

    forAll(psif, fi)
    {
        if(mag(psif[fi]) <= eps )
          sif[fi] = scalar(1) / (scalar(2) * eps ) * (scalar(1) + cos(pi * psif[fi] / eps ));
        else   
          sif[fi] = scalar(0);
    }
    
    forAll(mesh.boundary(), patchi)
    {
        fvsPatchField<scalar>& pf = sf.boundaryFieldRef()[patchi];
        forAll(mesh.boundary()[patchi], fi)
        {
            const scalar psiF = psif.boundaryField()[patchi][fi];
            if(mag(psiF) <= eps )
              pf[fi] = scalar(1) / (scalar(2) * eps ) * (scalar(1) + cos(pi * psiF / eps ));
            else   
              pf[fi] = scalar(0);
        }
    }
}

void Foam::interfacePropertiesLS::calHeaviside()
{
    const scalar eps = widthEpsilon().value();
    const scalar pi  = Foam::constant::mathematical::pi;

    volScalarField& vf = Heaviside_;
    scalarField& vif   = vf.ref();

    forAll(psi_, ci)
    {
             if( psi_[ci] < -eps )  vif[ci] = scalar(0);
        else if( psi_[ci] >  eps )  vif[ci] = scalar(1);
        else vif[ci] = scalar(1) / scalar(2) * ( scalar(1) + psi_[ci] / eps + sin(pi * psi_[ci] / eps ) / pi );
    }
    vf.correctBoundaryConditions();
}

void Foam::interfacePropertiesLS::calHeavisideFace()
{
    const fvMesh& mesh = psi_.mesh();
    const scalar eps = widthEpsilon().value();
    const scalar pi  = Foam::constant::mathematical::pi;
    const surfaceScalarField psif = fvc::interpolate(psi_);
    
    surfaceScalarField& sf = HeavisideFace_;
    scalarField& sif = sf.ref();

    forAll(psif, fi)
    {
             if( psif[fi] < -eps )  sif[fi] = scalar(0);
        else if( psif[fi] >  eps )  sif[fi] = scalar(1);
        else sif[fi] = scalar(1) / scalar(2) * ( scalar(1) + psif[fi] / eps + sin(pi * psif[fi] / eps ) / pi );
    }
    
    forAll(mesh.boundary(), patchi)
    {
        fvsPatchField<scalar>& pf = sf.boundaryFieldRef()[patchi];
        forAll(mesh.boundary()[patchi], fi)
        {
           const scalar psiF = psif.boundaryField()[patchi][fi];
                if( psiF < -eps )  pf[fi] = scalar(0);
           else if( psiF >  eps )  pf[fi] = scalar(1);
           else pf[fi] = scalar(1) / scalar(2) * ( scalar(1) + psiF / eps + sin(pi * psiF / eps ) / pi );
        }
    }
}

void Foam::interfacePropertiesLS::calScaledHeavisideFace()
{
    const fvMesh& mesh = psi_.mesh();
    const scalar eps = widthEpsilon().value();
    const scalar pi  = Foam::constant::mathematical::pi;
    const surfaceScalarField psif = fvc::interpolate(psi_);
    
    surfaceScalarField& sf = scaledHeavisideFace_;
    scalarField& sif = sf.ref();

    forAll(psif, fi)
    {
             if( psif[fi] < -eps )  sif[fi] = scalar(0);
        else if( psif[fi] >  eps )  sif[fi] = scalar(1);
        else 
        {
              sif[fi] = scalar(1) / scalar(2) * 
              (
                    pow(scalar(1) + psif[fi] / eps, 2) / scalar(2) 
                   -(cos(scalar(2)*pi*psif[fi] / eps) - scalar(1)) / (scalar(4)*pi*pi)
                   +(eps + psif[fi]) / (eps * pi)*sin(pi*psif[fi] / eps)
              );
        }
    }
    
    forAll(mesh.boundary(), patchi)
    {
        fvsPatchField<scalar>& pf = sf.boundaryFieldRef()[patchi];
        forAll(mesh.boundary()[patchi], fi)
        {
           const scalar psiF = psif.boundaryField()[patchi][fi];
                if( psiF < -eps )  pf[fi] = scalar(0);
           else if( psiF >  eps )  pf[fi] = scalar(1);
           else 
           {
               pf[fi] = scalar(1) / scalar(2) * 
                       (
                             pow(scalar(1) + psiF / eps, 2) / scalar(2) 
                            -(cos(scalar(2)*pi*psiF / eps) - scalar(1)) / (scalar(4)*pi*pi)
                            +(eps + psiF) / (eps * pi)*sin(pi*psiF / eps)
                       );
           }
        }
    }
}


void Foam::interfacePropertiesLS::calScaledDeltaFunc()
{
    scaledDeltaFace_ = scalar(2) * HeavisideFace_ * DeltaFace_;
}

void Foam::interfacePropertiesLS::calScaledHeaviside()
{
    const scalar eps = widthEpsilon().value();
    const scalar pi  = Foam::constant::mathematical::pi;

    volScalarField& vf = scaledHeaviside_;
    scalarField& vif   = vf.ref();

    forAll(psi_, ci)
    {
             if( psi_[ci] < -eps )  vif[ci] = scalar(0);
        else if( psi_[ci] >  eps )  vif[ci] = scalar(1);
        else 
        {
            vif[ci] = scalar(1) / scalar(2) * 
            (
                pow(scalar(1) + psi_[ci] / eps, 2) / scalar(2) 
               -(cos(scalar(2)*pi*psi_[ci] / eps) - scalar(1)) / (scalar(4)*pi*pi)
               +(eps + psi_[ci]) / (eps * pi)*sin(pi*psi_[ci] / eps)
            );
        }
    }
    vf.correctBoundaryConditions();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesLS::surfaceTensionForce() const
{
    if(params_.densityScaled)
    {
       // snGrad(scaledHevaviside) is blanced formulation because "laplacian(p_rgh)" is evaluated by div(snGrad(p_rgh))
       return fvc::interpolate(sigmaK()) * fvc::snGrad( scaledHeaviside_ );       
    }
    else
    {
//       return fvc::interpolate(sigmaK()) * fvc::snGrad(psi_) * DeltaFace_;
       return fvc::interpolate(sigmaK()) * fvc::snGrad( Heaviside_ );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesLS::nearInterface() const
{
    return pos0(alpha1_ - 0.01) * pos0(0.99 - alpha1_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacePropertiesLS::interfacePropertiesLS
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),

    nAlphaSmoothCurvature_
    (
        alpha1.mesh().solverDict(alpha1.name()).
            getOrDefault<int>("nAlphaSmoothCurvature", 0)
    ),
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),
    
    deltaN_
    (
        "deltaN",
        1e-8 / cbrt(average(alpha1.mesh().V())).value()
    ),

    alpha1_(alpha1),
    U_(U),
    
    nHatv_
    (
        IOobject
        (
            "nHatv",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, vector::zero)
    ),


    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),

    K_
    (
        IOobject
        (
            "kappa",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    
    psi_
    (
        IOobject
        (
            "psi",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimLength, Zero),
        alpha1_.boundaryField().types()
    ),
    
    psi0_
    (
        IOobject
        (
            "psi0",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimLength, Zero),
        alpha1_.boundaryField().types()
    ),

    Heaviside_
    (
        IOobject
        (
            "H",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, Zero),
        alpha1_.boundaryField().types()
    ),
    
    HeavisideFace_
    (
        IOobject
        (
            "Hf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, Zero)
    ),

    DeltaFace_
    (
        IOobject
        (
            "DeltaF",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless / dimLength, Zero)
    ),
    
    scaledHeaviside_
    (
        IOobject
        (
            "scaledH",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, Zero),
        alpha1_.boundaryField().types()
    ),
    
    scaledDeltaFace_
    (
        IOobject
        (
            "scaledDelta",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/ dimLength, Zero)
    ),
    
    scaledHeavisideFace_
    (
        IOobject
        (
            "scaledH",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    
    narrowBand_(psi_)
{
	Info << "interfacePropertiesLS: meshSize = " << meshSize().value() << endl;
    read();
    calculateK();
}


bool Foam::interfacePropertiesLS::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);
     
    const dictionary* sclsvofDict =  alpha1_.mesh().solverDict(alpha1_.name()).findDict("SCLSVOF");
   
    params_.correctPsi      = sclsvofDict->lookupOrDefault<Switch>("correctPsi",      params_.correctPsi);
    params_.widthFactor     = sclsvofDict->lookupOrDefault<scalar>("widthFactor",     params_.widthFactor);
    params_.densityScaled   = sclsvofDict->lookupOrDefault<Switch>("densityScaled",   params_.densityScaled);
    params_.initializeAtanh = sclsvofDict->lookupOrDefault<Switch>("initializeAtanh", params_.initializeAtanh);
    params_.denomDeltaTau   = sclsvofDict->lookupOrDefault<scalar>("denomDeltaTau",   params_.denomDeltaTau);
    params_.factorNumLoop   = sclsvofDict->lookupOrDefault<scalar>("factorNumLoop",   params_.factorNumLoop);
    params_.tolMagGradPsi   = sclsvofDict->lookupOrDefault<scalar>("tolMagGradPsi",   params_.tolMagGradPsi);
    params_.densityFunctionHeaviside = sclsvofDict->lookupOrDefault<Switch>("densityFunctionHeaviside",   params_.densityFunctionHeaviside);

    return true;
}

// ************************************************************************* //
