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
    Foam::interfacePropertiesLS

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfacePropertiesLS.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePropertiesLS_H
#define interfacePropertiesLS_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "narrowBand.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfacePropertiesLS
{
public:
        struct Parameters
        {
            bool   correctPsi;      //- switch whether calculate Level-Set function. otherwise set psi = alpha.
            scalar widthFactor;     //- parameter controls the interface width
            bool   densityScaled;   //- switch for density-scaled CSF model
            bool   initializeAtanh; //- true if initialize psi with the function atanh
            scalar denomDeltaTau; 
            scalar factorNumLoop;
            scalar tolMagGradPsi;
            
            bool   densityFunctionHeaviside;
            
            Parameters()
              : 
                correctPsi(true),
                widthFactor(1.5),
                densityScaled(true),
                initializeAtanh(true),
                denomDeltaTau(10.0),
                factorNumLoop(1.0),
                tolMagGradPsi(0.1),
                densityFunctionHeaviside(false)
            {}
            
            void operator=(const Parameters& p)
            {
                correctPsi      = p.correctPsi;
                widthFactor     = p.widthFactor;
                densityScaled   = p.densityScaled;
                initializeAtanh = p.initializeAtanh;
                denomDeltaTau   = p.denomDeltaTau;
                factorNumLoop   = p.factorNumLoop;
                tolMagGradPsi   = p.tolMagGradPsi;
                densityFunctionHeaviside  = p.densityFunctionHeaviside;
            }
        };

private:
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Number of iterations to smooth interfacial curvatures
        int nAlphaSmoothCurvature_;
        
        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& alpha1_;
        const volVectorField& U_;
        volVectorField     nHatv_;
        surfaceScalarField nHatf_;
        volScalarField K_;
        
        Parameters params_;

        
        //- Level-Set function
        volScalarField psi_;
        volScalarField psi0_;
        
        //- Heaviside functions
        volScalarField     Heaviside_;
        surfaceScalarField HeavisideFace_;
        
        surfaceScalarField DeltaFace_;
        
        //- scaled functions for density-scaled CSF
        volScalarField     scaledHeaviside_;
        surfaceScalarField scaledDeltaFace_;
        
        surfaceScalarField scaledHeavisideFace_;
        
        narrowBand narrowBand_;

        
    // Private Member Functions
        //- No copy construct
        interfacePropertiesLS(const interfacePropertiesLS&) = delete;

        //- No copy assignment
        void operator=(const interfacePropertiesLS&) = delete;

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
                  surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf
        ) const;

        
        //- initialize Level-Set function
        void initializePsi();

        //- Solve Level-Set re-initialization equation
        void solvePsi();
        
        tmp<volScalarField> calSignFunc(const volScalarField& psi);
        
        void psi_timeAdvanceRK3(volScalarField& psi, const dimensionedScalar& deltaTau);
        void psi_timeAdvanceRK3B(volScalarField& psi, const dimensionedScalar& deltaTau);
        
        tmp<volScalarField> psi_timeStepFunc(const volScalarField& psi);
        tmp<volScalarField> psi_timeStepFunc2(const volScalarField& psi);
        
        void K_timeAdvanceRK3(volScalarField& K, const dimensionedScalar& deltaTau);
        tmp<volScalarField> K_timeStepFunc(const volScalarField& K);
        
        //- calculate gradPsi with Godunov Hamiltonian scheme
        tmp<volVectorField> calGradPsi
        (
           const volScalarField& signPsi0, 
           const volScalarField& psi
        ) const;
        
        //- calculate magGradPsi with Godunov Hamiltonian scheme
        tmp<volScalarField> calMagGradPsi
        (
           const volScalarField& signPsi0, 
           const volScalarField& psi
        ) const;
        
        tmp<volScalarField> forcingTerm(const volScalarField& psi);

        // 
        scalar calResidualMagGradPsi
        (
           const volScalarField& magGradPsi
        ) const;

        //- interpolate the curvature for density-scaled CSF
        void curvatureInterpolation();

        //- Re-calculate the interface curvature
        void calculateK();
        
        void calHeaviside();
        void calHeavisideFace();
        void calDeltaFunc();
        void calScaledHeaviside();
        void calScaledDeltaFunc();
        void calScaledHeavisideFace();
        
public:


    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfacePropertiesLS
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );
        
        
    // Destructor
        virtual ~interfacePropertiesLS() = default;


    // Member Functions

        dimensionedScalar meshSize() const;

        dimensionedScalar widthEpsilon() const
        {
            return params_.widthFactor * meshSize();
        }
        
        scalar widthFactor() const
        {
            return params_.widthFactor;
        }

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }
        
        bool densityFunctionHeaviside() const
        {
            return params_.densityFunctionHeaviside;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }
        
        const volVectorField& nHatv() const
        {
            return nHatv_;
        }

        const volScalarField& K() const
        {
            return K_;
        }
        
        const volScalarField& psi() const
        {
            return psi_;
        }

        tmp<volScalarField> sigmaK() const
        {
            return sigmaPtr_->sigma()*K_;
        }
        
        
        //- Heaviside function
        const volScalarField& Heaviside() const
        {
            return Heaviside_;
        }

        //- Heaviside function
        const surfaceScalarField& HeavisideFace() const
        {
            return HeavisideFace_;
        }
        
        //- smoothed delta function
        const surfaceScalarField& DeltaFace() const
        {
            return DeltaFace_;
        }
        
        //- scaled delta function
        const surfaceScalarField& scaledDelta() const
        {
            return scaledDeltaFace_;
        }
        
        //- scaled Heaviside function
        const volScalarField& scaledHeaviside() const
        {
            return scaledHeaviside_;
        }
        
        //- scaled Heaviside function
        const surfaceScalarField& scaledHeavisideFace() const
        {
            return scaledHeavisideFace_;
        }
        
        
        tmp<surfaceScalarField> surfaceTensionForce() const;
                
        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- correct Level-Set and Heaviside functions
        void correctPsi()
        {
            initializePsi();
            
            narrowBand_.update(widthEpsilon().value());
            
            if(params_.correctPsi)
            {
                solvePsi();
            }
            calculateK();
            
            calHeaviside();
            calHeavisideFace();
            calDeltaFunc();

            if(params_.densityScaled)
            {
                calScaledHeaviside();
                calScaledHeavisideFace();
                calScaledDeltaFunc();
            }
            else
            {
//                calDeltaFunc();
            }
        }
        
        //- calculate curvature
        void correctCurvature()
        {
            if(params_.densityScaled)
            {
                curvatureInterpolation();
            }
        }

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
