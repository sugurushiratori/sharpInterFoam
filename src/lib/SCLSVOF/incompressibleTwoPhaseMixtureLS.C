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

#include "incompressibleTwoPhaseMixtureLS.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseMixtureLS, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::incompressibleTwoPhaseMixtureLS::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu() / (densityFunc_ * rho1_ + (scalar(1) - densityFunc_) * rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseMixtureLS::incompressibleTwoPhaseMixtureLS
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),
    
    kappa1_("kappa", dimViscosity, nuModel1_->viscosityProperties().getOrDefault<scalar>("kappa", 0.0)),
    kappa2_("kappa", dimViscosity, nuModel2_->viscosityProperties().getOrDefault<scalar>("kappa", 0.0)),
    
    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimViscosity, Zero),
        calculatedFvPatchScalarField::typeName
    ),

    densityFunc_
    (
        IOobject
        (
            "densityFunc",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    ),
    
    densityFuncFace_
    (
        IOobject
        (
            "densityFuncF",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    )

{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseMixtureLS::rho() const
{
    return tmp<volScalarField>::New
    (
        "rho",
         densityFunc_ * rho1_ + (scalar(1) - densityFunc_) * rho2_
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseMixtureLS::rhof() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "rhof",
            densityFuncFace_ * rho1_ + (scalar(1) - densityFuncFace_) * rho2_
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseMixtureLS::kappa() const
{
    return tmp<volScalarField>::New
    (
       "kappa",
       densityFunc_ * kappa1_ 
       + (scalar(1) - densityFunc_) * kappa2_ 
    );
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseMixtureLS::mu() const
{
    return tmp<volScalarField>::New
    (
       "mu",
       densityFunc_ * rho1_*nuModel1_->nu()
       + (scalar(1) - densityFunc_) * rho2_*nuModel2_->nu()
    );
}

Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseMixtureLS::mu(const label patchI) const
{
	return mu()().boundaryField()[patchI];
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseMixtureLS::muf() const
{
    return tmp<surfaceScalarField>::New
    (
        "muf",
         densityFuncFace_ * rho1_*fvc::interpolate(nuModel1_->nu())
        + (scalar(1) - densityFuncFace_ ) * rho2_*fvc::interpolate(nuModel2_->nu())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseMixtureLS::nuf() const
{
    return tmp<surfaceScalarField>::New
    (
        "nuf",
        (
            densityFuncFace_               * rho1_*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - densityFuncFace_) * rho2_*fvc::interpolate(nuModel2_->nu())
        ) / ( densityFuncFace_             * rho1_ 
          + (scalar(1) - densityFuncFace_) * rho2_)
    );
}


bool Foam::incompressibleTwoPhaseMixtureLS::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
            nuModel1_->viscosityProperties().readEntry("rho", rho1_);
            nuModel2_->viscosityProperties().readEntry("rho", rho2_);

            if(nuModel1_->viscosityProperties().found("kappa"))
            {
                 nuModel1_->viscosityProperties().readEntry("kappa", kappa1_);
            }
            if(nuModel2_->viscosityProperties().found("kappa"))
            {
                 nuModel2_->viscosityProperties().readEntry("kappa", kappa2_);
            }

            return true;
        }
    }

    return false;
}


// ************************************************************************* //
