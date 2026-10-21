/**********************************************************
*    Weighted THINC 
*                             Suguru Shiratori
***********************************************************/

{
    const fvMesh& mesh = alpha1.mesh();

    // Explicit THINC step from the old-time value of the (sub-)step: in PIMPLE
    // outer iterations other than the first, alpha1 holds the result of the
    // previous outer iteration and must not be used to compute the fluxes.
    alpha1 == alpha1.oldTime();
    alpha1.correctBoundaryConditions();

    scalarField& alpha1In = alpha1.ref();
    
    if(mesh.moving())
    {
        alpha1In *= (mesh.Vsc0() / mesh.Vsc());
        mixture.correctPsi();
        nHatTHINC = mixture.nHatv();
    }
    
    const dimensionedScalar deltaT = alpha1.time().deltaT();
    
    surfaceScalarField Hflux(thinc.getFlux(alpha1) * mesh.magSf());
    Hflux.setOriented();
    
    surfaceScalarField alphaPhi = Hflux / deltaT;

    // alpha1 holds the old-time value here, scaled by Vsc0/Vsc on moving
    // meshes (ALE: alpha^{n+1} V^{n+1} = alpha^n V^n - sum F)
    alpha1.ref() = alpha1() - fvc::div(Hflux)().internalField();
    
    alpha1.correctBoundaryConditions();

    rhoPhi = alphaPhi * (rho1 - rho2) + phi * rho2;

    alpha2 = 1.0 - alpha1;
    
    mixture.correct();

    // Normal of the updated alpha for the next sub-cycle
    nHatTHINC = mixture.nHatv();

//    rhoPhi = mixture.rhof() * phi;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}

