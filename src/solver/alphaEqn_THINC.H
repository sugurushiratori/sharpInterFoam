/**********************************************************
*    Weighted THINC 
*                             Suguru Shiratori
***********************************************************/

{
    const fvMesh& mesh = alpha1.mesh();

    scalarField& alpha1In = alpha1.ref();
    
    if(mesh.moving())
    {
        alpha1In *= (mesh.Vsc0() / mesh.Vsc());
        mixture.correctPsi();
    }
    
    const dimensionedScalar deltaT = alpha1.time().deltaT();
    
    surfaceScalarField Hflux(thinc.getFlux(alpha1) * mesh.magSf());
    Hflux.setOriented();
    
    surfaceScalarField alphaPhi = Hflux / deltaT;

    const volScalarField& alpha0 = alpha1.oldTime();
    alpha1 == alpha0 - fvc::div(Hflux);
    
    alpha1.correctBoundaryConditions();

    rhoPhi = alphaPhi * (rho1 - rho2) + phi * rho2;

    alpha2 = 1.0 - alpha1;
    
    mixture.correct();

//    rhoPhi = mixture.rhof() * phi;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}

