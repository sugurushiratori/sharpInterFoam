{
    // Energy equation in advective form,
    //     rho*Cp*(dT/dt + U & grad(T)) = div(lambda*grad(T)),
    // discretised as the conservative form minus the continuity error of
    // rho*Cp (consistent for both the alpha- and Heaviside-based mixture
    // properties). The face conductivity is the harmonic mean, which
    // ensures continuity of the heat flux across the interface.
    // Without "Cp" in the phase dictionaries rho*Cp = 1 and lambda = kappa.

    rhoCp == mixture.rhoCp();

    const surfaceScalarField rhoCpPhi(fvc::interpolate(rhoCp)*phi);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T, "div(phi,T)")
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
      - fvm::laplacian(mixture.lambdaf(), T)
     ==
        fvOptions(rhoCp, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);
}
