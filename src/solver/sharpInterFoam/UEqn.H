    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    // With the Heaviside-based density (densityFunctionHeaviside yes) rho is
    // not transported by rhoPhi, so ddt(rho) + div(rhoPhi) != 0. Remove this
    // continuity error from the momentum equation (non-conservative,
    // advective form) to keep momentum transport consistent; otherwise a
    // drop translating at uniform velocity is strongly disturbed at large
    // density ratios.
    if (mixture.densityFunctionHeaviside())
    {
        UEqn -= fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U);
    }

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  + mixture.marangoniForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
