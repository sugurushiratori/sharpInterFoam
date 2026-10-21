{
    fvScalarMatrix TEqn
    (
        fvm::ddt(T) + fvm::div(phi, T)
       -fvm::laplacian(mixture.kappa(), T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correct();
}
