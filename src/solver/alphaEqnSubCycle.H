const scalar beta0 = Foam::constant::mathematical::pi / mixture.widthFactor();
const dictionary* thincDict = mesh.solverDict(alpha1.name()).findDict("THINC");
if(thincDict == NULL)
      FatalError << "dictionary THINC is not found." <<nl;

THINC  thinc(mesh, phi, mixture.psi(), mixture.nHatv(), *thincDict, beta0);


if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );
    
    tmp<volScalarField> trSubDeltaT;

    if(LTS)
    {
          trSubDeltaT = fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
     }


    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        if(thinc.enabled())
        {
            #include "alphaEqn_THINC.H"
        }
        else
        {
            #include "alphaEqn_MULES.H"
        }
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
        if(thinc.enabled())
        {
            #include "alphaEqn_THINC.H"
        }
        else
        {
            #include "alphaEqn_MULES.H"
        }
}

//rho == alpha1*rho1 + alpha2*rho2;
rho == mixture.rho();
