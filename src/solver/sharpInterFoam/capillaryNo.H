/*---------------------------------------------------------------------------*\
Global
    capillaryNo

Description
    Capillary time-step constraint for the explicit treatment of surface
    tension (Brackbill et al. 1992; Denner & van Wachem 2015):

        deltaT_sigma = sqrt( rhoAvg*h^3/(2*pi*sigma) ),
        rhoAvg = (rho1 + rho2)/2,

    where h is the smallest cell-centre distance across faces next to the
    interface. Reports capillaryNo = deltaT/deltaT_sigma, which is limited
    to maxCapillaryNo (controlDict, default 1) when adjustTimeStep is on.

\*---------------------------------------------------------------------------*/

scalar maxCapillaryNo
(
    runTime.controlDict().getOrDefault<scalar>("maxCapillaryNo", 1)
);

scalar capillaryNo = 0;

{
    const surfaceScalarField nearInterfacef
    (
        fvc::interpolate(mixture.nearInterface())
    );
    const scalarField& dc = mesh.deltaCoeffs().primitiveField();

    scalar hMin = GREAT;
    forAll(dc, facei)
    {
        if (nearInterfacef[facei] > SMALL)
        {
            hMin = min(hMin, 1.0/dc[facei]);
        }
    }
    reduce(hMin, minOp<scalar>());

    const scalar sigmaMax = gMax(mixture.sigma()().primitiveField());

    if (hMin < 0.5*GREAT && sigmaMax > SMALL)
    {
        const scalar rhoAvg = 0.5*(rho1.value() + rho2.value());
        const scalar deltaTSigma =
            Foam::sqrt
            (
                rhoAvg*pow3(hMin)
               /(constant::mathematical::twoPi*sigmaMax)
            );

        capillaryNo = runTime.deltaTValue()/deltaTSigma;
    }
}

Info<< "Capillary time-step ratio (deltaT/deltaT_sigma): " << capillaryNo
    << endl;

if (!adjustTimeStep && capillaryNo > maxCapillaryNo)
{
    WarningInFunction
        << "deltaT exceeds the capillary time-step constraint: "
        << "deltaT/deltaT_sigma = " << capillaryNo
        << " > maxCapillaryNo = " << maxCapillaryNo << nl
        << "    Reduce deltaT or set adjustTimeStep yes." << endl;
}

// ************************************************************************* //
