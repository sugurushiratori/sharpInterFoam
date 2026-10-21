# Regression tests

Small cases (each below about 2 minutes on 2 cores) that check the main
properties of the solvers. They are run by the CI and should be run
locally before committing changes to the solvers.

```sh
source $WM_PROJECT_DIR/etc/bashrc   # OpenFOAM v2412 or v2512
./Allwmake
tests/regression/run.sh                  # all cases
tests/regression/run.sh zalesak          # selected cases
```

The cases are run in `./regression-work` (set `REGRESSION_WORK` to change).
`metrics.py` evaluates a finished case and `check.py` compares the result
with `reference.json`. Python 3 without extra packages is sufficient.

| Case | Solver | Checks |
|---|---|---|
| zalesak | sharpInterFoam | interface advection error E_alpha, volume conservation, boundedness |
| staticDroplet | sharpInterFoam | spurious velocity and pressure error with the density-scaled CSF |
| translatingDrop | sharpInterFoam | uniform translation at density ratio 1000 (mass-momentum consistency) |
| twoLayerConduction | sharpInterThermoFoam | heat flux through two layers (harmonic conductivity) |
| thermocapillaryCavity | sharpInterThermoFoam | Marangoni flow; serial and 2-process results agree |
| movingMesh | sharpInterFoam | volume conservation on a moving mesh |

The 2D cases with flow in the x-z plane are affected by a defect of some
OpenFOAM builds, see `transposeCheck` in the documentation; the
thermocapillary cavity is therefore set up in the x-y plane.
