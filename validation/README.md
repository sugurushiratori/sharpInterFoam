# Validation cases

Benchmark cases with evaluation scripts. Each case is run by
`./Allrun [resolution]` and evaluated by `evaluate.py` (called at the end of
`Allrun`); `./Allclean` restores the case. Results are collected in
[docs/validation.md](../docs/validation.md).

| Case | Solver | Reference | Resolution argument |
|---|---|---|---|
| [risingBubble](risingBubble) | sharpInterFoam | Hysing et al., Int. J. Numer. Meth. Fluids 60 (2009) 1259, test case 1 | cells per unit length (default 40) |
| [capillaryWave](capillaryWave) | sharpInterFoam | Prosperetti, Phys. Fluids 24 (1981) 1217; Popinet, J. Comput. Phys. 228 (2009) 5838 | cells per wavelength (default 32) |
| [thermocapillaryMigration](thermocapillaryMigration) | sharpInterThermoFoam | Young, Goldstein & Block, J. Fluid Mech. 6 (1959) 350 | cells per drop radius (default 12) |

`capillaryWave/evaluate.py` needs numpy and scipy; the other scripts need
Python 3 only.

The thermocapillary migration case is axisymmetric (wedge). On OpenFOAM
builds that fail `tests/transposeCheck` the laminar viscous stress of this
case is wrong; check the installation first.
