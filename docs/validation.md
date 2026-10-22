# Validation results

[日本語](ja/validation.md)

OpenFOAM v2412 (conda-forge, gcc 13), 2 cores. This OpenFOAM build fails
`tests/transposeCheck` (the OpenCFD Debian packages v2412 and v2512 pass);
cases in the x-z plane and axisymmetric cases were therefore run with the
laminar viscous stress computed with an out-of-place transpose, which gives
the result of a correct OpenFOAM build (not part of the solver). Status of the code: v2.0 development
version (see [CHANGELOG](../CHANGELOG.md)).

## Interface advection: Zalesak disc (paper Sec. 4.1)

One revolution, E_alpha = sum|alpha - alpha0| / sum alpha0.

| Mesh | Paper (OpenFOAM v2006) | v2.0 |
|---|---|---|
| structured 100 x 100 | 6.85e-2 | 7.09e-2 |
| unstructured 12785 cells | 8.38e-2 | 8.71e-2 |

The difference of 3-4 % to the paper is due to the OpenFOAM version.
Volume change 5e-10, min(alpha) >= 0.

## Static droplet (paper Sec. 4.2, density-scaled CSF)

Pressure error E_p (paper Eq. 49) and maximum spurious velocity after one step.

| N | U_max | E_p |
|---|---|---|
| 20 | 8.3e-4 | 3.75e-2 |
| 40 | 1.2e-3 | 1.90e-2 |
| 80 | 2.8e-3 | 1.01e-2 |

Axisymmetric (wedge) droplet, Laplace pressure: -0.4 % (N = 20), +0.8 % (N = 40).

## Oscillating drop (axisymmetric, mode l = 2)

`validation/oscillatingDrop`: e = 0.05, density ratio 1000, nu = 0.05,
zeta = Gamma/omega = 0.07. The viscosity is recovered from the measured
damping rate and frequency with the finite-viscosity correction (series
in zeta^(1/2)) and compared with the input value; shape measured by the
quadrupole moment.

| R/h | omega_d / omega_Lamb | nu recovered (correction) | nu from Lamb (R^2 Gamma/5) |
|---|---|---|---|
| 25 | 0.977 | -4.2 % | -16 % |
| 50 | 0.982 | -14 % | -24 % |

**Open issue:** the damping decreases with mesh refinement, also without
the density-scaled model (-14 % at R/h = 50) and with a 4 times smaller time
step at R/h = 25 (no change). The viscous operator itself is correct
(axisymmetric Poiseuille flow within 0.3 %).

## Rising bubble (Hysing et al. 2009, test case 1)

Reference: group TP2D, h = 1/320.

| Quantity | Reference | h = 1/40 | h = 1/80 |
|---|---|---|---|
| minimum circularity | 0.9013 | 0.8955 (-0.6 %) | 0.8960 (-0.6 %) |
| time of minimum circularity | 1.9000 | 1.8881 (-0.6 %) | 1.8810 (-1.0 %) |
| maximum rise velocity | 0.2417 | 0.2406 (-0.5 %) | 0.2412 (-0.2 %) |
| time of maximum rise velocity | 0.9213 | 0.9328 (+1.3 %) | 0.9233 (+0.2 %) |
| centroid height at t = 3 | 1.0813 | 1.0802 (-0.1 %) | 1.0798 (-0.1 %) |

The circularity is computed from the interface length integral of |grad(alpha)|,
which is 0.6 % too long for the initial circle.

## Capillary wave (Prosperetti 1981)

Equal densities and viscosities, La = 3000, initial amplitude 0.01 wavelength.
RMS error of a(t)/a0 against the exact initial-value solution for
0 <= omega0 t <= 25.

| Cells per wavelength | RMS error |
|---|---|
| 32 | 2.7e-2 |
| 64 | 1.2e-2 |

## Thermocapillary migration (Young, Goldstein & Block 1959)

Axisymmetric drop in a linear temperature field, mu' = k' = 1,
Re = Ma = 0.03, Ca = 0.07, container 32 R.

| R/h | drop velocity / V_YGB (t = 3) |
|---|---|
| 12 | 1.11 |
| 24 | 1.11 (t = 2.9) |

With the separate normal/tangential discretisation of the variable surface
tension used before v2.0 the drop moved at -2.0 V_YGB (constant sigma in the
normal force) or 1.6 V_YGB (sigma(T) in the normal force). The remaining
excess of about 10 % grows with the base surface tension and is under
investigation.

## Thermocapillary cavity (paper Sec. 4.3, Table 4)

Wall heights of the interface (normalised by the cavity width) at t = 5,
gas density 1e-3 and gas kinematic viscosity equal to the liquid
(assumed; to be confirmed against the paper's setup).

| Mesh, model | hot wall | cold wall | difference |
|---|---|---|---|
| Paper (400 x 160) | 0.187 | 0.214 | 0.027 |
| Sen & Davis (1982) | 0.188 | 0.213 | 0.025 |
| 50 x 20, density-scaled CSF, Heaviside density | 0.184 | 0.214 | 0.030 |
| 100 x 40, density-scaled CSF, Heaviside density | 0.1823 | 0.2184 | 0.036 |
| 50 x 20, standard CSF, alpha density | 0.1852 | 0.2107 | 0.0255 |

**Open issue:** with the density-scaled CSF model the deformation is larger
than the reference and grows with mesh refinement; with the standard CSF
the deformation agrees with Sen & Davis. Surface velocity at the cavity
centre: 0.1167 (Sen & Davis 0.1169).

## Other checks

| Case | Result |
|---|---|
| Translating drop, density ratio 1000 | max |U - U0| = 3e-15 |
| Two-layer heat conduction, conductivity ratio 10 | heat flux error < 1e-10 |
| Sessile drop, contact angle 60 / 90 deg | 61.3 / 90.3 deg |
| Moving mesh (moving wall) | volume change 2e-15 |
| Serial vs 2 processes (thermocapillary cavity) | U difference 2e-10 |
| Axisymmetric Poiseuille flow (wedge, 20 cells) | velocity profile within 0.3 % |
