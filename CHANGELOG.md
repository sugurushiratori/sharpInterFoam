# Changelog

## v2.0 (2026-10-22)

Changes from the version used for the paper (tag `v1.0-ijmsa2021`).
The numbers given are from the validation runs recorded in
[docs/validation.md](docs/validation.md).

### New
- `sharpInterThermoFoam`: energy equation (rho*Cp advective form, harmonic
  face conductivity, optional Cp) and Marangoni force; previously only in
  the laboratory version.
- Temperature-dependent surface tension `sigma(T) = sigma + sigmaTemp*(T - TRef)`:
  the deviation from sigma at the mean interface temperature is applied as
  the divergence of a surface stress (normal and tangential part together),
  so that the net force over a closed interface vanishes. The paper version
  used a constant sigma in the normal force; a drop in a temperature
  gradient then moved towards the cold side (Young-Goldstein-Block
  migration: -2.0 V_YGB before, 1.1 V_YGB after).
- Capillary time-step constraint (`maxCapillaryNo` in controlDict).
- Moving meshes: volume-conserving (ALE) THINC update.
- Contact-angle boundary conditions of alpha are applied to the level set.
- `Allwmake` / `Allwclean`; regression tests (`tests/regression`), CI for
  OpenFOAM v2412 and v2512, validation cases (`validation`), Japanese README.

### Fixed
- Curvature for the density-scaled CSF model: the interProcInterpolation-based
  extension (not usable with recent OpenFOAM versions) is replaced by a
  smoothing in the narrow band; no projection to the interface (it made the
  thermocapillary-cavity interface relax back to flat).
- Mass-momentum consistency with the Heaviside-based mixture density
  (translating drop at density ratio 1000: velocity error 4.4 -> 3e-15).
- THINC: fluxes from the old-time alpha in every PIMPLE outer iteration;
  face fluxes evaluated in the upwind-cell orientation; flux continuous for
  small steepness (results no longer depend on round-off).
- Parallel runs give the same results as serial runs (re-initialisation
  forcing term and curvature across processor boundaries).
- Cell size for 2D and axisymmetric meshes from the empty/wedge face area.
- Initial density from the mixture; duplicate field registration; duplicate
  calls of `mixture.correct()`.
- Tutorials: time step within the capillary limit, Euler time scheme,
  executable scripts.

### Performance
- Level-set re-initialisation: Godunov differences in one loop over the
  faces of the narrow band, constant terms computed once per
  re-initialisation (about 2x faster for the thermocapillary cavity).

### Other
- `test_case` renamed to `tutorials`; sources in UTF-8 with English comments.
