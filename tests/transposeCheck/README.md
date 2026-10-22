# transposeCheck

Checks whether `Foam::T()` applied to a temporary tensor field
(`T(fvc::grad(U))`) gives correct results with the installed OpenFOAM
build, and whether the precompiled laminar viscous term
`divDevReff(U) = -div(nuEff*dev2(T(grad(U)))) - laplacian(nuEff, U)` is
correct.

```sh
./Allrun     # prints two lines, each ending with "==> OK" or "==> BROKEN"
./Allclean
```

A coded function object is compiled, so the OpenFOAM development files
(wmake, compiler) are needed.

Results: OpenCFD Debian packages v2412 and v2512: both OK. Conda-forge
OpenFOAM v2412 (gcc 13.4, `-O2 -ftree-vectorize`): both BROKEN, every other cell receives the xz component in
place of zx (and similar). This affects the viscous stress of 3D cases
and of 2D cases in the x-z plane; 2D cases in the x-y plane are not
affected because the xz and zx components vanish.
