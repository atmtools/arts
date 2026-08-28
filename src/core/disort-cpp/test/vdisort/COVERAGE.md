# VDISORT coverage of the Fortran DISORT reference tests

`DISOTESTAUX.f` is the canonical numerical reference.  The target is for
VDISORT, under a strict scalar embedding, to pass the same 16 supported
problems as CPP-DISORT: Problems 1–15 and 17.  Problem 16 remains intentionally
unsupported because both solvers are kept plane-parallel.

A problem is only marked `ported` when VDISORT:

- uses the shared Fortran input and reference data;
- embeds scalar radiation as `[I, 0, 0, 0]` and scalar scattering in the M00
  Mueller component;
- reproduces every output currently asserted by the CPP-DISORT port; and
- verifies that Q, U, and V remain zero to numerical tolerance.

Agreement with CPP-DISORT for a simplified atmosphere, finite-output checks,
or the existing synthetic scalar-limit tests does not count as a port.

| Problem | Cases | Current state | Work required for the VDISORT port |
|---|---|---|---|
| 1 | a–f | ported | All user-angle radiances and direct, diffuse-down, and upward fluxes match the shared Fortran references. Exact conservative scattering is exercised by 1b and 1e, and Q, U, and V are asserted to vanish. |
| 2 | a–d | ported | The shared scalar Legendre-to-VDISORT Fourier embedding reproduces every Rayleigh user-angle radiance and flux. Cases 2b and 2d retain exact `omega = 1`, and Q, U, and V are asserted to vanish. |
| 3 | a–b | ported | Classical scalar delta-M preprocessing, physical/scaled optical-depth mapping, and arbitrary-angle TMS/IMS reproduce every corrected Henyey–Greenstein radiance and flux. Corrections affect only I; Q, U, and V are asserted to vanish. |
| 4 | a–c | ported | Haze-L delta-M transport and arbitrary-angle TMS/IMS reproduce every radiance at all requested depths and azimuths, together with all fluxes. Q, U, and V remain zero. |
| 5 | a–b | ported | All Cloud C.1 moments feed the 48-stream scalar delta-M transport and TMS/IMS correction path. Every original radiance and flux matches, including the strong forward aureole, and Q, U, and V remain zero. |
| 6 | a–h | ported | All active direct, diffuse-down, upward, and `DFDT` flux references are asserted. The port covers the transparent limit, absorption, Lambertian and Hapke reflection, top/bottom emission, and linear atmospheric thermal sources; Q, U, and V remain zero. Cases 6f–h use a 1% tolerance for the original single-precision non-Lambertian thermal-boundary integration, while all other cases use 0.1%. No active radiance references exist for this problem. |
| 7 | a–e | ported | All active direct, diffuse-down, upward, and `DFDT` fluxes are asserted for scattering plus linear atmospheric emission, beam and isotropic illumination, and black, perfectly reflecting Lambertian, and Hapke surfaces. Tiny case-7b fluxes use combined relative and absolute tolerances; the legacy Hapke case 7e uses 1%, and Q, U, and V remain zero. No active radiance references exist for this problem. |
| 8 | a–c | not ported | Reuse the shared two-layer inputs and references. Compare all fluxes and four user-angle radiances. |
| 9 | a–c | not ported | Port the multilayer thermal, beam, and Lambertian cases. Centralize scalar BRDF-to-M00 embedding for reuse by Problems 14 and 15. |
| 10 | a–b | not ported | Add a test-facing arbitrary-angle scalar formal-solution path and compare it with quadrature output. No public ARTS user-angle interface is required for this port. |
| 11 | a–b | not ported | Compare the homogeneous layer with its three-layer subdivision for radiances and fluxes. Add VDISORT `DFDT` so the complete reference output can be asserted. |
| 12 | a–b | not ported | Compare the one-layer and equivalent three-layer solutions without an absorption cutoff, including radiances, fluxes, and `DFDT`. |
| 13 | a–d | not ported | Reproduce both regular-solve albedo/transmission pairs. `IBCND=1` itself is unnecessary, as in CPP-DISORT. |
| 14 | a–d | not ported | Embed Hapke, Cox-Munk, RPV, and Ross-Li scalar Fourier modes in M00 and compare all radiances, fluxes, and `DFDT`. Use the same transparent-layer representation as CPP-DISORT. |
| 15 | a–d | not ported | Reuse all 600 aerosol moments, the scalar BRDF embedding, and scalar IMS/TMS corrections. Compare every original radiance, flux, and `DFDT` value. |
| 16 | a | intentionally unsupported | Pseudo-spherical direct-beam corrections are out of scope. This is the only Fortran problem not targeted. |
| 17 | a–b | not ported | For the scalar reduction, apply the validated scalar Gaussian delta-M-plus transformation before constructing the M00 VDISORT operator. Compare all 1,080 radiances; a general Mueller-valued delta-M formulation is not required for this scalar test. |

## Common scalar-reduction infrastructure

The ports must share one adapter rather than adding problem-specific
normalization code.  It should provide:

- directional M00 phase and beam matrices constructed from scalar Legendre
  moments with the exact VDISORT Fourier normalization;
- scalar boundary conditions, polynomial sources, direct beam, and BRDF modes
  embedded in the combined Stokes representation;
- VDISORT radiance and flux extraction with Q/U/V-zero assertions;
- the test-facing arbitrary-angle formal solution needed by the Fortran
  references;
- scalar-limit IMS/TMS corrections using the already validated CPP-DISORT
  conventions;
- `DFDT` from the combined flux evaluation; and
- scalar delta-M-plus preprocessing for Problem 17.

## Porting order

The intended order is:

1. Problem 1, establishing the shared adapter and normalization;
2. Problem 2, establishing exact conservative scattering;
3. Problems 8, 9, 11, 12, and 13, extending multilayer, source, BRDF, flux,
   and `DFDT` coverage;
4. Problems 3–5, adding scalar IMS/TMS and test-facing user angles;
5. Problems 6 and 7, completing thermal and boundary cases;
6. Problems 14 and 15, completing physical BRDF and aerosol coverage; and
7. Problem 17, completing scalar Gaussian delta-M-plus coverage.

Problem 10 should be enabled as soon as the test-facing arbitrary-angle path
exists.  Problem 16 will remain the sole intentionally unsupported problem.
