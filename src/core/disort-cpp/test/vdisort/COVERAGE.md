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
| 2 | a–d | not ported | Reuse the Problem 1 adapter. Exact-unity cases must exercise conservative scattering without changing the physical input. |
| 3 | a–b | not ported | Add scalar-limit IMS/TMS corrections at the original user angles and compare every corrected radiance. This does not require a general polarized correction yet. |
| 4 | a–c | not ported | Reuse arbitrary-angle scalar IMS/TMS for the Haze-L phase function and compare all fluxes and user-angle radiances. |
| 5 | a–b | not ported | Reuse the shared Cloud C.1 data and scalar IMS/TMS path. Compare the complete original output set. |
| 6 | a–h | not ported | Reuse the extracted references and the CPP-DISORT thermal/boundary normalization. Port all active flux assertions before adding any unavailable radiance references. |
| 7 | a–e | not ported | Reuse the extracted references and established thermal/boundary mapping; compare all active flux quantities. |
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
