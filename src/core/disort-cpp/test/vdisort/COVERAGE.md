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
| 8 | a–c | ported | Both inhomogeneous layers use the shared isotropic-scattering inputs. Direct, diffuse-down, and upward fluxes and all four arbitrary-angle radiances match at the top, exact layer interface, and bottom; Q, U, and V remain zero. This exercises multilayer continuity and formal integration across an interface. |
| 9 | a–c | ported | Six distinct layers reproduce all direct, diffuse-down, and upward fluxes and all 100 stored arbitrary-angle radiances. Coverage progresses from isotropic scattering through an eighth-order phase function to layer-dependent Henyey–Greenstein scattering with beam, three azimuths, top/bottom and internal thermal emission, and a Lambertian surface. Q, U, and V remain zero. |
| 10 | a–b | ported | The four-stream Problem 9c atmosphere is evaluated both through the arbitrary-angle formal solution at the rounded Gauss nodes and directly on the quadrature streams. Every Stokes component agrees, Q, U, and V vanish, and radiance evaluation leaves all combined flux results unchanged. |
| 11 | a–b | ported | The homogeneous layer and its identical three-layer subdivision agree at all four requested depths for every arbitrary-angle radiance and for direct, diffuse-down, upward, and `DFDT` fluxes. Q, U, and V remain zero in both representations. |
| 12 | a–b | ported | The optically thick homogeneous atmosphere and its identical three-layer subdivision use the same scalar delta-M transport without an absorption cutoff. All corrected arbitrary-angle radiances and direct, diffuse-down, upward, and `DFDT` fluxes agree; Q, U, and V remain zero. |
| 13 | a–d | ported | The regular scalar delta-M solves reproduce both single- and two-layer albedo/transmission pairs, covering the physical results of the `IBCND=1` shortcut cases without implementing that shortcut. An absorbing-atmosphere check independently verifies direct transmission and Lambertian beam reflection. |
| 14 | a–d | ported | Hapke, Cox-Munk, RPV, and Ross-Li raw BRDFs are transformed into 32 scalar M00 Fourier modes over the same transparent-layer representation as CPP-DISORT. Every surface-reflected radiance and all direct, diffuse-down, upward, and `DFDT` fluxes match; Q, U, and V remain zero. |
| 15 | a–d | ported | All 600 aerosol moments feed the two-layer scalar reduction of the full cached Mueller IMS/TMS correction, with an exact spectral convolution of the removed peak. Hapke, shadowed Cox-Munk, RPV, and Ross-Li cases reproduce every radiance and all direct, diffuse-down, upward, and `DFDT` fluxes; Q, U, and V remain zero. |
| 16 | a | intentionally unsupported | Pseudo-spherical direct-beam corrections are neither required by polarization nor part of the VDISORT geometry model. This is the only Fortran problem not targeted. |
| 17 | a–b | ported | The validated scalar Gaussian delta-M-plus transformation constructs the M00 VDISORT operator for the aerosol and cloud cases. All 1,080 radiances match and Q, U, and V remain zero. IMS/TMS is deliberately disabled, matching DISORT 4.0.99; this is not a general Mueller-valued delta-M-plus implementation. |

## Geometry scope

VDISORT is a polarized extension of the plane-parallel discrete-ordinate
solver.  Its Stokes and Mueller treatment does not require pseudo-spherical
geometry: polarization changes the transported quantities and scattering
operators, not the assumed geometry.  Both the diffuse field and direct beam
therefore remain plane-parallel.

Problem 16 tests DISORT's pseudo-spherical direct-beam shortcut and is
intentionally unsupported.  We do not intend to add a planet radius,
pseudo-spherical switch, or curved-beam calculation to the VDISORT core.

If ARTS later requires curved direct-beam paths, they should be supplied by
the ARTS geometry and ray-tracing layer, with path-dependent illumination
passed into VDISORT.  This keeps the transport solver independent of a
single-radius planetary approximation and is suitable for ARTS's
multiple-planet use cases.

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
