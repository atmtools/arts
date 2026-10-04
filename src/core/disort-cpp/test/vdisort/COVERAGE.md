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

## Polarized comparison with RT4

`cpp.fast.vdisort-rt4-test` (`vdisort-rt4-comparison.cpp`) is built only with
`ENABLE_RT4=ON`.  It compares VDISORT with Evans' RT4, a polarized
doubling-adding solver (`src/core/rt4`, `doc/arts/dev.rt4.rst`).  This is the
external reference that the scalar ports above cannot provide.  Both solvers
get the same discrete problem:

- RT4's double-Gauss streams, which are VDISORT's (asserted);
- the same azimuthally averaged phase matrix on those streams;
- the same scalar extinction, absorption and Planck source, linear in optical
  depth;
- the same sky and surface.

They then solve the same linear system of ODEs in optical depth.  VDISORT
solves it to round-off.  RT4 solves it with an error that is first order in
the thickness of its initial doubling layer.  So the difference must vanish
with `max_delta_tau`, and it does:

- it halves exactly when `max_delta_tau` halves, and falls by 8 per decade;
- at `max_delta_tau = 1e-7` it is 0.1 to 0.7 times the initial-layer
  thickness, at most 4e-8 relative to max I;
- against RT4 Richardson-extrapolated to `max_delta_tau = 0`, it is at most
  3.4e-10, and is limited by RT4's second-order error and round-off;
- gas-only atmospheres, which RT4 integrates analytically, agree to 5e-15.

Every level (top, each interface, bottom), every stream and both directions
are compared, for I and Q.

**What it validates.** Polarized m = 0 I/Q multiple scattering with thermal
sources:

- Rayleigh, including thick (tau = 20) and exactly conservative layers;
- a forward-peaked polarizing phase matrix (polarized Henyey-Greenstein,
  g = 0.7), from a numerical azimuthal average that is first validated
  against the Rayleigh closed form;
- a constructed non-reciprocal, Stokes-asymmetric phase matrix;
- multilayer atmospheres with omega = 0, 0.5, 0.95 and 1, with and without
  gas;
- black, Lambertian, Fresnel (n = 1.5 and 3+0.2i) and a custom non-specular,
  polarizing, non-reciprocal surface;
- 2 to 32 streams per hemisphere;
- RT4 `nstokes = 1` against the I-only embedding;
- the user-angle formal solution at off-node angles, for black and Lambertian
  surfaces.

**Not blind.** The test feeds VDISORT deliberately wrong inputs and asserts
that each misses RT4 by more than 100 times the tolerance:

- a stream, Stokes or full transpose of the phase matrix;
- an exchange of the same- and opposite-hemisphere quadrants;
- a stream or Stokes transpose of the surface.

The deviations are 6.7e-3 to 0.16 of max I, against tolerances below 1e-6.  Two
exchanges remain undetectable:

- RT4 requires mirror symmetry between the hemispheres, so exchanging the
  (down <- up) and (up <- down) quadrants changes nothing;
- the full transpose of a reciprocal phase matrix is the matrix itself,
  which is why the non-reciprocal case exists.

**What it does not validate.**

- U and V.  RT4 computes only [I, Q], so the sine (alpha = 1) system and the
  [U, V] reflection are never excited.  The RT3 comparison below covers U and
  V.
- Fourier modes m > 0, the solar beam, and delta-M, IMS and TMS.  RT4 is
  thermal-only and azimuthally symmetric.  The RT3 comparison below covers
  the beam and m > 0.
- Direction-dependent extinction or emission.  VDISORT cannot represent
  them.
- Upward user-angle radiances over a Fresnel surface.  This comparison is
  printed but not asserted, because it fails.  `fresnel_fourier_modes`
  reflects only into an outgoing mu equal to a quadrature node.  At mu = 0.35
  and 1 (8 streams per hemisphere, n = 3+0.2i), VDISORT's upward radiance at
  the surface is therefore the emission alone.  RT4 also includes the
  specular reflection R I_down, which is 22% and 13% of max I there.  The
  remaining difference is the barycentric interpolation of the
  angle-dependent emission from the nodes: 7e-4 and 1.4e-3 of B_s.  The
  downward user-angle radiances agree.

**MKL.** The top-level CMakeLists.txt enables Fortran only after LAPACK has
been found.  Enabled earlier, it makes CMake's FindBLAS link MKL's GNU layers
(`mkl_gf_lp64`, `mkl_gnu_thread`), whose multi-threaded `zgbsv` (MKL 2026.1)
segfaults for band widths KL = KU >= 65 independently of the matrix values.
VDISORT's boundary system has KL = KU = 6 NQuad - 1, so it crashed from
NQuad = 12 with two layers.

## Polarized comparison with RT3

`cpp.fast.vdisort-rt3-test` (`vdisort-rt3-comparison.cpp`) is built only with
`ENABLE_RT3=ON`.  It compares VDISORT with Evans' RT3, a polarized
doubling-adding solver with a solar beam and every Fourier azimuth mode
(`src/core/rt3`, `doc/arts/dev.rt3.rst`).  RT3 makes its own Fourier modes of
the phase matrix from the Legendre series of the six scattering-plane
elements, with its own rotations, FFT and beam pseudo-source.  Both solvers
get the same discrete problem:

- RT3's double-Gauss streams, which are VDISORT's (asserted);
- the same Legendre series.  The test builds VDISORT's ordinary Fourier
  coefficients C^m and S^m (no epsilon_m) from it by vector geometry
  (`lab-frame.h`, shared with the RT4 comparison), at RT3's azimuth samples.
  It combines the diffuse ones with `vdisort::combine_phase_matrices` and the
  beam ones with `vdisort::combine_beam_phase_matrices`, which V0 checks
  against the test's own `combine_beam()` (see the beam operator below);
- the same extinction, single-scattering albedo, solar beam
  (`beam_stokes = F / mu0`), Planck source linear in optical depth, sky and
  Lambertian surface.

VDISORT's radiance at azimuth `phi0 + psi` is RT3's at `psi`.  Every level
(top, each interface, bottom), every stream, both directions and the
azimuths 0, 30, 75, 90, 135, 180 and 250 deg are compared for I, Q, U and V.
The up- and downward fluxes of I (down including the direct beam) and of Q
are compared too.

RT3's error is first order in the thickness delta of its initial doubling
layer.  With a beam it is 0.41 to 0.68 delta / mu0, because RT3's
initial-layer beam source makes a relative error of delta / (2 mu0).
Without a beam it is 0.2 to 0.7 delta.  The tolerance is 10 times that.  The
difference, relative to max I (radiances) and max F (fluxes):

| Case | max_delta_tau = 1e-7: I / Q / U / V / F | RT3 Richardson-extrapolated: I / Q / U / V / F |
|---|---|---|
| R1 Rayleigh, tau 0.5, omega 0.95, mu0 0.6, A 0.1, nmu 8 | 5.5e-8 / 2.4e-8 / 3.0e-8 / 0 / 1.5e-8 | 3.4e-10 / 1.3e-10 / 1.3e-10 / 0 / 2.1e-12 |
| R2 Evans' mietest, tau 1, omega 0.99, mu0 0.2, A 0.1, nmu 8 | 2.0e-7 / 2.0e-8 / 1.3e-8 / 1.1e-10 / 8.9e-8 | 1.7e-10 / 7.4e-11 / 1.9e-11 / 2.8e-13 / 2.8e-11 |
| R2, nmu 12 | 1.9e-7 / 2.1e-8 / 1.3e-8 / 1.2e-10 / 8.9e-8 | 5.5e-10 / 1.7e-10 / 4.1e-11 / 5.1e-13 / 2.4e-11 |
| R3 Rayleigh / Mie / gas, solar + thermal, 3 um, A 0.25, nmu 8 | 8.0e-8 / 3.0e-8 / 4.3e-8 / 8.9e-11 / 2.2e-8 | 1.5e-10 / 5.4e-11 / 3.9e-11 / 2.0e-13 / 7.0e-12 |
| R3, thermal source only | 1.7e-8 / 1.2e-9 / 0 / 0 / 1.4e-8 | 2.3e-11 / 4.1e-12 / 0 / 0 / 1.2e-11 |
| R4 R2 with VDISORT phi0 = 1.1 | as R2 | as R2 |
| R5 R2 and R3 with nstokes 1, 2, 3 | at most 2.0e-7 (I) | at most 1.7e-10 (I) |
| R6 R2 and R3 with nmu 2, 4, 8, 16 | at most 2.0e-7 (I) | at most 1.1e-9 (I, R2 nmu 16) |
| R9 R3 with a forward peak, RT3 delta-M | 8.0e-8 / 3.0e-8 / 4.3e-8 / 8.1e-11 / 2.1e-8 | 1.3e-10 / 5.3e-11 / 3.9e-11 / 1.4e-13 / 1.2e-11 |

RT3's max |Q|, |U| and |V| / max |I| are up to 0.26, 0.50 and 6.3e-4, so
every component is exercised.  For R7 (R3, max_delta_tau 1e-5, 5e-6, 1e-6,
1e-7) the difference is 7.4e-6, 3.7e-6, 6.4e-7 and 8.0e-8 of max I.  It
halves exactly (2.000) and falls by 11.7 and 8.0 per decade.  The Richardson
residual is RT3's second-order error, 0.04 to 0.23 delta^2 / mu_min with
mu_min the smallest stream cosine: residual / delta^2 stays constant from
max_delta_tau = 1e-4 down to 1e-5.

The Fourier builder is checked first:

- its m = 0 Rayleigh matrix equals the closed form to 4e-16 (all 16
  elements, mu = 1 included);
- the single-scattering beam column synthesised from `combine_beam()` with
  VDISORT's field convention equals an independent dipole construction to
  7e-16, for phi0 = 0 and 1.1;
- for Evans' Mie series, which is regular at Theta = 0 and 180 deg only to
  1e-8, RT3's 32 azimuth samples and 1024 midpoints give coefficients that
  differ by 2.6e-10 relative.  That is why the builder uses RT3's samples.

**What it validates.**

- The solar beam: its pseudo-source in every Fourier mode, its attenuation,
  `beam_stokes` as the irradiance normal to the beam, `mu0`, `phi0` and its
  Lambertian reflection.
- Fourier modes m > 0 (up to m = 11) of the diffuse operator, and
  `combine_phase_matrices` (Lin et al. Eq. 81).
- U and V: the sine (alpha = 1) system, U at non-principal-plane azimuths
  (up to 0.5 of max I), and V generated through F34.
- Rayleigh and Evans' Mie matrix (F34 != 0); multilayer atmospheres with gas
  and a gas-only layer; solar and thermal sources together, with sky and
  surface emission.
- nstokes 1, 2 and 3 against the leading block of the phase matrix, and 2 to
  16 streams per hemisphere.
- VDISORT given RT3's delta-M scaled problem (scaled tau, omega, Legendre
  series, beam attenuation).

**Not blind.** Each deliberately wrong input or mapping misses RT3 by
39000 to 400000 times the tolerance:

| Mistake (R2 / R4, nmu 8) | I / Q / U / V |
|---|---|
| S^m sign flipped (diffuse and beam) | 2.0e-7 / 2.0e-8 / 0.118 / 4.1e-4 |
| Z at -phi' in the Fourier transform | 2.0e-7 / 2.0e-8 / 0.118 / 4.1e-4 |
| VDISORT at phi0 - psi (azimuth sense reversed) | 2.0e-7 / 2.0e-8 / 0.118 / 4.1e-4 |
| VDISORT at psi - phi0 (phi0 sign) | 0.870 / 0.120 / 0.123 / 4.0e-4 |
| diffuse C^m, S^m (m > 0) times epsilon = 2 | 0.135 / 0.081 / 0.115 / 9.0e-4 |
| beam C^m, S^m (m > 0) times epsilon = 2 | 0.630 / 0.070 / 0.059 / 2.1e-4 |
| cosine and sine systems swapped | 1.19 / 0.168 / 0.115 / 3.0e-4 |

The first three are the same error, U, V -> -U, -V (the similarity transform
diag(1, 1, -1, -1) of the combined systems), so U and V carry their
detection.  The fluxes catch none of them except the swap.

**Beam operator.**  The comparison sets VDISORT up through
`vdisort::combine_beam_phase_matrices`, and V0 asserts that it equals the
test's own derivation exactly.  The beam is a delta in azimuth at phi0,
`(1 / 2 pi) sum_m eps_m cos m(phi0 - phi)`, i.e. cosine terms only for every
Stokes component, so the combined beam matrices are:

- for the cosine system (I^c, Q^c, U^s, V^s), rows I, Q of C^m(mu, -mu0) and
  rows U, V of S^m;
- for the sine system (I^s, Q^s, U^c, V^c), rows I, Q of S^m and rows U, V of
  C^m;

in all four columns.  This differs from the diffuse combination of Eq. 81,
which would put C^m_{IQ,I} and -S^m_{UV,I} into the sine system for an
unpolarized beam; that error is antisymmetric in azimuth, leaves m = 0 and the
fluxes unchanged, and reaches 0.18 to 0.51 of max I in R1 to R3.

**What it does not validate.**

- IMS and TMS, and delta-M-plus.  The delta-M case gives VDISORT RT3's
  scaled problem; VDISORT's own corrected delta-M is not compared.
- Fresnel, Cox-Munk and the other BRDFs with a beam.  RT3 allows only a
  Lambertian surface with a beam.  Thermal Fresnel at m = 0 is covered by
  RT4.
- Off-node user angles (`ungridded_u_user`): RT3's extra angles need its
  gauss quadrature, which VDISORT does not have.
- The absolute sign of V.  V agrees, but it is generated through F34, so
  this tests that the two solvers use F34 consistently, not V's sign.
  Neither RT3 nor this test pins V against an external derivation.
- A polarized beam (`beam_stokes` Q, U, V), because RT3's beam is
  unpolarized.  Only column I of the beam matrices acts.
- Phase matrices that are not of the six-element form of randomly oriented
  particles with a plane of symmetry.

## The three solvers on shared problems, and CI

`tests/core/disort/vdisort-polradtran.rt3.rt4.py` runs VDISORT, RT4 and RT3
from `pyarts3` on ARTS atmospheres through the same path builders as
`cpp.fast.vdisort-arts-comparison` (`vdisort.main_data_from_path`,
`rt4.problem_from_path`, `rt3.problem_from_path`).  Run without
`ARTS_HEADLESS`, it draws the solutions' plots (`pyarts3.plots.cppvdisort`,
`RT4Result` and `RT3Result`) on shared axes.  All solvers get the same double-Gauss streams.  RT4 runs
the thermal problems with nstokes <= 2.  The cases are:

| Case | Solvers |
|---|---|
| thermal, Rayleigh + Mie drops, Fresnel 3+0.2i, nstokes 2 | VDISORT, RT4, RT3 |
| thermal, Henyey-Greenstein, Fresnel 3+0.2i, nstokes 2 | VDISORT, RT4, RT3 |
| thermal, Rayleigh + Mie drops, Lambertian 0.3, nstokes 2 | VDISORT, RT4, RT3 |
| thermal, isotropic + Henyey-Greenstein (unpolarized), Lambertian 0.3, nstokes 1 | VDISORT, RT4, RT3 |
| solar mu0 = 0.6 + thermal, Rayleigh + Henyey-Greenstein + Mie drops, Lambertian 0.3, 8 modes, nstokes 4, 6 azimuths | VDISORT, RT3 |

Every pair must agree to 10 max_delta_tau / mu0 of max I (1e-6, and 1.7e-6
with the beam).  The measured values are 2e-7 to 4e-7 for VDISORT against
RT3 or RT4, and 1e-11 to 1.4e-8 for RT4 against RT3.  The nstokes 1 case is
unpolarized, so VDISORT's I is the reference for the scalar I of RT3 and
RT4.  The Henyey-Greenstein case over the polarizing Fresnel surface needs a
scattering matrix that is regular at forward and backward scattering (see
`doc/arts/dev.rt3.rst`, azimuth sampling), as ARTS's is.  The C++
comparisons above carry the convergence, Richardson and non-blindness
evidence.

RT3 and RT4 are built by default when a Fortran compiler is found.  CI
(`.github/workflows/build-test.yml`) still sets `-DENABLE_RT4=ON
-DENABLE_RT3=ON` explicitly, through the `polradtran` matrix key, for every
job with a Fortran compiler (all Linux jobs, including LGPL, and macOS).  A
job that loses its compiler then fails to configure instead of silently
skipping these tests.
`check` then runs `cpp.fast.rt3-test`, `cpp.fast.rt4-test`,
`cpp.fast.vdisort-rt3-test`, `cpp.fast.vdisort-rt4-test`, the closed-form
Python tests of both bindings, and this three-solver test.  Windows has no
Fortran compiler in its environment and runs none of them.

## The three solvers on ARTS data

`cpp.fast.vdisort-arts-test` (always built) and
`cpp.fast.vdisort-arts-comparison` (`ENABLE_RT3` and `ENABLE_RT4`) test the
ARTS-native input builders `vdisort_arts.h`, `rt3_arts.h` and `rt4_arts.h`
(see `doc/arts/dev.disort.rst`, `dev.rt3.rst`, `dev.rt4.rst`).

| test | reference | measured |
|---|---|---|
| V1: Rayleigh `GasScatterer` C^m, S^m, m = 0..3, diffuse and beam, mu = +-1 included | dipole Jones-matrix closed forms | 1.6e-15 |
| V2: ARTS's laboratory-frame Z of a polarizing Mie particle, 300 directions | vector geometry, mu = cos(za), phi = -aa, same F | 6.7e-15 (other azimuth sense: 2; F34 negated: 0.35) |
| V2: exact and snapped near-forward pairs | Z = F on the diagonal (forward-scattered Q) | 1.2e-8 |
| V3: `vdisort::scattering_optics` of that particle | Fourier modes of ARTS's laboratory-frame Z | 1.9e-15 |
| thermal, 89 GHz, Rayleigh + Mie cloud + gas, Lambertian / Fresnel | RT3 vs VDISORT | 5.4e-7 / 5.2e-7 of max I (tolerance 1e-6) |
| same | RT4's layer phase matrices vs VDISORT's (different input routes) | 1.2e-15 relative (tolerance 1e-12) |
| same | RT4 vs RT3 (identical doubling) | 5.3e-10 / 5.0e-10 (tolerance 1e-8) |
| same | RT4 vs VDISORT | 5.4e-7 / 5.2e-7 (tolerance 1e-6) |
| solar mu0 = 0.6 + thermal, 8 modes, nstokes 4 | RT3 vs VDISORT | 5.8e-7 (tolerance 1.7e-6) |

The RT3 - VDISORT difference falls by 11.7 per decade of `max_delta_tau`.
