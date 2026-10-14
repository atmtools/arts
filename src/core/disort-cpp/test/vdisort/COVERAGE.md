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
| 6 | a–h | ported | All active direct, diffuse-down, upward, and `DFDT` flux references are asserted. The port covers the transparent limit, absorption, Lambertian and Hapke reflection, top/bottom emission, and linear atmospheric thermal sources; Q, U, and V remain zero. The thermally emitting Hapke cases 6e–h are checked against cDISORT's references (see the scalar port's coverage); VDISORT derives the emissivity from its streams, cDISORT from a fine quadrature, so 6f–h use a 1% tolerance (6h's DFDT at the surface differs by 0.4%), while all other cases use 0.1%. No active radiance references exist for this problem. |
| 7 | a–e | ported | All active direct, diffuse-down, upward, and `DFDT` fluxes are asserted for scattering plus linear atmospheric emission, beam and isotropic illumination, and black, perfectly reflecting Lambertian, and Hapke surfaces. Tiny case-7b fluxes use combined relative and absolute tolerances; the Hapke case 7e is checked against cDISORT's references at 0.2%, as the other cases, and Q, U, and V remain zero. No active radiance references exist for this problem. |
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

`cpp.fast.vdisort-rt4-test` (`vdisort-rt4-comparison.cpp`) compares VDISORT with Evans' RT4, a polarized
doubling-adding solver (`src/core/polradtran/rt4`, `doc/arts/dev.rt4.rst`).  This is the
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
- the user-angle formal solution at off-node angles, for black, Lambertian
  and Fresnel surfaces.  A Fresnel surface is the specular part
  `BDRF::specular` = R(mu) of its modes, without a reflection kernel.  VDISORT
  applies it between equal streams and, at an upward user angle mu, to the
  downward user-angle radiance at -mu at the surface.  An upward user angle
  therefore needs its downward partner, and the test asserts that VDISORT
  refuses one without it.  The upward boundary is the surface's Planck
  radiance, constant in mu, so its interpolation to a user angle is exact,
  and VDISORT transmits it with 1 - R(mu) at the user angle itself: the
  upward user-angle radiances must agree with RT4 to its doubling error at
  every level.

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
  the beam and m > 0, and the polarized delta-M test covers IMS and TMS.
- Direction-dependent extinction or emission.  VDISORT cannot represent
  them.

**MKL.** The top-level CMakeLists.txt enables Fortran (for the optional
T-matrix code) only after LAPACK has been found.  Enabled earlier, it makes CMake's FindBLAS link MKL's GNU layers
(`mkl_gf_lp64`, `mkl_gnu_thread`), whose multi-threaded `zgbsv` (MKL 2026.1)
segfaults for band widths KL = KU >= 65 independently of the matrix values.
VDISORT's boundary system has KL = KU = 6 NQuad - 1, so it crashed from
NQuad = 12 with two layers.

## Polarized comparison with RT3

`cpp.fast.vdisort-rt3-test` (`vdisort-rt3-comparison.cpp`) compares VDISORT with Evans' RT3, a polarized
doubling-adding solver with a solar beam and every Fourier azimuth mode
(`src/core/polradtran/rt3`, `doc/arts/dev.rt3.rst`).  RT3 makes its own Fourier modes of
the phase matrix from the Legendre series of the six scattering-plane
elements, with its own rotations, FFT and beam pseudo-source.  Both solvers
get the same discrete problem:

- RT3's double-Gauss streams, which are VDISORT's (asserted);
- the same Legendre series.  The test builds VDISORT's ordinary Fourier
  coefficients C^m and S^m (no epsilon_m) from it with ARTS's polarization
  geometry (`lab-frame.h`, shared with the RT4 comparison; the geometry is
  tested against electric fields, see below), at RT3's azimuth samples.
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
  scaled problem; VDISORT's own corrected delta-M is not compared with RT3.
  The polarized delta-M test below checks it against VDISORT without
  delta-M.  Polarized delta-M-plus is not covered.
- Fresnel, Cox-Munk and the other BRDFs with a beam.  RT3 allows only a
  Lambertian surface with a beam.  Thermal Fresnel at m = 0 is covered by
  RT4, and the reflection of a beam by the polarized Cox-Munk surface by the
  field-level test below.
- The absolute sign of V.  V agrees, but it is generated through F34, so
  this tests that the two solvers use F34 consistently, not V's sign.  The
  field-level tests below pin V = 2 Im(E_v E_h*) in the meridional basis for
  a polarized beam and for Cox-Munk reflection; which circular polarization
  that makes positive is not tested.
- A polarized beam (Q, U, V of the beam), because RT3's beam is
  unpolarized.  Only column I of the beam matrices acts.  The field-level
  test below covers it.
- Phase matrices that are not of the six-element form of randomly oriented
  particles with a plane of symmetry.

### User angles between the streams (R10)

RT3 takes zero-weight extra angles with its double-Gauss rule, so VDISORT's
formal solution at user angles can be compared on identical streams.  R10
gives RT3 seven extra angles between the streams (mu = 0.05, 0.13, 0.31, 0.5,
0.77, 0.93, 0.99) and evaluates VDISORT at the same angles, upward and
downward, at every level and at the azimuths of the main comparison.
Relative to max I:

| Case | max_delta_tau = 1e-7: I / Q / U / V | RT3 Richardson-extrapolated: I / Q / U / V |
|---|---|---|
| R2 (Mie layer) | 2.1e-7 / 2.4e-8 / 1.4e-8 / 8.3e-11 (0.69 delta / mu0) | 7.6e-11 / 2.8e-11 / 8.2e-12 / 2.0e-13 |
| R3 (Rayleigh / Mie / gas, solar + thermal) | 8.2e-8 / 2.5e-8 / 4.3e-8 / 9.0e-11 (0.49 delta / mu0) | 3.8e-11 / 1.6e-11 / 1.1e-11 / 5.8e-14 |

Between the streams the formal solution agrees with RT3 to RT3's doubling
error, as on the streams, with the same tolerances.

### Evans' benchmark settings (E)

Evans' original programs from the tar (rt3.f, rt4.f, scatcnv.f) reproduce the
tables of all four of his scripts in `3rdparty/polradtran`: the RT3 ones to one
unit in the last printed digit, the RT4 ones (brightness temperatures with two
decimals) exactly.  ARTS keeps none of the Fortran.  `cpp.fast.polradtran-rt4-arts` gives ARTS's RT4 the
optics rt4.f reads for runtestc (cl340d14.dda) and reproduces its table to
0.5 of the last printed digit (0.005 K); runtestr, whose scattering file
scatcnv made, is solved from its Mie series by
`tests/core/disort/evans-benchmarks.rt3.rt4.py`.

Evans' two RT3 scripts, `runmietest` (the Mie case of Evans and Stephens
1991: tau 1, omega 0.99, mu0 0.2, A 0.1, nmu 8, aziorder 8) and `runtesta`
(Rayleigh over Mie with gas, solar and thermal at 3 um, A 0.25, nmu 4,
aziorder 4), are run from the scripts themselves
(`src/tests/polradtran/evans-scripts.h` reads their here-documents).  Evans'
original program reproduced both tables to one unit in the last printed
digit; the values that are zero by symmetry are REAL*4 round-off below 1e-8
of max I in the tables.  Part E of `cpp.fast.vdisort-rt3-test` gives RT3
the same settings and VDISORT the same physical problem.  VDISORT has
double-Gauss streams only, so it is evaluated at Evans' Gauss angles with
its formal solution (`ungridded_u_user`).

Measured, relative to max I (radiances) and max F (I fluxes):

| | runmietest I / Q / U / V / F | runtesta I / Q / U / V / F |
|---|---|---|
| RT3 at Evans' settings vs his table | 6.3e-7 / 5.5e-7 / 6.5e-8 / 5.0e-10 / 7.5e-7 | 5.5e-7 / 5.0e-7 / 3.8e-7 / 5.5e-10 / 6.1e-7 |
| VDISORT 16 vs 32 streams per hemisphere, at Evans' angles | 1.4e-7 / 5.7e-8 / 6.9e-9 / 1.8e-10 / 2.2e-8 | 5.6e-8 / 2.8e-8 / 8.6e-9 / 1.3e-10 / 8.2e-9 |
| VDISORT at Evans' angles vs his table | 8.9e-3 / 2.4e-3 / 5.0e-4 / 7.0e-6 / 1.5e-3 | 1.3e-2 / 5.8e-3 / 3.4e-3 / 4.2e-6 / 3.9e-3 |
| RT3 double-Gauss nmu 16 vs VDISORT on the same streams | 1.9e-7 / 1.7e-8 / 5.5e-9 / 3.2e-11 / 8.9e-8 | 8.7e-8 / 4.1e-8 / 4.1e-8 / 5.0e-11 / 3.0e-8 |

So VDISORT is not Evans' tables: they differ by up to 0.9 and 1.3 % of max
I, at grazing upwelling angles at the top (mu = 0.095 and 0.183).  That is
the error of the Gauss quadrature of the tables, which handles the
discontinuity of the radiance at the horizon poorly; double-Gauss, which
VDISORT and DISORT use, does not have it.  RT3 itself shows this: with
Gauss quadrature its distance from VDISORT (at RT3's nodes) falls steadily
with nmu, about like 1 / nmu,

| RT3 Gauss nmu | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|
| runmietest, max I difference | 2.4e-2 | 1.3e-2 | 8.9e-3 | 5.1e-3 | 3.5e-3 |
| runtesta, max I difference | 1.3e-2 | 7.8e-3 | 5.4e-3 | 3.3e-3 | 2.3e-3 |

and at Evans' nmu it is the table's distance, since RT3 there is the table.
RT3 with double-Gauss quadrature agrees with VDISORT to its doubling error.
The test asserts the first, second and fourth rows of the first table, and
that the Gauss sequence decreases.

Evans' RT4 script `runtestr` (a rain layer of spherical drops at 85 GHz
over water, Fresnel n = 3.17 - 1.75i, 8 Gauss streams, V and H brightness
temperatures) is a randomly oriented problem, so RT3 (aziorder 0) and
VDISORT solve it from the Mie Legendre series that scatcnv converts for
RT4:

VDISORT is evaluated at Evans' angles by its formal solution; the upward
radiance reflects the downward one at the same angle (`BDRF::specular`).

| runtestr, every level, angle and direction | |
|---|---|
| RT3 at Evans' settings vs his RT4 table | 0.005 K (asserted, 0.01 K) |
| VDISORT, 16 vs 32 streams, at Evans' angles | 0.0001 K (asserted, 0.01 K) |
| VDISORT (32 streams) at Evans' angles vs his table | 0.73 K |
| RT3 Gauss nmu 4, 8, 16 vs VDISORT at RT3's nodes | 1.48, 0.73, 0.36 K (asserted to decrease) |
| RT3 double-Gauss nmu 8 and 16 vs VDISORT on the same streams | 3.7e-8, 5.3e-8 of max I (asserted, 2e-6) |

The difference from the table is again its Gauss quadrature: RT3 with Gauss
quadrature halves its distance from VDISORT when nmu doubles, and at Evans'
nmu 8 that distance (0.73 K) is the table's.  It is largest in the
upwelling H radiance at grazing angles, where the Fresnel reflection of the
downwelling radiance dominates.

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

RT3 and RT4 are C++ ports of Evans' Fortran, which ARTS no longer keeps, so
they are always built and every CI job runs `cpp.fast.rt3-test`,
`cpp.fast.rt4-test`, `cpp.fast.vdisort-rt3-test`,
`cpp.fast.vdisort-rt4-test`, the closed-form Python tests of both bindings,
and this three-solver test.

## The three solvers on ARTS data

`cpp.fast.vdisort-arts-test` (always built) and
`cpp.fast.vdisort-arts-comparison` test the
ARTS-native input builders `vdisort_arts.h`, `rt3_arts.h` and `rt4_arts.h`
(see `doc/arts/dev.disort.rst`, `dev.rt3.rst`, `dev.rt4.rst`).

| test | reference |
|---|---|
| V1: Rayleigh `GasScatterer` C^m, S^m, m = 0..3, diffuse and beam, mu = +-1 included | dipole Jones-matrix closed forms |
| V2: ARTS's laboratory-frame Z of a polarizing Mie particle, generic directions | the polarization geometry with mu = cos(za), phi = -aa and the same F; the other azimuth sense and a negated F34 must miss.  Both sides use the shared geometry (tested against fields below), so this pins the direction convention and the sign of F34 |
| V2: exact and near-forward pairs | Z = F on the diagonal (forward-scattered Q) |
| V3: VDISORT's coefficients of that particle from ARTS's Fourier modes | Fourier modes of the geometry's Z (`lab-frame.h`) |
| V4: the path builder | its inputs |
| V5: gridded ARO data (Rayleigh's laboratory frame on a grid) | `GasScatterer`'s exact modes, to the interpolation error of the grid |
| V6: Henyey-Greenstein, m = 0..3, diffuse and beam | the addition theorem |
| thermal, 89 GHz, Rayleigh + Mie cloud + gas, Lambertian / Fresnel | RT3 vs VDISORT: their doubling error |
| same | RT4's layer phase matrices vs VDISORT's (both the species' Fourier modes): rounding |
| same | RT4 vs RT3 (identical doubling; Fourier modes against Legendre series): rounding |
| same | RT4 vs VDISORT: RT4's doubling error |
| solar beam + thermal, 8 modes, nstokes 4 | RT3 vs VDISORT: RT3's doubling error with a beam |

The tests print the measured values and enforce the tolerances.  The RT3 - VDISORT difference must fall with `max_delta_tau`, as RT3's first-order doubling error does.

## Polarized delta-M with IMS and TMS

`cpp.fast.vdisort-delta-m-test` (`vdisort-delta-m-test.cpp`) checks VDISORT's
classical delta-M with the cached polarized TMS and IMS corrections in all
four Stokes components.  There is no external polarized reference, so the
reference is VDISORT itself without delta-M, on enough streams to carry the
whole phase-matrix series.  It shares the transport solver, which the RT3 and
RT4 comparisons validate, but not the scaling, the removed peak or the
corrections.

**Problem.**  One layer (omega 0.9) under an unpolarized beam (mu0 0.5) with
a forward-peaked, polarizing phase matrix: Henyey-Greenstein F11 (g = 0.75)
with F12 and F34 of a few tenths of F11 away from the forward and backward
directions.  The matrix is built in angle space so that it is regular at
forward and backward scattering (F12 = F34 = 0, F22 = F33 forward, F22 = -F33
backward).  Otherwise its laboratory-frame Q and U elements jump at the
forward peak and no discrete-ordinate solution converges there.  The delta-M
split is DISORT's classical one applied to the whole matrix, exact for every
Legendre degree.  User directions lie near and away from the forward peak
(upward at the top, downward at the bottom, five azimuths).

**Reference convergence.**  64 against 48 streams: 1.0e-5 in I and below
1e-7 in Q, U and V, relative to max I.  The reference has max |Q|, |U|, |V| /
max I of 1.9e-2, 1.8e-2 and 4.2e-5, so every component is exercised.

**Results** (optical depth 1), maximum deviation from the reference relative
to max I:

| NQuad | f | delta-M: I / Q / U / V | with TMS and IMS: I / Q / U / V |
|---|---|---|---|
| 8 | 0.10 | 4.2e-1 / 3.2e-3 / 1.2e-3 / 6.8e-6 | 1.6e-3 / 1.6e-3 / 3.9e-4 / 6.8e-6 |
| 12 | 3.2e-2 | 2.2e-1 / 1.2e-3 / 5.5e-4 / 2.5e-6 | 3.7e-4 / 3.8e-4 / 4.8e-5 / 2.5e-6 |
| 16 | 1.0e-2 | 1.1e-1 / 4.4e-4 / 2.4e-4 / 5.4e-7 | 7.7e-5 / 9.1e-5 / 8.1e-6 / 5.4e-7 |
| 24 | 1.0e-3 | 2.0e-2 / 4.5e-5 / 1.8e-5 / 1.3e-8 | 5.2e-6 / 6.2e-6 / 8.6e-7 / 1.3e-8 |

The test requires the corrections to reduce the error of I by more than 100
times and those of Q and U by more than 1.5 times, not to increase that of V,
and the corrected error to fall with NQuad.  V's error comes from multiple
scattering through the small F34, which the single-scattering corrections do
not address.

In a thin layer (optical depth 0.02, reference on 48 streams), where single
scattering dominates, TMS alone must remove most of delta-M's error: by more than 100 times in I and 10
times in Q and U.  Measured at NQuad 8 and 16: I 5.5e-1 and 1.5e-1 become
2.1e-3 and 4.0e-4, Q 4.6e-3 and 7.0e-4 become 1.4e-4 and 2.5e-5, and U
2.0e-3 and 4.1e-4 become 1.6e-4 and 2.8e-5.

**TMS.**  TMS subtracts the beam's single scattering as the solver carries
it, the NFourier-mode synthesis of the transport beam phase matrices, and
adds the exact one.  For F11 that is the whole transport matrix.  The
polarized elements also contain the rotations between the scattering and
meridional planes, which are not band-limited in azimuth; subtracting the
full transport matrix made Q 2 to 4 times worse than delta-M without
corrections.

**Cost.**  The corrections are cheap: VDISORT's cache arithmetic takes about
15 ms for the whole test and applying them 0.07 s.  The test takes about
20 s, most of it the reference solves (64 and 48 streams at optical depth 1,
48 for the thin layer).  With one OpenMP thread it takes 10 s; the difference
is the threaded BLAS waking its threads for every small call inside the
eigensolutions.

**Not covered.**  Delta-M-plus, multilayer atmospheres, surfaces and thermal
sources with the corrections.

## Polarized beam and surface against electric fields

`tests/core/disort/vdisort-polarized-beam-surface.py` covers through
`pyarts3` what the RT3 and RT4 comparisons cannot: a polarized beam, and a
polarized surface reflecting a beam.  Its references are built from electric
fields with ARTS's polarization primitives (the meridional basis and the
Mueller matrix of a Jones matrix, tested against fields below), not from
VDISORT's phase matrices or rtepack's reflection matrices.  Relative to max I:

| Part | Reference | Deviation | Not blind |
|---|---|---|---|
| A: polarized beam, single scattering in a thin Rayleigh layer | dipole field projected on the meridional bases | 4.2e-11 to 6.8e-11 in I, Q, U, V | an unpolarized beam misses by 0.48, the mirrored azimuth by 1.38 |
| B: VDISORT's raw polarized Cox-Munk reflection, n = 1.34, 1.5+0.3i, 7+3i | Fresnel fields of the specular facet in its plane of incidence, divided by M00 | 1.2e-15 to 1.4e-15 | |
| C: beam reflected by Cox-Munk under a transparent atmosphere, unpolarized and polarized beam | mu0 F times the facet reflection, Fourier-truncated as the solver's surface modes are | at most 5.2e-11 | the mirrored relative azimuth misses by 0.10 and 0.11 |

In C the reflected max |Q|, |U| and |V| / max I reach 0.80, 0.17 and 0.44, so
all components are exercised.  The scalar part of Cox-Munk (slope
distribution and shadowing) is DISORT's, validated by the Fortran problems;
B tests the polarization structure.

## Shared polarization geometry

ARTS has one polarization geometry (`doc/arts/concept.polarization.rst`): the
meridional basis of a ray, the Stokes rotation between two bases, the Mueller
matrix of a Jones matrix, the rotation of a scattering-plane matrix into the
laboratory frame, and Fresnel reflection by a facet.  The VDISORT tests'
laboratory-frame phase matrices (`lab-frame.h`, used by the RT3, RT4, ARTS
and delta-M tests), ARTS's conversion of randomly oriented scattering data
to the laboratory frame, VDISORT's Cox-Munk surface and ARTS's non-specular
surface reflection all use it.

`cpp.fast.test_rtepack_polarization` tests it against electric fields only.
For random pure states it compares the Mueller matrix applied to the Stokes
vector of a field with the Stokes vector of the transformed field, both from
the field components:

- the basis is orthonormal, right-handed with the ray, has v on the upward
  side, and at the poles is the limit along the given meridian;
- the rotation between two bases of a ray, and refusal of a basis of the
  other handedness;
- the Mueller matrix of random Jones matrices;
- the laboratory-frame phase matrix of scattering amplitudes diag(S2, S1)
  in the scattering plane, real (Rayleigh) and complex (giving F34), at
  generic directions and exactly forward and backward, at the poles too;
- Fresnel reflection by a facet for a complex refractive index, and a
  perfect mirror.

All agree to 1e-13.  The mirrored phase matrix, and the mirror with V
negated once more, miss by more than 1e-2.

Before ARTS's scattering data moved to this geometry, the two were compared
on a dense grid of directions including the poles, the principal plane and
nearly forward and backward pairs.  Away from the poles they agree to
rounding.  Near the poles the previous angle formulas had errors up to 7e-9,
against a 50-digit reference that the vector geometry matches to 1e-15.  For
pairs within 1e-6 rad of forward or backward scattering near the poles, the
previous formulas set the phase matrix to the scattering matrix and dropped
the rotation between the two meridional bases, up to 1e-2 in Q-U mixing in
the case checked.

## Reusing one solver

`cpp.fast.vdisort-reuse-test` (`vdisort-reuse-test.cpp`) reuses one solver
for several problems, as a thread does for many frequencies, and requires
every radiance and flux to equal those of a fresh solver for the same
problem bit for bit.  The update allocates nothing and resets only what is
accumulated or written conditionally.  The problems change every input
through the solver's accessors and cover each conditional path:

- a layer without scattering: no scattered beam source, but emission;
- an exactly conservative layer: the Jordan pair, and no emission;
- a polarized beam, and then no beam;
- a surface in the first Fourier mode only.

Removing the reset of the beam terms of a layer without scattering makes the
reused radiances miss by 0.14 of their scale.

## The half-size eigenproblem and the real boundary system

VDISORT diagonalizes a mirror-symmetric layer and Fourier mode through the
half-size product problem (except in the cosine system of order zero), and
solves the boundary system in real arithmetic (see
`doc/arts/concept.disort.rst`).  Against the previous solver, which used the
full eigenproblem and a complex boundary system everywhere, every comparison
in this document is unchanged to its printed digits or at the rounding level:
the radiances of the reference capture by at most 1.3e-14 of their maximum
and the fluxes by 4e-16, and the RT3 deviations in their third or fourth
digit at the 1e-10 level.  Rayleigh's V, which is zero, has round-off of
1.8e-14 instead of 2.9e-15 of max I.

`cpp.fast.vdisort-mirror-test` (`vdisort-mirror-test.cpp`) checks that an
asymmetric medium takes the full problem.  A problem driven only by boundary
radiance reads, in the half-size problem, nothing but the upward rows of the
operator.  Changing the downward rows of the phase matrix by a relative 1e-9
in every mode but the cosine system of order zero must therefore change the
radiances, which it does by 7.9e-13 of their maximum: taking the half-size
problem would have reproduced the symmetric result bit for bit, as a
symmetry check that accepts everything shows.  The agreement to that level
is also that of the half-size and full eigensolutions.

## Surface emission by Kirchhoff's law

Both solvers take the radiance beneath the surface as their upward boundary
and transmit it with 1 - R, R the surface reflection on the streams with its
specular part; for the Planck radiance this is Kirchhoff's emission (see
`doc/arts/concept.disort.rst`).  `cpp.fast.vdisort-kirchhoff-test`
(`kirchhoff-test.cpp`) puts an isothermal, scattering and emitting atmosphere
between a sky and a surface at its temperature.  The radiance must be B e_I
on every stream, depth and azimuth and at user angles in both hemispheres:

- DISORT over Lambertian, Hapke, shadowed Cox-Munk, RPV and Ross-Li surfaces
  and a Cox-Munk/Lambertian mixture, with Henyey-Greenstein scattering: 2.2e-15;
- VDISORT over Lambertian, Fresnel (specular), polarized shadowed Cox-Munk,
  Fresnel/Lambertian and depolarizing Hapke surfaces, with polarized Rayleigh
  scattering: I to 6.7e-16 and Q, U, V to 1.1e-16 of B.

Given the emission (1 - R) B as their boundary instead, the convention before,
the solvers miss B by 0.21 (DISORT, Lambertian) and 0.20 (VDISORT, Fresnel).

The Fortran references of the thermally emitting Hapke surface (problems 6e-h
and 7e) are cDISORT's: DISORT 4.0.99's correspond to the emissivity
1 - rho / pi, rho the hemispherical reflectance of the Hapke BRDF that reflects
the beam in 6d, where the two agree.  Both solvers reproduce cDISORT's to
4e-4, the difference of the stream and fine quadratures of the emissivity.
