RT4 reference solver
====================

RT4 is Evans' polarized doubling-adding solver from the polradtran package.
ARTS 3 keeps it as an external reference for validating other solvers, in
particular the polarized discrete-ordinate solver VDISORT.  It solves the
thermal radiative transfer equation in a plane-parallel, azimuthally
symmetric medium, for the Stokes components [I] or [I, Q].  Its equations
are in :doc:`concept.rt4`; its sibling RT3 is described in :doc:`dev.rt3`.

The C++ is in ``src/core/polradtran/rt4``, and what RT4 shares with RT3 is
in ``src/core/polradtran``.  The Python bindings are the submodules
``pyarts3.arts.rt4`` and ``pyarts3.arts.polradtran``.  There are no workspace
methods, variables or agendas.  The headers and the docstrings of the
bindings are the reference for names, fields, defaults and array layouts.
This page describes the design and the conventions.

Provenance
----------

* K. F. Evans' polradtran is distributed from
  https://nit.coloradolinux.com/polrad.html under the MIT licence.  RT4 is
  described briefly by :cite:t:`Evans1995`.  ``3rdparty/polradtran`` keeps
  Evans' licence, his ``README``, his benchmark scripts and the scattering
  file one of them reads, unchanged; its ``README.ARTS`` gives their origin.
* ARTS keeps no Fortran: RT4 is a C++ port of Evans' ``RADTRANO`` and of the
  routines it calls.
* ARTS 2.6 extended ``RADTRANO`` (J. Mendrok and others), and the port keeps
  these extensions: the optics are passed in memory, with an optics set per
  layer; the ground can be an external reflection operator and emission; the
  radiances are output at every level; and zero-weight extra angles follow
  the quadrature streams.

Design of the port
------------------

* **Evans' routines, one to one.**  Each Fortran subroutine is a C++
  function with its name in lower case, in a file named after Evans' file.
  The C++ follows the Fortran step by step.  The counts that the Fortran
  passes are the extents of the arrays, its ``STOP`` checks throw, and its
  static array limits are gone: every array is sized to the problem.
* **Shared with RT3.**  Some routines of RT4 are RT3's line by line.  Others
  are RT3's in the azimuthal mode 0 with the thermal source alone: the
  doubling, the adding, the internal radiance, the non-scattering layer, the
  sky radiance and the Lambertian and Fresnel grounds.  Each of these exists
  once, in the ``polradtran`` namespace, and RT4 calls it in mode 0.  The
  level loop and the choice of the initial doubling sublayer are shared in
  the same way.
* **ARTS numerics in place of Evans'.**

  * The quadratures are ARTS's rules.
  * The Planck function is ARTS's, in SI units, with :math:`e^x-1`
    evaluated without cancellation.
  * The Fresnel amplitudes are ARTS's.  They use the complex Snell's law,
    so absorbing grounds are exact.
  * Evans' matrix helpers are matpack, BLAS and LAPACK.  Each product is a
    GEMM or GEMV call whose scaling factors absorb the identity, subtraction
    or addition around it, and the inverse is LAPACK's.
  * Matrices keep the Fortran's column-major storage.  A matpack matrix
    therefore holds the transpose of Evans' matrix, and the products are
    written for that.

  These replacements change the results by rounding, amplified by the
  conditioning of the doubling (see `Validation of the port`_).
* **SI units.**  The solver takes the frequency and returns radiances per
  hertz.  There is no wavelength or per-micrometre unit anywhere.
* **The ground is data.**  Every ground of RT4 makes the same surface layer
  except for its reflection back up and its own radiance.  The solver core
  therefore takes those two arrays, and a separate step makes them for each
  kind of ground.
* **Work data.**  All scratch arrays of the solver, and of the routines it
  calls, live in one work-data object, which RT3's extends.  The solver sizes
  it and allocates only where an array grows, so keeping one over repeated
  calls with the same streams, as over frequency, allocates nothing.  The
  solver keeps no other state, so concurrent solves are safe.

Validation of the port
----------------------

* The port was compared with the Fortran, routine by routine and as a whole,
  on random inputs that exercise every branch.
* A step meant to leave the numbers alone was checked bit for bit.  A step
  that replaced numerics (quadratures, Planck function, LAPACK's inverse,
  fused BLAS products) had its effect measured.
* Comparisons of two implementations allow a few units of rounding, times
  what the computation amplifies rounding by.  That is :math:`2^n` for a
  layer doubled :math:`n` times and the condition number for an inverse.
  Any replaced numerics add a separate term.  Mathematically equivalent
  evaluation orders (FMA, vectorised libm functions, other BLAS kernels)
  therefore pass on every platform, while an error of the port still shows
  up in the cases with little amplification.
* The Fortran outputs of cases that together cover every branch are kept as
  constants in ``src/core/polradtran/rt4/test``, and the port is compared
  with them in the same way.
* Evans' benchmarks:

  * ``runtestc``: given the optics that Evans' program reads for this script
    (cirrus of horizontally oriented ice columns at 340 GHz, from the DDA
    scattering file), ARTS's RT4 reproduces his table.
  * ``runtestr``: a rain layer over water is solved from its Mie Legendre
    series by the Python benchmark test in ``tests/core/disort``.

Build
-----

RT4 is C++ and is always built.  Its MIT licence is registered as bundled
code (see :doc:`dev.licenses`), so it is also allowed in LGPL builds.

Tests
-----

Every reference is derived outside the code under test:

* closed forms derived in the tests:

  * quadrature exactness;
  * an isothermal blackbody;
  * gas-only layers;
  * doubling convergence;
  * a layout test built to detect transpositions, with lower-triangular,
    stream-dependent K, Stokes-asymmetric absorption and a forward-only
    phase matrix;
  * every kind of ground;
  * isothermal Kirchhoff with Rayleigh scattering and a non-reciprocal
    variant, also beyond the Fortran's array sizes;
  * the error paths.
* the Fortran outputs described above;
* Evans' benchmark tables.

The C++ tests are in ``src/core/polradtran/rt4/test`` and
``src/tests/polradtran``.  The Python tests are in ``tests/core/rt4``: they
repeat the closed forms through the bindings and are the maintained usage
examples.  The comparison with VDISORT is in
``src/core/disort-cpp/test/vdisort``, whose ``COVERAGE.md`` lists its cases;
see `Mapping to VDISORT inputs`_.

Interface
---------

A problem (``rt4.h``) describes:

* the streams: the number of quadrature nodes per hemisphere, the quadrature
  rule, and optional zero-weight extra angles;
* the layers, top-down: the interface heights, the interface temperatures,
  the scalar gas extinction of each layer, and for each layer either an
  optics set or none (gas only);
* the optics sets: the extinction matrix, the absorption vector and the
  azimuthally averaged phase matrix on the streams (see `Conventions`_);
* the boundaries: the sky temperature, and the ground with its temperature;
* the numerics: the largest optical thickness of the initial doubling
  sublayer, and the tolerance of the energy-conservation check;
* the frequency.

The result holds the streams with their weights, and the upward and downward
radiance at every level and stream.  The Python classes take every field as
a keyword argument with the C++ default, and they accept numpy arrays and
lists.  The ground is returned as a copy, so assign a new one to change it.
The solve releases the GIL.

Conventions
-----------

**Stokes basis.**

* [I, Q], with the meridional plane as the reference: "vertical"
  polarization lies in the plane of the ray and the z axis.
* I = I_v + I_h and Q = I_v - I_h.
* The same basis is used in both hemispheres, so Q does not change sign
  between up and down.  A warm dielectric surface emits Q > 0 at oblique
  angles.
* U and V are not computed.

**Hemispheres.**

* The downward hemisphere (index 0) is radiation propagating downward,
  toward increasing optical depth: RT4's "+".
* The upward hemisphere (index 1) is radiation propagating upward: RT4's "-".
* The upward radiance at a level on stream mu is what a sensor there sees
  looking down at the nadir angle acos(mu).

**Streams.**

* A stream is ``mu = |cos(zenith)|`` in (0, 1], with the same values in both
  hemispheres.
* The quadrature nodes come first, ascending; none of the rules includes
  mu = 0.  The zero-weight extra angles follow, in the order given, each in
  (0, 1].
* The weights are for the integral over mu in [0, 1] and sum to 1.  The
  2 pi azimuth factor is not included.
* The rules:

  * double Gauss ('D'): an nmu-point Gauss-Legendre rule mapped to [0, 1];
  * Gauss ('G'): the positive half of a 2 nmu-point Gauss-Legendre rule on
    [-1, 1];
  * Lobatto ('L'): the positive half of a 2 nmu-point Lobatto rule; it
    includes mu = 1.
* The extra angles receive scattering and reflection but contribute nothing
  to the angular integrals.  The one exception is the discrete ground, whose
  columns are applied without weights.

**Layers and levels.**

* Both are ordered top-down.  Level 0 is the top of the atmosphere, and the
  last level is just above the ground.
* Only the absolute differences of the interface heights are used.  Any
  length unit works if the extinctions are given in its reciprocal.
* The interface temperatures must be positive.  The scalar gas extinction of
  a layer must be non-negative.

**Optics**, all per unit length:

* The extinction matrix K is given for propagation in each hemisphere and
  stream, and the absorption vector a likewise.  RT4 multiplies a by the
  Planck function.
* The gas extinction is scalar and unpolarized.  RT4 adds it to every
  Stokes diagonal of K and to the I component of a.
* The phase matrix is indexed by the outgoing hemisphere, the incoming
  hemisphere, the outgoing and incoming streams, and the outgoing and
  incoming Stokes components.  The four (out, in) hemisphere blocks are
  RT4's quadrants: (down, down) is 1 (+ <- +), (down, up) is 2 (+ <- -),
  (up, down) is 3 (- <- +) and (up, up) is 4 (- <- -).
* The phase matrix is the azimuthal mean ``(1 / 2 pi) int Z dDelta-phi`` of
  the phase matrix in the meridional basis, per unit length and per
  steradian.  It includes the number density, is not normalised to 4 pi,
  and contains no quadrature weights.
* The scattering source into (h, mu_i) is
  ``2 pi sum_j sum_h' w_j Z(h <- h')(i, j) I(h', mu_j)``.
* Energy conservation on the streams reads
  ``K11(h, mu_j) = a1(h, mu_j) + 2 pi sum_i w_i [Z(up <- h) + Z(down <- h)](1, i; 1, j)``.
  The solver checks it on the quadrature streams to the tolerance given in
  the problem, and rejects optics that break it; an infinite tolerance
  disables the check.  Nothing renormalises the phase matrices (ARTS 2.6
  rescaled them).
* For an isothermal Kirchhoff test (I = B and Q = 0 on the streams), use
  ``a_s(h, mu_j) = K_sI - 2 pi sum_{i, h'} w_i Z_sI(h <- h')(j, i)``.

**Planck source.**

* Within a layer, the Planck function B, not the temperature, is linear in
  optical depth: B is linear in height, and the layer is homogeneous.
* Gas-only layers are solved analytically, with ``exp(-tau / mu)`` and the
  exact linear-in-tau source integral.
* Scattering layers are doubled from a thin initial sublayer, so the source
  is piecewise constant per sublayer.  The error is first order in that
  sublayer's slant thickness, which is at most ``max_delta_tau / mu_min``.
* RT4 chooses the number of doublings from the I-I extinction of the first
  downward stream plus the gas extinction.

**Units.**

* The frequency is in Hz.  The radiances, including the emission of a
  discrete ground, are in W m-2 Hz-1 sr-1.
* The sky is an isotropic, unpolarized blackbody at the sky temperature.  The
  Planck function is 0 at 0 K, and negative temperatures are rejected.

**Grounds.**

* Lambertian (RT4 'L'): reflection ``2 A mu_j w_j`` into every stream, I to I
  only (depolarizing), and emission ``[(1 - A) B_s, 0]``.
* Fresnel (RT4 'F'): the medium above has index 1, and the reflection is
  specular, stream by stream.  The reflection is
  ``R = [[R1, R2], [R2, R1]]``, with ``R1 = (|r_v|^2 + |r_h|^2) / 2`` and
  ``R2 = (|r_v|^2 - |r_h|^2) / 2``.  The emission is
  ``[(1 - R1) B_s, -R2 B_s]``.  The sign of ``Im n`` does not matter for
  [I, Q].
* Specular (RT4 'S'): a fixed ``R(out, in)`` of nstokes x nstokes, applied
  specularly to every stream.  The emission is
  ``[(1 - R(I, I)) B_s, -R(Q, I) B_s]``.
* Discrete (RT4 'A'): ``I_up(i) = sum_j R(i, j) I_down(j) + e(i)``.  R must
  include all quadrature factors; a Lambertian ground is ``2 A mu_j w_j``.
  Its columns of the extra angles are applied without weights, so set them
  to zero to keep the extra angles as pure outputs.  The emission e is given
  directly, and the ground temperature is not used.

Inputs from ARTS data
---------------------

``rt4_arts.h`` builds RT4 inputs from ARTS-native data.  It adds no physics.

**The optics at an atmospheric point** are ARTS's bulk scattering
properties, in the laboratory frame:

* They are the azimuthal Fourier modes that the scattering species give at
  the propagation zenith angles of the streams.  RT4 takes the mode 0.
* RT4's downward stream mu propagates at the ARTS propagation zenith angle
  ``180 - acos(mu)`` deg, and the upward stream at ``acos(mu)``.  ARTS's
  angular grids ascend strictly, so the stream cosines must be distinct.
* Azimuthally randomly oriented (ARO) species work as well as totally
  randomly oriented (TRO) ones.  RT4 is the only one of RT4, RT3 and VDISORT
  that can take ARO particles.
* The optics are as accurate as the species' modes:

  * exact for the closed-form scatterers (Rayleigh gas scattering,
    Henyey-Greenstein) and for Legendre series;
  * gridded ARO data give their modes only at their own zenith angles, with
    no interpolation in zenith;
  * gridded TRO data must first be converted to a Legendre series.
* ARTS's laboratory-frame phase matrix and RT4 use the same meridional
  basis and the same Q = I_v - I_h.  The [I, Q] block of the azimuthal mean
  does not depend on the sense of the azimuth.
* A vertical ray (mu = 1) has the meridional plane of its azimuth label,
  the limit along its meridian, so vertical rays need no special treatment.
* The solver then requires of these optics what it requires of any optics:

  * mirror symmetry, so species whose orientation is not symmetric under
    reflection in the horizontal plane are refused;
  * energy conservation on the streams, so forward-peaked optics whose peak
    falls between the streams need more streams.

**A problem from a propagation path** follows the conventions of the DISORT
workspace methods:

* The ray path, the atmospheric path and the gas propagation matrices have
  one entry per level, top first, with strictly decreasing altitudes.  The
  altitudes are the heights in metres, and only they are used of the ray
  path.
* The level temperatures come from the atmospheric path.
* The gas propagation matrices must be unpolarized.  The gas extinction of a
  layer is the mean of A at its two levels.
* The frequency is selected from the frequency grid by index.
* A layer's optics are the mean of its two levels' optics.  A layer is
  gas-only when that mean is zero.
* The solver settings, the ground and the boundary temperatures are passed
  through unchanged.
* The heights, temperatures and gas extinction come from a path-layer
  builder that RT3 uses as well.

**Tests:**

* In ``src/core/polradtran/rt4/test``:

  * the Rayleigh closed form of the m = 0 [I, Q] block, in all four
    quadrants, with extra angles and mu = 1;
  * Henyey-Greenstein against the addition theorem;
  * the path builder against its inputs and the level optics;
  * the error paths.
* The three solvers on one ARTS atmosphere, through their path builders, in
  the VDISORT test directory.
* The Python tests in ``tests/core/disort``, which run the builders through
  ``pyarts3``.  Run without ``ARTS_HEADLESS``, they overlay the per-type
  plots of the three solutions.

Limitations
-----------

* **Scope.**  Only [I] or [I, Q].  RT4 is thermal-only (no beam source),
  plane-parallel and limited to azimuthally symmetric media, and it returns
  only the m = 0 azimuthal mode.  With one Stokes component, the Q coupling
  of a polarizing ground or medium is dropped.
* **Mirror symmetry.**  RT4's doubling is hard-coded as symmetric: it
  discards the "-" reflection and transmission after the first step.  The
  medium must therefore be mirror symmetric between the hemispheres:

  * equal extinction in both hemispheres;
  * equal absorption in both hemispheres;
  * equal same-hemisphere phase blocks;
  * equal opposite-hemisphere phase blocks.

  The solver rejects optics that break this beyond rounding.  Because the two
  opposite-hemisphere blocks must be equal, no test can detect an exchange
  of them.
* **Energy conservation** is enforced, not repaired (see `Conventions`_).
* **Lambertian ground with Gauss or Lobatto streams.**  The Lambertian ground
  conserves energy on the streams only with double Gauss, where
  ``2 sum mu w = 1``.
* **No array limits.**  The Fortran's static arrays limited the number of
  streams times Stokes components, the number of layers and their product.
  The port sizes every array to the problem.
* **Singular inverse.**  The inverse of ``1 - R R`` throws for an exactly
  singular matrix, which is not expected for physical inputs.
* **No beam, no m > 0 modes, no U and V.**  RT3 (:doc:`dev.rt3`) is the
  reference for those paths.

Mapping to VDISORT inputs
-------------------------

Feed both solvers the same layer optics, so that a comparison tests the
solvers and not two preprocessing chains.  The comparison test in
``src/core/disort-cpp/test/vdisort`` implements this mapping, and it is the
reference for the code.

VDISORT can represent only a subset of RT4's inputs:

* an extinction that is a scalar times the identity, the same for all
  streams;
* an absorption vector ``[a1, 0]``, the same for all streams;
* scatterers whose m = 0 [I, Q] block does not couple to [U, V].

It cannot represent a direction-dependent or dichroic K(mu), or a
direction-dependent a(mu).  Such cases must be rejected, not approximated.

For a layer with particle extinction ``k_e``, absorption ``a1``, scattering
``k_s = k_e - a1`` (for optics that conserve energy) and gas extinction
``k_g``:

* **Quadrature.**

  * Use RT4's double Gauss with half as many nodes as VDISORT has streams.
    The nodes and weights are then VDISORT's, to rounding.
  * VDISORT orders its upward streams first, then the downward ones with the
    same cosines.
* **Extra angles.**

  * RT4's extra angles are VDISORT's user angles, ``+mu`` upward and ``-mu``
    downward.  The user-angle phase matrices come from RT4's phase rows of
    the extra angles, which must hold the physical phase matrix (their
    columns carry no weight).
  * VDISORT's user-angle evaluation is part of its C++ interface
    (``vdisort.h``); the Python bindings evaluate on the streams.
  * Over a Fresnel ground, VDISORT reflects the downward radiance at ``-mu``
    into the upward user angle ``mu``, so an upward user angle needs its
    downward partner, and VDISORT refuses it without one.
  * The ground's emission at a user angle is interpolated from the streams,
    unless the caller gives the exact boundary radiance.  With the
    interpolation removed, the reflected radiance agrees with RT4 to RT4's
    doubling error.
* **Optical depth.**  The cumulative ``(k_e + k_g) dz`` at the layer
  bottoms.  RT4's level l is VDISORT's ``tau_l``, with ``tau_0 = 0``.
* **Single-scattering albedo.**  ``omega = k_s / (k_e + k_g)``, and 0 for
  gas-only layers.
* **Phase matrix.**

  * VDISORT's m = 0 [I, Q] block is ``4 pi Z / k_s``, with Z the RT4
    quadrant that matches the hemispheres of the VDISORT streams.
  * It contains no omega, no weights and no ``2 - delta_m0``.  VDISORT's
    ``(omega / 2) sum_j W_j P`` is then RT4's ``2 pi sum_j w_j Z`` per unit
    length.
  * The [U, V] block does not affect [I, Q] at m = 0.
* **Source.**

  * VDISORT emits ``(1 - omega) B`` per unit optical depth, which is RT4's
    ``(a1 + k_g) B`` per unit length.
  * The Planck function is linear in the global optical depth.  The slope is
    ``(B_bottom - B_top) / dtau_l``, and the offset is
    ``B_top - slope * tau_top``.
* **Top boundary.**  The isotropic, unpolarized ``B(T_sky)``.
* **Lambertian ground.**  VDISORT's Lambertian modes, whose discrete
  operator ``(2 A / pi) pi W mu`` is RT4's ``2 A mu w``.  The emission is
  ``[(1 - A) B_s, 0, 0, 0]``.
* **Fresnel ground.**  VDISORT's Fresnel modes, whose specular part gives
  RT4's R at the nodes.  The emission is ``[(1 - R1) B_s, -R2 B_s, 0, 0]``.
* **Discrete ground.**  A custom VDISORT BRDF that returns
  ``rho(mu_i, mu_j) = R(i, j) / (pi w_j mu_j)`` at the cosines it is given,
  because VDISORT adds ``pi W_j mu_j`` itself at m = 0.  Use the same
  emission for both solvers.
* **One Stokes component.**  Put only ``4 pi Z_II / k_s`` in VDISORT's M00
  element and compare I.
* **Units.**  Both solvers are linear in B.  Pass B per hertz to VDISORT to
  compare directly.
* **Stokes basis.**  Both solvers use [I, Q] with Q = I_v - I_h in the
  meridional plane.

RT4 has the opposite limits: no U and V, no m > 0 modes, no beam source, and
mirror symmetry between the hemispheres is required.

What the comparison shows:

* Gas-only atmospheres agree to rounding.
* Otherwise, RT4's first-order initial-sublayer error dominates.  It is
  proportional to ``max_delta_tau``, halving exactly when ``max_delta_tau``
  halves, and it does not grow with the number of streams.
* Richardson extrapolation of RT4 to ``max_delta_tau = 0`` removes that
  error, and the agreement then reaches the comparison's tolerance.
