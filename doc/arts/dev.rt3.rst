RT3 reference solver
====================

RT3 is Evans' polarized doubling-adding solver for randomly oriented
particles, from the polradtran package.  ARTS 3 keeps it as an external
reference for validating other solvers, in particular the paths of the
polarized discrete-ordinate solver VDISORT that RT4 (:doc:`dev.rt4`) cannot
check: the solar beam, the Fourier modes m > 0, and U and V.  RT3 solves the
plane-parallel radiative transfer equation with a direct (solar) beam and
thermal sources, for every Fourier azimuth mode and the Stokes components
[I], [I, Q], [I, Q, U] or [I, Q, U, V].  Its equations are in
:doc:`concept.rt3`.

The C++ is in ``src/core/polradtran/rt3``, and what RT3 shares with RT4 is
in ``src/core/polradtran``.  The Python bindings are the submodules
``pyarts3.arts.rt3`` and ``pyarts3.arts.polradtran``.  There are no workspace
methods, variables or agendas.  The headers and the docstrings of the
bindings are the reference for names, fields, defaults and array layouts.
This page describes the design and the conventions.

Provenance
----------

* K. F. Evans' polradtran is distributed from
  https://nit.coloradolinux.com/polrad.html under the MIT licence.  RT3 is
  described by :cite:t:`Evans1991`.  ARTS 2 never shipped RT3.
  ``3rdparty/polradtran`` keeps Evans' licence, his ``README`` and his
  benchmark scripts, unchanged; its ``README.ARTS`` gives their origin.
* ARTS keeps no Fortran: RT3 is a C++ port of Evans' ``RADTRAN`` and of the
  routines it calls.
* As for RT4, the scattering properties are passed in memory: scattering sets
  (extinction, scattering coefficient and Legendre series), with a set index
  per layer, in place of Evans' scattering files.  Each set is processed
  once.  With delta-M, the coefficients above the end of a series are zero.

Design of the port
------------------

The port follows the same paradigms as RT4's (see :doc:`dev.rt4`):

* Evans' routines are ported one to one, as functions named after them in
  lower case, in files named after Evans' files.
* What RT3 and RT4 share exists once, in the ``polradtran`` namespace.  This
  covers the doubling, the adding, the internal radiance, the non-scattering
  layer, the sky radiance, the Lambertian and Fresnel grounds, the external
  surface layer, the level loop, the initial sublayer and the work data.
  RT4 uses the shared code in the azimuthal mode 0, with the thermal source
  alone.
* The quadratures, the Planck function (SI, per hertz) and the Fresnel
  amplitudes are ARTS's, and the matrix algebra is BLAS and LAPACK.
* The ground is data: for every azimuth mode, the solver core takes the
  ground's reflection, its own radiance, and the radiance it reflects from
  the direct beam per unit of direct flux.  The direct flux that reaches the
  ground is computed inside the core, which adds that part itself.
* All scratch arrays live in RT3's work data, which extends RT4's.  The
  solver sizes it, so repeated calls with the same sizes allocate nothing.
  Its values between calls are unspecified, except the table of the FFT.  No
  other state is kept, so concurrent solves are safe.

What is particular to RT3:

* **Scattering.**  A Legendre series of the six scattering-plane elements is
  a vector of rtepack's compact planar Mueller matrices, with the elements
  named, not in the column order of RT3's files.

  * The phase matrices in the meridional planes are full Mueller matrices,
    of which only the leading nstokes x nstokes block is transformed.
  * The Stokes rotations are rtepack's closed-form rotation of the compact
    matrix.  Evans' symmetry operation on the off-diagonal 2 x 2 blocks is
    rtepack's mirror.
  * The Legendre polynomials are those of ARTS's Legendre module, made once
    for all six series.  The cosine of the scattering angle is clamped to
    [-1, 1].
  * Every scattering set is precomputed for every azimuth mode.
* **Delta-M.**  A division by zero in the delta-M scaling (a non-positive
  extinction, ``1 - f = 0`` or ``1 - albedo f = 0``) throws, where the
  Fortran returned NaN or infinity.
* **Normalisation check.**  RT3's discrete normalisation check throws
  instead of stopping, and it also catches NaN.
* **FFT.**  Evans' real FFT is self-contained, in its own source file, and is
  used by nothing outside the scattering routines.  It needs a power of two
  of values but has no other size limit.  Its phase table grows as needed
  and lives in the work data.  Another FFT can replace it, if it gives the
  same packed format and scaling, which the header documents.

Validation of the port
----------------------

* The port was compared with the Fortran, routine by routine and as a whole,
  on random inputs that cover every branch: 1 to 4 Stokes components, every
  quadrature with extra angles, azimuth orders, delta-M, every source, both
  grounds, gas, thin, thick and shared layers, and each of RT3's Legendre
  summation cases.
* The tolerance paradigm is RT4's: a few units of rounding, times the
  amplification of the computation.  For the doubling, that is :math:`2^n`
  times the condition number of ``1 - R R``.  Any replaced numerics add a
  separate term.
* The Fortran outputs of the branch-covering cases are kept as constants in
  ``src/core/polradtran/rt3/test``, and the port is compared with them in
  the same way.  The test also runs every case with one work data shared over
  all cases, and with one that starts out as NaN, and requires identical
  bits.
* The diagonal of the initial sublayer's transmission keeps the Fortran's
  form, ``1 - f (1 - albedo P)``, rounded once.  The sublayer's optical
  thickness f is as small as ``max_delta_tau``, so the last bit of the
  diagonal is a relative ``eps / max_delta_tau`` of f, and the doubling
  carries it to the radiances.  Another rounding order would move the
  radiances by that much.

Build
-----

RT3 is C++ and is always built, like RT4.  Its arrays are sized to the
problem, so the large static buffers of the Fortran are gone.

Tests
-----

Every reference is either Evans' own benchmark output or a closed form
derived in the test.  None is an output of this build.  The C++ tests are in
``src/core/polradtran/rt3/test``.  The Python tests in ``tests/core/rt3``
repeat several of them through the bindings and serve as the usage examples.

* **Quadrature exactness** of the three rules.
* **Evans' Mie benchmark** (``runmietest``, the Mie case of
  :cite:t:`Evans1991`).  It has a Mie layer of optical depth 1 and
  single-scattering albedo 0.99 over a Lambertian ground, under a grazing
  sun.  It uses Gauss streams and every Stokes component, and gives I, Q, U,
  V and the fluxes at both levels, at 0, 90 and 180 deg.  The Python test
  reads the scattering file and the expected output from the script itself.
* **Evans' general benchmark** (``runtesta``): a Rayleigh layer over a Mie
  layer with gas absorption, solar and thermal sources at 3 um, and a
  Lambertian ground.
* **Gas-only atmospheres**, with the thermal source and 1 to 4 Stokes
  components, at microwave and infrared frequencies.  They run over a
  Lambertian ground with double Gauss streams and over a Fresnel ground with
  Gauss streams and extra angles.  The reference is the exact solution with
  the Planck function linear in optical depth, and every m > 0 mode and U,
  V must be exactly 0.
* **Single scattering** by a thin conservative Rayleigh layer over a black
  ground, lit by the beam, with 1 to 4 Stokes components.  The reference
  Stokes column is built from the vector geometry of dipole scattering,
  independently of RT3's rotation formulas.  It fixes RT3's azimuth origin
  and sense, Q, and the sign of U.  The neglected multiple scattering must
  shrink with the optical depth.
* **Delta-M with f = 0** is the identity.
* **The direct beam** with the delta-M scaled optical depth.
* **Extra angles** with Lobatto and double Gauss leave the radiances on the
  nodes unchanged, and an extra angle on a node reproduces it.
* **Isothermal Kirchhoff**: I = B and Q = U = V = 0 in an isothermal
  Rayleigh medium over a Fresnel ground, also beyond the Fortran's array
  sizes.
* **The error paths.**

Evans' tables were printed by his ``OUTPUT_FILE``, which sums the Fourier
series in single precision and prints six significant digits.  The tolerance
of each entry is therefore half a unit in its last printed digit, plus a
bound on that single-precision summation error.  The ``runtesta`` table was
made with Evans' five-digit Planck constants.  RT3 evaluates the Planck
function only at the interface, ground and sky temperatures, so the test
reproduces the five-digit Planck values exactly: it passes the temperatures
at which ARTS's Planck function equals the five-digit one.

The comparison of VDISORT with RT3 is in
``src/core/disort-cpp/test/vdisort``, whose ``COVERAGE.md`` lists its cases.
Its part E solves the problems of Evans' two scripts, read from the
scripts, with VDISORT.  See `Mapping to VDISORT inputs`_.

Interface
---------

A problem (``rt3.h``) describes:

* the streams: the number of quadrature nodes per hemisphere, the quadrature
  rule, and optional zero-weight extra angles, with any rule;
* the highest Fourier azimuth mode;
* the sources: the direct beam (its flux on the horizontal plane at the top,
  and its zenith cosine), and whether the layers and a Lambertian ground
  emit thermally;
* delta-M scaling on or off;
* the layers, top-down: the interface heights, the interface temperatures,
  the scalar gas extinction of each layer, and for each layer a scattering
  set or none (gas only);
* the scattering sets: the particle extinction and scattering coefficients,
  and the Legendre series of the scattering-plane matrix;
* the boundaries: the sky temperature, and the ground with its temperature;
* the largest optical thickness of the initial doubling sublayer;
* the frequency.

The result holds the streams with their weights; the Fourier coefficients of
the upward and downward radiance at every level, mode and stream; and the
upward and downward fluxes.  A function in ``rt3.h`` sums the coefficients
at given azimuths, as Evans' ``OUTPUT_FILE`` does, in double precision.  The
Python classes take every field as a keyword argument with the C++ default.
The ground is returned as a copy, and the solve releases the GIL.

Conventions
-----------

**Geometry and azimuth.**

* A right-handed frame with z pointing up, out of the medium.  A ray
  propagating at zenith cosine mu_z (> 0 upward) and azimuth phi has the
  direction ``k = (sin(theta) cos(phi), sin(theta) sin(phi), mu_z)``.
* The direct beam propagates downward toward phi = 0.  So phi is the
  azimuth of the propagation direction of a ray, counted from x toward y,
  relative to the propagation direction of the beam.  For upwelling
  radiation, phi = 0 is forward scattering in azimuth.
* The single-scattering test pins this frame, together with the sign of U,
  against the vector construction.  It is the frame of VDISORT's polarized
  tests, with the beam azimuth phi0 = 0.

**Stokes basis.**

* [I, Q, U, V], with the meridional plane as the reference.
  ``h = k x z / |k x z|`` is horizontal, and ``v = h x k`` lies in the plane
  of k and z.
* I = I_v + I_h, Q = I_v - I_h and U = 2 Re(E_v E_h*).  The same basis is
  used in both hemispheres.  A warm dielectric surface emits Q > 0 at
  oblique angles.
* RT3 transports V through F34 (and the Fresnel R4), so V has the sign
  convention of the scattering data.  Evans' Mie series is in the convention
  of :cite:t:`Evans1991`, in which RT3 reproduced the V of Garcia and
  Siewert (1989).  For the same spheres, ARTS's Mie code gives F34 of the
  opposite sign (see `Inputs from ARTS data`_), so RT3 run on ARTS data gives
  -1 times V in the convention of that paper.  Which of the two conventions
  makes V positive for left-hand circular polarization is not tested.

**Streams and hemispheres.**

* A stream is ``mu = |cos(zenith)|`` in (0, 1], with the same values in both
  hemispheres.  The quadrature nodes come first, ascending, followed by the
  zero-weight extra angles in the order given.
* The rules are those of RT4 (see :doc:`dev.rt4`): Gauss, double Gauss and
  Lobatto.  The weights are for the integral over [0, 1] and sum to 1.
* The extra angles are allowed with every rule; Evans' RT3 had them with the
  Gauss rule only.  They receive scattering and reflection but carry no
  weight, so they do not change the radiances on the quadrature nodes.
* The downward hemisphere is radiation propagating downward, toward
  increasing optical depth: RT3's "+", printed with mu > 0 by Evans'
  program.  The upward hemisphere is RT3's "-".
* :doc:`concept.rt3` labels the hemispheres as DISORT does, "+" upward, and
  names the operators of a slab by their outgoing and incident hemispheres.
  The code's (Evans') R, T and S of the "+" part are its :math:`R^{-+}`,
  :math:`T^{--}` and :math:`S^-`.  Those of the "-" part are its
  :math:`R^{+-}`, :math:`T^{++}` and :math:`S^+`.  The same holds for RT4.

**Fourier series.**

* The radiance at azimuth phi is ``sum_m c_m cos(m phi)`` for I and Q, and
  ``sum_m c_m sin(m phi)`` for U and V, for m = 0 up to the highest mode.
  These are RT3's cosine modes of I and Q and sine modes of U and V.  The
  m = 0 coefficients of U and V are 0.
* At phi = 90 deg only the odd modes of U and V contribute, and at 0 and 180
  deg U and V vanish.  Evans' tables at 0, 90 and 180 deg therefore do not
  test the even modes of U and V; the single-scattering test does.
* With only thermal sources, every m > 0 mode and U and V are exactly 0.

**Layers, levels and units.**

* Layers and levels are top-down.  Level 0 is the top, and the last level is
  just above the ground.  Only the differences of the interface heights are
  used; the extinctions are in the reciprocal unit.
* The frequency is in Hz.  The radiances are in W m-2 Hz-1 sr-1, the fluxes
  and the direct flux in W m-2 Hz-1.  The problem is linear in its sources,
  so a solar-only result does not depend on the frequency.
* The fluxes are ``2 pi sum_i w_i mu_i c_0``; only I and Q can be non-zero.
  The downward flux includes the direct beam ``F exp(-tau / mu0)`` in I,
  with the delta-M scaled optical depth when delta-M is on.  The downward
  radiance does not include the direct beam.

**Scattering sets.**

* The extinction and scattering coefficients are those of the particles,
  per unit length.  The scalar, unpolarized gas extinction is added per
  layer, and the single-scattering albedo of a layer is
  ``scattering / (extinction + gas)``.
* The Legendre series is that of the scattering-plane phase matrix

  ::

    [[F11, F12,   0,   0],
     [F12, F22,   0,   0],
     [  0,   0, F33, F34],
     [  0,   0, -F34, F44]]

  as a plain Legendre series in cos(Theta),
  ``F = sum_l c_l P_l(cos(Theta))``.  The basis is that of the scattering
  plane, with Q = I_par - I_perp, so Rayleigh scattering has
  F12 = -3/4 sin^2(Theta).
* The coefficients include the factor 2 l + 1, so Henyey-Greenstein has
  ``(2 l + 1) g^l``.  The phase function is normalised to 1 over 4 pi: the
  F11 coefficient of degree 0 must be 1.
* RT3 sums only F11 for one Stokes component.  Trailing all-zero
  coefficients are dropped before the solve.
* With delta-M, every set is scaled with M = 2 times the number of streams,
  including the extra angles, because RT3 passes its total stream count:

  * ``f`` is the F11 coefficient of degree M divided by 2 M + 1;
  * the extinction becomes ``(1 - omega f) k``;
  * the albedo becomes ``(1 - f) omega / (1 - omega f)``;
  * the series becomes ``(2 l + 1) (c_l / (2 l + 1) - f id) / (1 - f)``,
    with ``id`` the identity, to degree M - 1, so F12 and F34 become
    ``c_l / (1 - f)``.

  The beam is attenuated with the scaled extinction.  There is no correction
  of the radiances for the truncated peak.

**Sources.**

* The direct beam is on when its flux is positive, and it requires a
  Lambertian ground.
* The thermal switch turns the emission of the layers on or off:
  ``(1 - omega) B``, with B linear in optical depth.  It does the same for
  the emission of a Lambertian ground, ``(1 - A) B``.
* The sky and a Fresnel ground always emit, as in RT3; set their
  temperatures to 0 K to remove them.
* The Planck function is 0 at 0 K, and negative temperatures are rejected.
  With thermal emission, every interface temperature must be positive.

**Grounds.**

* Lambertian (RT3 'L'): reflection ``2 A mu_j w_j`` of the m = 0 mode into
  every stream, I to I only, and direct-beam reflection
  ``A F_direct(ground) / pi``.  The diffuse reflection conserves energy on
  the streams only with double Gauss (``2 sum mu w = 1``).
* Fresnel (RT3 'F'): specular reflection in every mode, under a medium of
  index 1.  The reflection is ``[[R1, R2, 0, 0], [R2, R1, 0, 0], [0, 0, R3,
  -R4], [0, 0, R4, R3]]``, with ``R1 = (|r_v|^2 + |r_h|^2) / 2``,
  ``R2 = (|r_v|^2 - |r_h|^2) / 2``, ``R3 = Re(r_v r_h*)`` and
  ``R4 = Im(r_v r_h*)``.  The emission is ``[(1 - R1) B, -R2 B, 0, 0]``.
  Without a beam the field is azimuthally symmetric, so R3 and R4 never
  act.

Inputs from ARTS data
---------------------

``rt3_arts.h`` builds RT3 inputs from ARTS-native data.

**The scattering set at an atmospheric point** comes from the Legendre
series that the scattering species give themselves (their TRO spectral bulk
properties):

* The species' coefficients on the orthonormal Y_l0 are converted to plain
  Legendre coefficients, ``c_l = a_l sqrt((2 l + 1) / 4 pi)``.
* The elements are ARTS's, by name, without any sign change.  The series is
  normalised so that its F11 coefficient of degree 0 is 1.
* The extinction is K11, and the scattering coefficient is K11 - a1, per
  metre.
* The phase-function integral ``4 pi c_0(F11)`` must equal K11 - a1, to a
  tolerance relative to K11.  Otherwise it is an error: RT3 normalises the
  series and takes the albedo from the scattering coefficient, so a mismatch
  would change the scattered energy silently.  An infinite tolerance
  disables the check.
* Without particles, the set has zero extinction and scattering and an
  isotropic, depolarizing series.
* Only species that can give a Legendre series can be used.  Gridded TRO
  data must first be converted to a Legendre series.  ARO data cannot be
  used, because RT3 takes randomly oriented particles only.

**Sign conventions.**  The mapping without sign changes is the one under
which RT3 transports ARTS's Stokes vector:

* ARTS's laboratory-frame phase matrix, for the propagation directions
  (za, aa), equals the vector-geometry phase matrix of the same F in RT3's
  meridional basis, element by element and F34 included.  The angles map as
  ``mu = cos(za)`` and ``phi = -aa``: ARTS's azimuth runs clockwise seen
  from above, RT3's counterclockwise.  This is tested against vector
  geometry in VDISORT's ARTS-input test.
* Against an external reference, the ARTS-input test of RT3 reproduces
  Evans' runmietest series from ARTS's Mie code.  These are spheres of
  refractive index 1.44 at 0.951 um, with a gamma distribution of effective
  radius 0.2 um and effective variance 0.07.  F11, F12 and F33 (and
  F22 = F11, F44 = F33) agree with Evans' table, and F34 agrees with the
  opposite sign.
* So F12, the meridional U and the relative sign of F34 are pinned, and RT3
  run on ARTS data gives V with the opposite sign to that convention.

**A problem from a propagation path** follows the path conventions of RT4's
builder (see :doc:`dev.rt4`), and uses the same path-layer builder:

* A layer's scattering set has the mean extinction and scattering of its two
  levels, and the scattering-weighted mean of their normalised series, which
  is the normalised series of the mean phase matrix.
* A layer with zero mean extinction is gas-only.
* By default, the Legendre degree is RT3's limit for the quadrature (see
  `Limitations`_).  With delta-M it is raised to 2 times the number of
  streams, so that RT3 can read the delta-M fraction.
* The problem has thermal emission and no beam.  Set the beam and the
  thermal switch on the problem for other sources.

**Tests:**

* In ``src/core/polradtran/rt3/test``:

  * ARTS's Rayleigh gas scattering against the closed-form series of
    Evans' Rayleigh scattering file;
  * Evans' runmietest series, as above;
  * the path builder against its inputs, and the scattering-weighted layer
    mean of Rayleigh and Henyey-Greenstein scattering whose mix changes with
    height;
  * the error paths.
* The three solvers on one ARTS atmosphere, in the VDISORT test directory.
  RT3 and VDISORT get the same inputs by different routes, and they agree to
  their doubling error, which is proportional to ``max_delta_tau``.
* The Python tests in ``tests/core/disort``, which run the builders through
  ``pyarts3`` and overlay the per-type plots of the three solutions.

Limitations
-----------

* **Legendre truncation is an error.**  RT3 silently truncates a series,
  after delta-M scaling, to its highest kept degree:

  * Gauss: ``4 nmu - 3``;
  * double Gauss: ``2 nmu - 3``;
  * Lobatto: ``4 nmu - 5``;
  * at least 1, with nmu the quadrature nodes only.

  The solver raises an error instead when that would drop a non-zero
  coefficient, and it drops zeros silently.  Delta-M gives the degree
  2 nmu_total - 1, so delta-M is rejected, unless the scaled series vanishes
  there, with:

  * double Gauss (nmu >= 2);
  * Gauss with nmu or more extra angles.
* **Normalisation.**  The F11 coefficient of degree 0 must be 1, after
  delta-M scaling, to a tight tolerance.  For a series within RT3's degree
  limit, RT3's discrete normalisation equals it up to rounding.
* **No array limits.**  Every array is sized to the problem.  Every
  scattering set is precomputed, including an unused one, so memory grows as
  the number of sets times the number of modes times ``2 N^2``, with N the
  number of streams times Stokes components.
* **Accuracy.**  The initial doubling sublayer is first order in its slant
  thickness ``max_delta_tau / mu_min``.  Its transmission is stored as
  ``1 - O(max_delta_tau)``, so the results carry a rounding floor of about
  ``eps / max_delta_tau``.  The thermal source of a scattering layer is
  linear in optical depth.
* **Cost.**  For each output level, RT3 adds all layers above and below anew,
  so the cost grows as the square of the number of layers, per azimuth mode.
  Every level is returned.
* **Singular inverse.**  The inverse of ``1 - R R`` throws for an exactly
  singular matrix, which is not expected for physical inputs.
* **Grounds.**  Only a Lambertian ground with a beam, and no BRDF.
* **Azimuth sampling.**  RT3 samples the azimuth for its Fourier modes as
  densely as its Legendre degree requires.

  * That is exact when the scattering matrix is regular at forward and
    backward scattering, as a physical matrix is: F12 and F34 vanish there,
    F22 + F33 has a double zero at 180 deg, and F22 - F33 has one at 0 deg.
  * For a matrix that is not regular there, the laboratory-frame matrix is
    not a trigonometric polynomial in the azimuth.  RT3's modes then carry
    an aliasing error that does not fall with ``max_delta_tau``.
  * VDISORT and RT4 take the species' own Fourier modes, and do not sample.

RT3 complements RT4: it has the beam, the m > 0 modes and U and V.  But it
takes only randomly oriented particles (a scattering-plane matrix with six
elements), a scalar unpolarized gas extinction, and no dichroic extinction.

Mapping to VDISORT inputs
-------------------------

The comparison test in ``src/core/disort-cpp/test/vdisort`` implements this
mapping, and it is the reference for the code.  Its cases are listed in that
directory's ``COVERAGE.md``.  Feed both solvers the same scattering sets, so
that a comparison tests the solvers and not two preprocessing chains.  Give
VDISORT as many Fourier modes as RT3 has.  For a layer with the set
(extinction k, scattering sigma, series) and gas extinction ``k_g``:

* **Quadrature.**

  * Use RT3's double Gauss with half as many nodes as VDISORT has streams.
    The nodes and weights are then VDISORT's, to rounding.  RT3's double
    Gauss rule keeps Legendre degrees up to ``2 nmu - 3`` only.
  * VDISORT's upward stream i is RT3's upward stream i.  VDISORT's downward
    streams follow its upward ones, in the same order.
  * RT3 takes extra angles with every rule, so its extra angles can match
    VDISORT's user angles on identical quadrature streams.  The comparison
    test compares them at seven angles between the streams, at every level
    and azimuth, for the Mie and the multilayer cases.
* **Optical depth.**  The cumulative ``(k + k_g) dz`` at the layer bottoms.
  RT3's level l is VDISORT's ``tau_l``, with ``tau_0 = 0``.
* **Single-scattering albedo.**  ``omega = sigma / (k + k_g)``, and 0 for
  gas-only layers.
* **Phase matrix.**

  * Build the ordinary Fourier coefficients on the signed streams (> 0
    upward), without ``epsilon_m``:
    ``C^m(mu_o, mu_i) = (1 / 2 pi) int Z(mu_o, 0; mu_i, phi') cos(m phi') dphi'``,
    and ``S^m`` with ``sin(m phi')``.
  * ``Z = L_out^T F(Theta) L_in`` is the laboratory-frame phase matrix of the
    series in the meridional basis of `Conventions`_, from ARTS's
    polarization geometry (:ref:`Sec Polarization scattering`), which its
    own test checks against electric fields.  The test's ``lab-frame.h``
    applies it.
  * Combine the coefficients into VDISORT's diffuse operator without omega,
    weights or ``epsilon_m``.
  * For an identical discrete problem, sample at RT3's own azimuths.  For a
    series that is regular at Theta = 0 and 180 deg, the samples do not
    matter.
  * With fewer than four Stokes components, zero the rows and columns of Z
    from nstokes on before the transform.  RT3 transports the leading block,
    and VDISORT's other components then stay exactly 0.
* **Beam.**

  * VDISORT's beam Stokes vector is ``[F / mu0, 0, 0, 0]``, the irradiance
    normal to the beam, with ``mu0`` the beam cosine, which must not be a
    node.
  * The beam's phase matrix comes from ``C^m(mu_i, -mu0)`` and
    ``S^m(mu_i, -mu0)``.  The beam is a delta in azimuth, so it has cosine
    terms only.
  * The cosine system (I^c, Q^c, U^s, V^s) gets rows I and Q of ``C^m`` and
    rows U and V of ``S^m``.  The sine system (I^s, Q^s, U^c, V^c) gets rows
    I and Q of ``S^m`` and rows U and V of ``C^m``.  VDISORT applies
    ``epsilon_m`` itself.
  * VDISORT's beam combination does exactly this.  The diffuse combination
    does not apply to the beam (see ``COVERAGE.md``).
* **Azimuth.**  VDISORT's radiance at ``phi0 + psi`` is RT3's at ``psi``.
  Both are azimuths of the propagation direction: VDISORT's beam propagates
  toward ``phi0``, and RT3's toward 0.
* **Sources.**

  * The Planck function is linear in the global optical depth.  The slope is
    ``(B_bottom - B_top) / dtau_l``, and the offset is
    ``B_top - slope * tau_top``.  VDISORT emits ``(1 - omega) B``, as RT3
    does.
  * The top boundary is ``B(T_sky)``, always, because RT3 always includes
    the sky.
  * With thermal emission, the ground emits ``[(1 - A) B_s, 0, 0, 0]``.
* **Lambertian ground.**  VDISORT's Lambertian modes.  Their diffuse operator
  equals RT3's ``2 A mu w``, and their beam reflection equals RT3's
  ``A F_direct / pi``.
* **Fluxes.**  VDISORT's upward flux is RT3's upward I flux.  Its downward
  diffuse plus direct flux is RT3's downward I flux.  The Q fluxes are
  ``2 pi sum_i W_i mu_i`` times the m = 0 Q radiance over each hemisphere.
* **Delta-M.**

  * RT3's delta-M corresponds to VDISORT given the scaled set:
    ``k' = (1 - w f) k`` with ``w = sigma / k``, ``sigma' = (1 - f) sigma``,
    and the scaled series (see `Conventions`_).  The direct beam is then
    attenuated with the scaled optical depth in both.
  * This is not VDISORT's own delta-M with IMS and TMS corrections.
  * With double Gauss, the scaled series must vanish at the degrees
    ``2 nmu - 2`` and ``2 nmu - 1``.  The test builds such a series from Mie
    scattering plus a forward peak truncated at that degree.
* **Other grounds.**  RT3 allows only a Lambertian ground with a beam, so
  VDISORT's Fresnel and BRDF beam paths are not compared.
* **V.**  VDISORT reproduces RT3's V, which is generated through F34.  That
  shows that both use F34 consistently.  The relative sign of ARTS's F34 is
  in `Inputs from ARTS data`_.
* **Units.**  Both solvers are linear in their sources.  The results are per
  hertz, so pass the Planck function per hertz and the per-hertz direct flux
  to VDISORT unchanged.
