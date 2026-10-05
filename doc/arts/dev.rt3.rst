RT3 reference solver
====================

RT3 is Evans' polarized doubling-adding solver for randomly oriented
particles, from the polradtran package.  ARTS 3 keeps it only as an external
reference for validating other solvers, in particular the solar-beam, m > 0
and U, V paths of the polarized discrete-ordinate solver VDISORT, which RT4
(:doc:`dev.rt4`) cannot check.  RT3 solves the plane-parallel radiative
transfer equation with a direct (solar) beam and thermal sources, for every
Fourier azimuth mode and the Stokes components [I], [I, Q], [I, Q, U] or
[I, Q, U, V].  It has a core C++ interface (``src/core/rt3``) and a
low-level Python interface (``pyarts3.arts.rt3``).  It has no workspace
methods, variables or agendas.

Provenance
----------

* **Original code.** K. F. Evans, polradtran (RT3/RT4), distributed from
  https://nit.coloradolinux.com/polrad.html under the MIT licence
  (``3rdparty/polradtran/LICENSE``).  RT3 is described in Evans and Stephens
  (1991), J. Quant. Spectrosc. Radiat. Transfer 46, 413-423.  ARTS 2 never
  shipped RT3.
  All of ``PolRadTran.tar`` (sha256 ``7b0eff79...a6f7cff9d``) is in the one
  folder ``3rdparty/polradtran``, as in the tar: where ARTS modifies a file,
  the modified file has the tar's name and the tar's file is kept beside it
  as ``.orig``; every other file (``README``, ``rt3.f``, ``rt4.f``,
  ``scatcnv.f``, the four test scripts and ``cl340d14.dda``) is unchanged.
  ``README.ARTS`` there lists all changes.
* **ARTS 3 changes**, each marked ``c ARTS3:`` in the source and listed in
  ``3rdparty/polradtran/README.ARTS``:

  * The scattering properties are passed in memory as scattering sets
    (extinction, scattering coefficient, Legendre series) with a set index
    per layer, instead of scattering files.  ``READ_SCAT_FILE`` is replaced
    by ``GET_SCAT_SET``, which copies a set and applies the delta-M scaling
    of ``READ_SCAT_FILE`` unchanged.  Each set is processed once.
  * With delta-M, RT3 scales the coefficients up to l = 2 NUMMU - 1 even
    when the series is shorter, and ``READ_SCAT_FILE`` read entries there
    that the file never set.  ``GET_SCAT_SET`` sets them to zero.
  * The routines and COMMON blocks that RT4 also defines (``CHECK_NORM``,
    ``COMBINE_LAYERS``, the three quadratures, ``PLANCK_FUNCTION``,
    ``/SCRATCH1/``, ``/SCRATCH2/`` and ten more) carry an ``RT3_`` prefix.
    Otherwise the two would share symbols and scratch memory in the ARTS
    link; ``librt3`` and ``librt4`` export no common symbol.
  * Planck constants computed from the exact SI h, c and k, as for RT4.  The
    original 5-digit constants give a bias of -3e-5 in the Rayleigh-Jeans
    limit, growing to -2e-4 and -3e-4 at 3 um and 300 and 200 K.
  * A new ``rt3_c_interface.f90`` with ``ISO_C_BINDING`` entry points.

  ``rt3.f``, the original main program, is built with the ``.orig`` files
  as the program ``rt3-evans``.  ``cpp.fast.polradtran-runmietest`` and
  ``cpp.fast.polradtran-runtesta`` run Evans' two RT3 scripts with it as
  they are (``src/tests/polradtran/polradtran-scripts.cpp`` executes their
  here-documents, so csh is not needed) and compare the output with his
  tables numerically: every value to one unit in its last printed digit
  (measured: 1 and 0), and the values that are zero by symmetry, REAL*4
  round-off of ``OUTPUT_FILE``, below 1e-7 of max I (measured: 5e-9).
  ``OUTPUT_FILE`` also shows how RT3 sums its Fourier series.  The RT4
  scripts are run likewise, see :doc:`dev.rt4`.

Build
-----

RT3 is built when CMake finds a Fortran compiler, independently of RT4, and
``-DENABLE_RT3=OFF`` turns it off.  As for RT4 (see :doc:`dev.rt4` for details):

* the Fortran language is enabled only after LAPACK has been found, so the
  BLAS/LAPACK choice does not depend on ``ENABLE_RT3``;
* the legacy ``.f`` sources are compiled with
  ``-std=legacy -fdefault-real-8 -fdefault-double-8`` (GNU) or ``-r8``
  (Intel), and default INTEGER is not promoted;
* configure with ``-DCMAKE_Fortran_FLAGS=""`` in a conda environment; the
  flags in use are printed at configure time;
* the library is built shared on macOS with GNU Fortran.

Evans' matrix routines, ``3rdparty/polradtran/radmat.f``, are built once as
``polradtran_radmat`` and linked by both ``rt4`` and ``rt3``.  They keep no
state; they are compiled so that their local arrays are on the stack
(``-frecursive`` for GNU, ``-auto`` for Intel), because RT3 and RT4 are
serialised by separate mutexes and may call them concurrently.

RT3 has about 230 MB of static arrays (the scattering-matrix buffer
``SCATBUF`` is 210 MB), which are part of every binary linked with it; the
pages are only touched as far as a problem needs them.

The C++ wrapper ``arts_rt3`` is always built.  When RT3 is disabled,
``rt3::available()`` returns false and ``rt3::get_quadrature()`` and
``rt3::solve()`` throw.  The Python module exists in both cases and raises
``RuntimeError`` the same way.  Python test files whose names contain
``.rt3.`` are collected only with ``ENABLE_RT3=ON``.

Tests
-----

There are two tests, ``cpp.fast.rt3-test``
(``src/core/rt3/test/rt3-test.cpp``) and
``tests/core/rt3/closed-form.rt3.py``.  Every reference is either Evans' own
benchmark output of the original program or a closed form derived in the
test; none is an output of this build.  The C++ test covers:

* **(a) Quadrature exactness** of the three rules (moment errors 1e-14).
* **(b) runmietest**, the Mie case of Evans and Stephens (1991): tau = 1,
  omega = 0.99, an 11-term Mie series, Lambertian albedo 0.1, solar flux
  0.2 pi on the horizontal at mu0 = 0.2 (rt3.f's own conversion of the
  zenith angle 78.46304097 deg, with its truncated pi / 180), gauss with 8
  nodes, aziorder 8, nstokes 4; I, Q, U, V and the fluxes at both levels and
  at 0, 90 and 180 deg (400 entries).
* **(c) runtesta**: a Rayleigh layer over a Mie layer with gas absorption,
  solar and thermal sources at 3 um, gauss with 4 nodes, aziorder 4,
  Lambertian albedo 0.25 at 300 K, three levels (312 entries).
* **(d) Gas-only atmospheres**, thermal source, nstokes 1 to 4, at 50 GHz
  and 30 THz, over a Lambertian surface (double_gauss) and a Fresnel surface
  (gauss with two extra angles), against the exact solution with the Planck
  function linear in optical depth.  Every m > 0 mode and U, V must be
  exactly 0.
* **(e) Single scattering** by a conservative Rayleigh layer of tau = 1e-5
  and 1e-6 over a black surface, lit by the beam, nstokes 1 to 4, at 8
  azimuths, against the exact single-scattering solution.  The reference
  Stokes column is built from the vector geometry of dipole scattering
  (incoherent sum over two incident polarizations, projection on the
  meridional v and h of the outgoing ray), independently of RT3's rotation
  formulas.  It fixes RT3's azimuth origin and sense, Q, and the sign of U.
* **(f) Delta-M with f = 0** is the identity (a regression test of the
  zero-filling in ``GET_SCAT_SET``).
* **(g) The direct beam**, ``F exp(-tau / mu0)``, with the delta-M scaled
  tau ``(1 - omega f) k`` for a Henyey-Greenstein set.
* **(h) 31 error paths**, and the boundary case just inside the FFT limit.

The tables of (b) and (c) were printed by ``OUTPUT_FILE``, which sums the
Fourier series in single precision (REAL*4 ``PHI``, cosine and running sum)
and prints 6 significant digits.  The tolerance of each entry is therefore
half a unit in its last printed digit plus a bound on that REAL*4
summation error,
``sum_m |c_m| (m |phi_f - phi| + 2^-24 m phi + 2^-23)`` plus
``2^-24 sum_m |c_m|`` per term.  The table of (c) was made with the 5-digit
Planck constants.  RT3 evaluates the Planck function only at the interface,
surface and sky temperatures, so the test reproduces the 5-digit Planck
values exactly by passing the temperatures T' at which the exact Planck
function equals the 5-digit one at T.

Measured results:

* (b) runmietest: max deviation 0.94 of the tolerance.  In print half-units
  I agrees to 1.02 and Q to 9.25; the 9.25 is a Q of 1.5e-4 summed from
  terms of 1e-2, i.e. Evans' REAL*4 rounding.  Emulating ``OUTPUT_FILE``'s
  REAL*4 sum reproduces every entry to 1.2 print half-units, and the 32 U, V
  entries at 180 deg, which vanish by symmetry and are pure single-precision
  noise (1e-9 to 1e-11), to 5.6e-5 relative.
* (c) runtesta: max deviation 0.84 of the tolerance; I and Q agree to 1.05
  print half-units, the REAL*4 emulation to 0.99, the noise entries to
  1.1e-5 relative.  With the true temperatures, i.e. the exact Planck
  constants, I differs from the table by up to 1.3e-4 relative.
* (d) Gas-only: 1.2e-14 relative (50 GHz) and 1.1e-13 (30 THz), fluxes the
  same.
* (e) Single scattering: max deviation 1.7 tau relative to max I for both
  tau, i.e. the neglected multiple scattering; the tolerance is 10 tau, and
  the deviation must shrink with tau.
* (f) 6.2e-11 relative, below RT3's round-off floor of eps / max_delta_tau
  = 2.2e-10 (the initial sublayer stores T = 1 - O(max_delta_tau), so
  1 - T has that relative precision).
* (g) Direct beam: 1.3e-16 and 4.2e-16 relative.

The Python test reads the scattering file and the expected output of (b)
directly from ``runmietest`` and compares with the same tolerance (0.94), and
repeats the quadrature, single-scattering (1.7e-6 at tau = 1e-6), gas-only
Fresnel (4.4e-15) and error-path checks through the bindings.

The comparison of VDISORT against RT3, ``cpp.fast.vdisort-rt3-test``, is also
built only with ``ENABLE_RT3=ON``.  See `Mapping to VDISORT inputs`_.  Its
part E solves the problems of Evans' two scripts, read from the scripts, with
VDISORT; see ``src/core/disort-cpp/test/vdisort/COVERAGE.md``.  The
tests of the inputs from ARTS data are listed in `Inputs from ARTS data`_.

Interface
---------

C++ (``#include <rt3.h>``, namespace ``rt3``):

.. code-block:: cpp

  bool available();
  enum class quadrature_type { gauss, double_gauss, lobatto };   // RT3 'G' ('E'), 'D', 'L'
  struct quadrature { Vector mu; Vector weights; };
  quadrature get_quadrature(Index nmu, quadrature_type type);
  Index max_legendre_degree(Index nmu, quadrature_type type);    // RT3's NLEGLIM
  struct scattering_set { Numeric extinction; Numeric scattering; Matrix legendre; };  // [nleg + 1, 6]
  struct lambertian_surface { Numeric albedo; };
  struct fresnel_surface { Complex refractive_index; };
  using surface = std::variant<lambertian_surface, fresnel_surface>;
  struct problem {
    Index nstokes{4}; Index nmu{8}; quadrature_type quad{quadrature_type::gauss};
    Vector extra_mu; Index aziorder{0}; Numeric max_delta_tau{1e-6}; bool delta_m{false};
    Numeric direct_flux{0.0}; Numeric direct_mu{1.0}; bool thermal{true};
    Numeric frequency; Vector height; Vector temperature; Vector gas_extinction;
    std::vector<scattering_set> scattering_sets; ArrayOfIndex layer_scattering_index;
    Numeric sky_temperature; Numeric surface_temperature; surface ground;
  };
  struct result { Vector mu; Vector weights; Tensor4 up; Tensor4 down; Matrix up_flux; Matrix down_flux; };
  result solve(const problem& p);
  Tensor4 azimuth_radiance(const Tensor4& coefficients, const Vector& phi);  // phi in radians

  // #include <rt3_arts.h>: inputs from ARTS data
  scattering_set scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                   const AtmPoint& atm_point, Numeric frequency, Index degree,
                                   Index scattering_angle_count, Numeric normalisation_tolerance);
  struct path_settings {
    Index nstokes{4}; Index nmu{8}; quadrature_type quad{quadrature_type::gauss}; Vector extra_mu;
    Index aziorder{0}; Numeric max_delta_tau{1e-6}; bool delta_m{false}; Index legendre_degree{-1};
    Index scattering_angle_count{512}; Numeric normalisation_tolerance{1e-3};
  };
  problem problem_from_path(const ArrayOfPropagationPathPoint& ray_path,
                            const ArrayOfAtmPoint& atm_path,
                            const ArrayOfPropmatVector& spectral_propmat_path,
                            const AscendingGrid& freq_grid, Index freq_index,
                            const ArrayOfScatteringSpecies& scattering_species,
                            const path_settings& settings, const surface& ground,
                            Numeric surface_temperature, Numeric sky_temperature);

Python (``pyarts3.arts.rt3``) mirrors this: ``available()``,
``QuadratureType``, ``get_quadrature(nmu, type)``,
``max_legendre_degree(nmu, type)``, ``ScatteringSet(extinction, scattering,
legendre)`` and ``ArrayOfScatteringSet``, ``LambertianSurface(albedo)``,
``FresnelSurface(refractive_index)``, ``Problem(...)`` (every field as a
keyword argument with the C++ default), ``RT3Result`` (read-only ``mu``,
``weights``, ``up``, ``down``, ``up_flux``, ``down_flux``), ``solve(problem)``,
``azimuth_radiance(coefficients, phi)``, ``scattering_optics(scattering_species,
atm_point, frequency, degree, scattering_angle_count=512,
normalisation_tolerance=1e-3)``, ``PathSettings(...)`` and
``problem_from_path(...)`` (the arguments of the C++ function).
``Problem.ground`` returns a copy.  ``solve`` releases the GIL.

.. code-block:: python

  import numpy as np
  from pyarts3 import arts

  rt3 = arts.rt3
  rayleigh = np.array([[1.0, -0.5, 0.0, 0.0, 1.0, 0.0],    # F11 F12 F33 F34 F22 F44, l = 0
                       [0.0, 0.0, 1.5, 0.0, 0.0, 1.5],     # l = 1
                       [0.5, 0.5, 0.0, 0.0, 0.5, 0.0]])    # l = 2
  p = rt3.Problem(nstokes=4, nmu=8, aziorder=2,
                  direct_flux=1.0, direct_mu=0.6, thermal=False,
                  frequency=6e14,
                  height=[1.0, 0.0], temperature=[0.0, 0.0], gas_extinction=[0.0],
                  scattering_sets=[rt3.ScatteringSet(0.1, 0.1, rayleigh)],
                  layer_scattering_index=[0],
                  ground=rt3.LambertianSurface(0.0))
  r = rt3.solve(p)
  phi = np.radians([0.0, 45.0, 90.0])
  top_up = np.asarray(rt3.azimuth_radiance(r.up, phi))[0]   # [phi, nmu_total, 4]

Conventions
-----------

**Geometry and azimuth.**

* A right-handed frame with z pointing up, out of the medium.  A ray
  propagating at zenith cosine mu_z (> 0 upward) and azimuth phi has the
  direction ``k = (sin(theta) cos(phi), sin(theta) sin(phi), mu_z)``.
* The direct beam propagates downward toward phi = 0.  So phi is the azimuth
  of the propagation direction of a ray, counted from x toward y, relative to
  the propagation direction of the beam.  For upwelling radiation, phi = 0 is
  forward scattering in azimuth: in (b) the reflected radiance at the
  grazing mu = 0.095 is 0.81 at 0 deg and 0.13 at 180 deg.
* Test (e) pins this, together with the sign of U, against the vector
  construction below.  It is the frame of VDISORT's own polarized tests
  (``_rayleigh_stokes_column`` in ``tests/core/disort/vdisort-interface.py``)
  with the beam azimuth phi0 = 0.

**Stokes basis.**

* [I, Q, U, V], with the meridional plane as reference:
  ``h = k x z / |k x z|`` is horizontal, ``v = h x k`` lies in the plane of
  k and z.
* I = I_v + I_h, Q = I_v - I_h and U = 2 Re(E_v E_h*).  The same basis is
  used in both hemispheres.  A warm dielectric surface emits Q > 0 at
  oblique angles.
* RT3 transports V through F34 (and the Fresnel R4), so V has the sign
  convention of the scattering data.  Evans' runmietest series is in the
  convention of Evans and Stephens (1991), in which RT3 reproduced the V of
  Garcia and Siewert (1989).  For the same physical spheres, ARTS's Mie code
  gives F34 of the opposite sign (see `Inputs from ARTS data`_), so RT3 run
  on ARTS data gives -1 times V in the convention of that paper.  Which of
  the two conventions makes V positive for left-hand circular polarization
  is not tested.

**Streams and hemispheres.**

* ``mu = |cos(zenith)|`` in (0, 1], the same values in both hemispheres.
  The first ``nmu`` are the quadrature nodes, ascending, then the
  zero-weight ``extra_mu`` in the order given.
* ``gauss`` ('G'): the positive half of a 2 nmu-point Gauss-Legendre rule;
  ``double_gauss`` ('D'): an nmu-point Gauss-Legendre rule on [0, 1];
  ``lobatto`` ('L'): the positive half of a 2 nmu-point Lobatto rule,
  including mu = 1.  Weights are for the integral over [0, 1] and sum to 1.
* ``extra_mu`` is allowed only with ``gauss``; RT3 then uses its 'E' type,
  which is Gauss-based.  The extra angles receive scattering and reflection
  but carry no weight.
* ``down`` is radiation propagating downward, toward increasing optical
  depth (RT3's "+", printed with mu > 0 by rt3.f); ``up`` propagates upward
  (RT3's "-", printed with mu < 0).

**Fourier series.**

* ``result.up`` and ``result.down`` are ``[level, m, mu, stokes]``, with
  ``m = 0 .. aziorder``.  The radiance at azimuth phi is
  ``sum_m c_m cos(m phi)`` for I and Q and ``sum_m c_m sin(m phi)`` for U
  and V.  These are RT3's cosine modes of I, Q and sine modes of U, V; the
  m = 0 coefficients of U and V are 0.
* ``azimuth_radiance(coefficients, phi)`` sums the series as ``OUTPUT_FILE``
  does (in double precision), with phi in radians, and returns
  ``[level, phi, mu, stokes]``.
* At phi = 90 deg only the odd modes of U and V contribute, and at 0 and
  180 deg U and V vanish, so Evans' tables (0, 90, 180 deg) do not test the
  even modes of U and V; test (e) does.
* With only thermal sources every m > 0 mode and U, V are exactly 0.

**Layers, levels and units.**

* Layers and levels are top-down.  Level 0 is the top and level ``nlay`` is
  just above the surface.  ``height`` holds ``nlay + 1`` interfaces, of
  which only differences are used; extinctions are in the reciprocal unit.
* ``frequency`` is in Hz.  RT3 is given the wavelength ``1e6 c / f`` um,
  used only by its Planck function.  RT3 works per micrometre; the wrapper
  divides ``direct_flux`` by ``lambda[um] / f`` and multiplies the radiances
  and fluxes by it, so that they are in W m-2 Hz-1 sr-1 and W m-2 Hz-1.  The
  problem is linear in its sources, so a solar-only result does not depend
  on ``frequency``.
* ``up_flux[l, s]`` and ``down_flux[l, s]`` are
  ``2 pi sum_i w_i mu_i c_0`` (only I and Q can be non-zero), and
  ``down_flux`` includes the direct beam ``F exp(-tau / direct_mu)`` in I.

**Scattering sets.**

* ``extinction`` and ``scattering`` are the particle coefficients per unit
  length.  The gas extinction is added per layer; it is scalar and
  unpolarized.  The single-scattering albedo of a layer is
  ``scattering / (extinction + gas)``.
* ``legendre[l, c]`` is ``[nleg + 1, 6]``: the scattering-plane phase matrix

  ::

    [[F11, F12,   0,   0],
     [F12, F22,   0,   0],
     [  0,   0, F33, F34],
     [  0,   0, -F34, F44]]

  with each element a plain Legendre series in cos(Theta),
  ``F_c = sum_l legendre[l, c] P_l(cos(Theta))``, and the columns in the
  order of RT3's scattering files: c = 0 F11, 1 F12, 2 F33, 3 F34, 4 F22,
  5 F44 (``SUM_LEGENDRE`` in ``radscat3.f``).
* The basis is that of the scattering plane with Q = I_par - I_perp, so
  Rayleigh scattering has F12 = -3/4 sin^2(Theta) (the ``rayleigh`` array
  in the example).  The coefficients include the factor 2 l + 1
  (Henyey-Greenstein is ``(2 l + 1) g^l``), and the phase function is
  normalised to 1 over 4 pi: ``legendre[0, 0]`` must be 1.
* RT3 sums F22 and F44 only when they differ from F11 and F33 for some l,
  F34 only for nstokes 4, and only F11 for nstokes 1.  Trailing all-zero
  rows are dropped before the call.
* With ``delta_m``, every set is scaled with M = 2 nmu_total (including the
  extra angles, because RT3 passes its NUMMU): f = legendre[M, 0] / (2 M +
  1), extinction (1 - omega f) k, albedo (1 - f) omega / (1 - omega f),
  diagonal series ``(2 l + 1) (c_l / (2 l + 1) - f) / (1 - f)``,
  off-diagonal (F12, F34) ``c_l / (1 - f)``, degree M - 1.  The beam is
  attenuated with the scaled extinction.  There is no correction of the
  radiances for the truncated peak.

**Sources.**

* The direct beam is on when ``direct_flux > 0``; it requires a Lambertian
  surface.
* ``thermal`` switches the emission of the layers ((1 - omega) B, with B
  linear in optical depth) and of a Lambertian surface ((1 - A) B).  The sky
  (``sky_temperature``) and a Fresnel surface always emit, as in RT3; set
  their temperatures to 0 K to remove them.  The Planck function is 0 for
  T <= 0.  With ``thermal``, every interface temperature must be > 0.

**Surfaces.**

* ``LambertianSurface`` (RT3 'L'): reflection ``2 A mu_j w_j`` of the m = 0
  mode into every stream, I to I only; direct-beam reflection
  ``A F_direct(surface) / pi``.  The diffuse reflection conserves energy on
  the streams only with ``double_gauss`` (``2 sum mu w = 1``).
* ``FresnelSurface`` (RT3 'F'): specular reflection for every mode under a
  medium of index 1, ``[[R1, R2, 0, 0], [R2, R1, 0, 0], [0, 0, R3, -R4],
  [0, 0, R4, R3]]`` with ``R1 = (|r_v|^2 + |r_h|^2) / 2``,
  ``R2 = (|r_v|^2 - |r_h|^2) / 2``, ``R3 = Re(r_v r_h*)``,
  ``R4 = Im(r_v r_h*)``; emission ``[(1 - R1) B, -R2 B, 0, 0]``.  Without a
  beam the field is azimuthally symmetric, so R3 and R4 never act.

**Fortran buffers** (for maintainers): ``SCAT_COEF(6, LDCOEF, set)`` is the
row-major ``legendre`` of each set; ``SCATLAYERS`` is the 1-based set or 0;
``OUTLEVELS`` lists every level; ``UP_RAD``/``DOWN_RAD(s, mu, m, level)`` is
exactly the row-major ``[level, m, mu, s]`` of the result, and
``UP_FLUX``/``DOWN_FLUX(s, level)`` the row-major ``[level, s]``.  For the
'E' type the first ``nmu`` entries of ``MU_VALUES`` are passed as 0.

Inputs from ARTS data
---------------------

``src/core/rt3/rt3_arts.h`` builds RT3 inputs from ARTS scattering species,
atmospheric points and propagation paths.  RT3 takes the scattering matrix of
totally randomly oriented particles as Legendre series, so only species with
TRO data (``ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_gridded``)
can be used.

**scattering_optics** returns the ``scattering_set`` of the species at one
atmospheric point:

* The TRO scattering matrix is evaluated at the nodes of an n-point
  Gauss-Legendre rule in cos(Theta), ``n = scattering_angle_count``, and
  projected, ``c_l = (2 l + 1) / 2 int F(x) P_l(x) dx`` for
  ``l = 0 .. degree``, by that rule.  The projection is exact when every
  element of F is a polynomial of degree ``<= 2 n - 1 - degree`` (Rayleigh:
  ``n >= 3`` for degree 2); otherwise it has the error of the rule, and the
  series is the truncation of F at ``degree``.  Particle habits interpolate F
  linearly from their own (ascending) scattering-angle grid to the nodes.
* ARTS's elements ``[F11, F12, F22, F33, F34, F44]`` are reordered to RT3's
  columns (F11, F12, F33, F34, F22, F44) without any sign change, and the
  series is normalised by ``c_0(F11)`` so that ``legendre[0, 0] = 1``.
* ``extinction`` is K11 and ``scattering`` is K11 - a1, per metre.  The
  phase-function integral ``4 pi c_0(F11)`` must equal K11 - a1 to
  ``normalisation_tolerance`` times K11, otherwise it is an error: RT3
  normalises the series and takes the albedo from the scattering
  coefficient, so a mismatch would change the scattered energy silently.
* Without particles the set has zero extinction and scattering and the
  isotropic, depolarizing series ``[1, 0, 0, 0, 0, 0]``.

**Sign conventions.** The mapping without sign changes is the one under
which RT3 transports ARTS's Stokes vector.  ARTS's laboratory-frame phase
matrix (``to_lab_frame`` in ``scattering/phase_matrix.h``, for the
propagation directions (za, aa)) equals the vector-geometry phase matrix of
the same F in RT3's meridional basis (Q = I_v - I_h, U = 2 Re(E_v E_h*))
with ``mu = cos(za)`` and ``phi = -aa``, element by element, F34 included
(``cpp.fast.vdisort-arts-test``, 6.7e-15 relative over 300 generic
directions; the other azimuth sense misses by 2 and a negated F34 by 0.35
of ``max |Z|``).  ARTS's azimuth runs clockwise seen from above, RT3's
counterclockwise.  RT3 itself uses exactly that vector-geometry matrix
(``cpp.fast.vdisort-rt3-test``).  Against an external reference,
``cpp.fast.rt3-arts-test`` reproduces Evans' runmietest series, the Mie
case of Evans and Stephens (1991) and of Garcia and Siewert (1989): spheres
of refractive index 1.44 at 0.951 um with a gamma distribution of effective
radius 0.2 um and effective variance 0.07, integrated with ARTS's Mie code.
F11, F12 and F33 (and F22 = F11, F44 = F33) agree with Evans' table to
4.8e-9, below half a unit in its 8th decimal, and F34 agrees with the
opposite sign at every degree (4.4e-9; the same sign misses by 0.094).  So
F12, the meridional U and the relative sign of F34 are pinned: RT3 run on
ARTS data gives V of the opposite sign to that convention.  Which convention
makes V positive for left-hand circular polarization is not tested.

**problem_from_path** follows the conventions of ``rt4::problem_from_path``
(see :doc:`dev.rt4`) and of the DISORT workspace methods: one entry per
level, top first, strictly decreasing altitudes (the heights [m]),
temperatures from ``atm_path``, gas extinction from the mean of A of the
unpolarized gas propagation matrix at the two levels, frequency
``freq_grid[freq_index]``.  The scattering set of a layer has the mean
extinction and scattering of its two levels' ``scattering_optics`` and the
scattering-weighted mean of their normalised series, which is the normalised
series of the mean phase matrix; a layer with zero mean extinction is
gas-only.  ``legendre_degree < 0`` selects ``max_legendre_degree(nmu,
quad)``, raised to ``2 nmu_total`` with ``delta_m`` so that RT3 can read the
delta-M fraction.  The problem has thermal emission and no beam; set
``direct_flux``, ``direct_mu`` and ``thermal`` on it for other sources.

Tests:

* ``cpp.fast.rt3-arts-test`` (``src/core/rt3/test/rt3-arts-test.cpp``):
  ARTS's Rayleigh ``GasScatterer`` against RT3's ``rayleigh.sca`` series
  ``[[1, -1/2, 0, 0, 1, 0], [0, 0, 3/2, 0, 0, 3/2], [1/2, 1/2, 0, 0, 1/2,
  0]]`` (1.3e-15); runmietest as above; the path builder against its inputs
  and the scattering-weighted layer mean of Rayleigh and Henyey-Greenstein
  scattering whose mix changes with height (1.1e-16); error paths.
* ``cpp.fast.vdisort-arts-comparison``: RT4, RT3 and VDISORT on one ARTS
  atmosphere (see :doc:`dev.rt4`).  RT3 and VDISORT get the same inputs by
  different routes (Legendre series against vector geometry at the exact
  scattering angles) and agree to their doubling error at
  ``max_delta_tau = 1e-7``: 5.4e-7 (Lambertian) and 5.2e-7 (Fresnel) of
  max I for thermal emission, and 5.8e-7 for a solar beam at mu0 = 0.6 with
  8 Fourier modes and all four Stokes components (tolerance
  ``10 max_delta_tau / mu0``); the difference falls by 11.7 when
  ``max_delta_tau`` falls by 10.
* ``tests/core/disort/arts-native-inputs.rt3.rt4.py`` repeats the closed
  forms through ``pyarts3``.
* ``tests/core/disort/vdisort-polradtran.rt3.rt4.py`` runs VDISORT, RT4 and
  RT3 through their path builders from ``pyarts3`` on ARTS atmospheres
  (Rayleigh, Mie and Henyey-Greenstein species, Lambertian and Fresnel
  surfaces, nstokes 1, 2 and 4, a solar beam).  Run without
  ``ARTS_HEADLESS``, it draws the solutions' plots
  (:func:`pyarts3.plots.cppvdisort.plot`,
  :func:`pyarts3.plots.RT4Result.plot` and
  :func:`pyarts3.plots.RT3Result.plot`) on shared axes.

Limitations
-----------

* **Legendre truncation is an error.** RT3 silently truncates a series (after
  delta-M scaling) to its NLEGLIM, ``max_legendre_degree()``: gauss
  4 nmu - 3, double_gauss 2 nmu - 3, lobatto 4 nmu - 5, at least 1, with nmu
  the quadrature nodes only.  ``solve()`` raises instead when that would drop
  a non-zero coefficient.  Dropping zeros is allowed; since trailing zero
  rows are removed before the call, that only happens with ``delta_m``, and
  RT3 then prints its truncation notice.  Because delta-M gives degree
  2 nmu_total - 1,
  ``delta_m`` with ``double_gauss`` (nmu >= 2), or with ``gauss`` and nmu or
  more extra angles, is rejected unless the scaled series vanishes there.
* **Normalisation.** ``legendre[0, 0]`` must be 1 to 1e-9 after delta-M
  scaling.  RT3's ``CHECK_NORM`` stops the process when the discrete
  normalisation is off by more than 1e-7; for a series within NLEGLIM that
  discrete normalisation equals ``legendre[0, 0] - 1`` up to round-off.
* **Array limits**, all checked before the Fortran call because a Fortran
  ``STOP`` would end the host process (N = nstokes * nmu_total, A =
  aziorder):

  * N <= 64, nlay <= 200, (nlay + 1) N^2 <= 101 * 4096;
  * at most 200 scattering sets, and sets * (A + 1) * 2 N^2 <= 26214400
    (every set is precomputed, also an unused one);
  * with a beam, (A + 1) * 2 N * max(nlay, sets) <= 409600;
  * 2 A + 1 <= 512 with a beam, <= 1024 without (azimuth basis buffers);
  * nleg <= 1023 per set (after dropping trailing zero rows);
  * with A > 0, the degree RT3 sums, ``min(degree, NLEGLIM)``, must be
    <= 251 (its FFT holds 512 azimuth samples).  This only binds for
    nstokes 1 with 64 gauss nodes.
* **Accuracy.** The initial doubling sublayer is first order in its slant
  thickness ``max_delta_tau / mu_min``.  Its transmission is stored as
  1 - O(max_delta_tau), so results carry a round-off floor of about
  eps / max_delta_tau (2e-10 at the default 1e-6).  The thermal source of a
  scattering layer is linear in optical depth.
* **Cost.** For each output level, RT3 adds all layers above and below anew,
  so the cost grows as nlay^2 per azimuth mode; ``solve()`` returns every
  level.
* **Not reentrant.** RT3 uses COMMON blocks, SAVEd FFT tables and large
  static arrays.  Every Fortran call is serialised by one mutex, separate
  from RT4's; concurrent calls are safe but do not run in parallel.
* **Zero pivot.** The zero-pivot ``STOP`` in ``MINVERT`` (``radmat.f``)
  cannot be checked beforehand; it needs an exactly singular 1 - R R.
* **Surfaces.** Only Lambertian with a beam; no BRDF.
* **Azimuth sampling.** RT3 samples the azimuth for its Fourier modes as
  densely as its Legendre degree requires.  That is exact when the
  scattering matrix is regular at forward and backward scattering, as a
  physical matrix is: F12 and F34 vanish there, F22 + F33 has a double zero
  at 180 deg and F22 - F33 one at 0 deg.  For a matrix that is not, the
  laboratory-frame matrix is not a trigonometric polynomial in the azimuth,
  and RT3's modes carry an aliasing error that does not fall with
  ``max_delta_tau``, while VDISORT and RT4 converge with their
  ``azimuth_count``.

RT3 complements RT4: it has the beam, the m > 0 modes and U, V, but only
randomly oriented particles (a scattering-plane phase matrix with six
elements), a scalar unpolarized gas extinction, and no dichroic extinction.

Mapping to VDISORT inputs
-------------------------

``cpp.fast.vdisort-rt3-test``
(``src/core/disort-cpp/test/vdisort/vdisort-rt3-comparison.cpp``) compares
VDISORT with RT3 by this mapping.  Its cases and results are listed in
``src/core/disort-cpp/test/vdisort/COVERAGE.md``.  Feed both solvers the same
scattering sets, so that a comparison tests the solvers and not two
preprocessing chains.  Use ``vdisort::main_data`` (or the low-level
``pyarts3.arts.cppvdisort``) with ``NFourier = aziorder + 1``.  For layer l
with the set (``k = extinction``, ``sigma = scattering``, ``legendre``) and
the gas extinction ``k_g``:

* **Quadrature.**

  * Use RT3 ``double_gauss`` with ``nmu = NQuad / 2``.  The nodes and
    weights are then VDISORT's double-Gauss streams, to round-off.  RT3's
    double-Gauss rule takes Legendre degrees up to ``2 nmu - 3`` only.
  * VDISORT stream ``i`` (mu > 0, upward) is RT3 ``up[.., i]`` and stream
    ``N + i`` is ``down[.., i]``.  RT3's extra angles need its ``gauss``
    rule, which VDISORT does not have, so off-node user angles cannot be
    compared on identical streams.

* **Optical depth.** ``tau_arr`` is the cumulative ``(k + k_g) dz`` at the
  layer bottoms.  RT3 level ``l`` is VDISORT ``tau = tau_l``, with
  ``tau_0 = 0``.
* **Single-scattering albedo.** ``omega = sigma / (k + k_g)``, and 0 for
  gas-only layers.
* **Phase matrix.**

  * Build the ordinary Fourier coefficients on the signed streams (> 0
    upward), without ``epsilon_m``:
    ``C^m(mu_o, mu_i) = (1 / 2 pi) int Z(mu_o, 0; mu_i, phi') cos(m phi') dphi'``
    and ``S^m`` with ``sin(m phi')``.  ``Z = L_out^T F(Theta) L_in`` is the
    lab-frame phase matrix of the series in the meridional basis of
    `Conventions`_, from vector geometry
    (``src/core/disort-cpp/test/vdisort/lab-frame.h``).
  * Then ``phase_matrix = vdisort::combine_phase_matrices(C, S)``.  Do not
    include omega, weights or ``epsilon_m``.
  * For an identical discrete problem, sample at RT3's azimuths
    ``phi' = 2 pi k / NUMPTS``, ``k = 0 .. NUMPTS - 1``, with
    ``NUMPTS = 2 * 2^int(log2(L + 4) + 1)`` for the summed degree L and
    aziorder > 0.  For a series that is regular at Theta = 0 and 180 deg the
    samples do not matter.  Evans' Mie series is regular only to 1e-8, and
    its coefficients change by 2.6e-10 relative.
  * With ``nstokes < 4``, zero the rows and columns of Z from ``nstokes`` on
    before the transform.  RT3 transports the leading block, and VDISORT's
    other components then stay exactly 0.

* **Beam.**

  * ``beam_stokes = [direct_flux / direct_mu, 0, 0, 0]``, the irradiance
    normal to the beam, and ``mu0 = direct_mu``.  ``mu0`` must not be a
    node.
  * ``beam_phase_matrix`` comes from ``C^m(mu_i, -mu0)`` and
    ``S^m(mu_i, -mu0)``.  The beam is a delta in azimuth, so it has cosine
    terms only.  The cosine system (I^c, Q^c, U^s, V^s) gets rows I, Q of
    ``C^m`` and rows U, V of ``S^m``.  The sine system (I^s, Q^s, U^c, V^c)
    gets rows I, Q of ``S^m`` and rows U, V of ``C^m``.  VDISORT applies
    ``epsilon_m`` itself.
  * ``vdisort::combine_beam_phase_matrices`` builds exactly this; the
    diffuse combination of Eq. 81 does not apply to the beam.  See
    ``COVERAGE.md``.

* **Azimuth.** VDISORT's radiance at ``phi0 + psi`` is RT3's at ``psi``.  Both
  are azimuths of the propagation direction; VDISORT's beam propagates
  toward ``phi0`` and RT3's toward 0.
* **Sources.**

  * ``s_poly_coeffs[l] = [c0, c1]`` (Stokes ``[c, 0, 0, 0]``) in the global
    optical depth, ``c1 = (B_bottom - B_top) / dtau_l`` and
    ``c0 = B_top - c1 tau_top``.  VDISORT emits ``(1 - omega) B``, as RT3.
  * ``b_neg[0, 0, i] = [B(T_sky), 0, 0, 0]``, always (RT3 always includes
    the sky).
  * With ``thermal``, ``b_pos[0, 0, i] = [(1 - A) B_s, 0, 0, 0]``.

* **Lambertian surface.** ``vdisort::brdf::lambertian_fourier_modes(A,
  NFourier)``.  Its diffuse operator equals RT3's ``2 A mu w``, and its beam
  reflection RT3's ``A F_direct / pi``.
* **Fluxes.** VDISORT's ``flux()``: ``up`` is ``up_flux[l, 0]``, and
  ``down_diffuse + down_direct`` is ``down_flux[l, 0]``.  The Q fluxes are
  ``2 pi sum_i W_i mu_i u0[i][1]`` over each hemisphere.
* **Delta-M.**

  * RT3's ``delta_m`` corresponds to VDISORT given the scaled set:
    ``k' = (1 - w f) k`` with ``w = sigma / k``, ``sigma' = (1 - f) sigma``,
    and the scaled series (see `Conventions`_).  The direct beam is then
    attenuated with the scaled tau in both.
  * This is not VDISORT's own IMS/TMS-corrected delta-M.
  * With ``double_gauss`` the scaled series must vanish at the degrees
    ``2 nmu - 2`` and ``2 nmu - 1``.  The test uses ``(1 - f)`` Mie plus
    ``f`` times a forward peak truncated at ``M = 2 nmu``, with the dyadic
    ``f = 1/4``.

* **Other surfaces.** RT3 allows only a Lambertian surface with a beam, so
  VDISORT's Fresnel and BRDF beam paths are not compared.
* **V.** VDISORT reproduces RT3's V, which is generated through F34.  That
  shows that both use F34 consistently.  The relative sign of ARTS's F34 is
  in `Inputs from ARTS data`_.
* **Units.** Both are linear in the sources.  The wrapper's results are per
  Hz, so pass ``planck(f, T)`` and the per-Hz ``direct_flux`` to VDISORT
  unchanged.
