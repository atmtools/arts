RT4 reference solver
====================

RT4 is Evans' polarized doubling-adding solver from the polradtran package.
ARTS 3 keeps it only as an external reference for validating other solvers,
in particular the polarized discrete-ordinate solver VDISORT.  It solves the
thermal-only radiative transfer equation in a plane-parallel, azimuthally
symmetric medium, for the Stokes components [I] or [I, Q].  It has a core C++
interface, the namespace ``polradtran::rt4`` (``src/core/polradtran/rt4``;
``rt4::`` below), and a low-level Python interface (``pyarts3.arts.rt4``).  It
has no workspace methods, variables or agendas.  Its equations are in
:doc:`concept.rt4`.

Provenance
----------

* **Original code.** K. F. Evans, polradtran (RT3/RT4), 1996, distributed
  from https://nit.coloradolinux.com/polrad.html under the MIT licence
  (``3rdparty/polradtran/LICENSE``).  RT4 is briefly described by
  :cite:t:`Evans1995`.  See Evans' ``3rdparty/polradtran/README``.
  Of ``PolRadTran.tar`` (sha256 ``7b0eff79...a6f7cff9d``),
  ``3rdparty/polradtran`` keeps the licence, ``README``, the four test
  scripts and ``cl340d14.dda``, unchanged.  The Fortran was removed once the
  port to C++ was complete; the history of ``README.ARTS`` there lists the
  ARTS changes to it.
* **ARTS 2.6 changes** to the Fortran, by J. Mendrok and others:

  * optical properties are passed in memory instead of read from files;
  * new surface types ``'S'`` (specular with a fixed reflection matrix) and
    ``'A'`` (an externally supplied reflection operator and emission);
  * output at all levels;
  * a choice of optics set per layer;
  * larger array limits;
  * a hard stop for ``NSTOKES > 2``;
  * extra zero-weight angles at the end of the stream list.
* **ARTS 3 changes** to the Fortran, each marked ``c ARTS3:`` in the source,
  before it was ported:

  * Planck constants computed from the exact SI h, c and k.  The original
    5-digit constants give a bias of about 3e-5, roughly 8 mK at 250 K.
  * ``NUUMMU`` declared ``INTEGER``; it was implicitly typed before.
  * ``UP_RAD``/``DOWN_RAD`` declared with the extent that is actually
    written.
  * A new ``rt4_c_interface.f90`` with ``ISO_C_BINDING`` entry points, so
    the build needs neither ``-fdefault-integer-8`` nor the hidden
    ``CHARACTER`` length ABI.
  * A rewritten ``CMakeLists.txt``.  ``radmat.f`` is built once as
    ``polradtran_radmat`` for both RT4 and RT3 (:doc:`dev.rt3`).

  ``SYMMETRIC`` is still hard-coded to ``.TRUE.``.

* **Port to C++.** RT4 is ported to C++, one routine at a time, on matpack
  and rtepack types, with the method unchanged; ``rt4::solve`` calls no
  Fortran.  Over the whole port the radiances of 91 test problems changed
  by at most 2.1e-12 of I (median 1.2e-14), the most in an optically thick,
  strongly scattering layer through LAPACK's inverse.  The steps:

  * ``RADTRANO`` is ``rt4::radtrano`` (``src/core/polradtran/rt4/radtran4.h``).  It
    follows the Fortran step by step and calls the same subroutines.
    All of them are ported, so it calls no Fortran and keeps no state.
    It takes no counts: ``NSTOKES``, ``NUMMU``, ``NUUMMU``, ``NUM_LAYERS``
    and ``NSL`` are the extents of its arrays, with the extra angles as an
    input of their own.  ``MZERO`` is ``= 0.0``, ``MCOPY`` is ``=`` and
    ``MIDENTITY`` is ``identity`` (which sets a square matrix to a
    multiple of the identity and returns it, so that ``1 - R R`` is
    ``mult(identity(y), R, R, -1.0, 1.0)``).  Its work arrays are sized to
    the problem; its ``STOP`` checks throw, except those of the static-array
    sizes (``MAXV``, ``MAXLAY``, ``MAXLM``), which the port does not have.
  * The quadratures (``DOUBLE_GAUSS_QUADRATURE``,
    ``GAUSS_LEGENDRE_QUADRATURE``, ``LOBATTO_QUADRATURE``) are ARTS's: the
    positive half of ``scattering::DoubleGaussQuadrature``,
    ``GaussLegendreQuadrature`` or ``LobattoQuadrature`` of degree
    ``2 nmu`` (``polradtran::get_quadrature``, which ``rt4::radtrano`` and
    ``rt4::solve`` share).  They are RT4's rules and differ from Evans'
    routines by rounding: for ``nmu`` up to 64 the nodes by at most
    4.4e-16 and the weights by at most 2.4e-12 relative.  Against a
    40-digit reference, Evans' Gauss weights are off by up to 2e4 ulps (his
    Newton iteration takes P' at the last iterate but one), ARTS's by at
    most 4; ARTS's smallest ``gauss`` node, cos(theta) near theta = pi/2,
    is off by up to 4e-16 absolute, Evans' by 3 ulps.
  * ``rt4::radtrano`` works in SI: it takes the frequency [Hz] and its
    radiances are W m-2 Hz-1 sr-1.  The layers' Planck function
    (``PLANCK_FUNCTION`` in ``RADTRANO``) is ARTS's ``planck()``, as are
    those of the ground and the sky as their routines were ported.
    Both functions use the exact SI h, c and k, but ``PLANCK_FUNCTION``
    evaluates ``exp(x) - 1``: against a 40-digit reference it is off by
    4.5e-13 at 1 GHz and 250 K, ``planck()`` (``expm1``) by 2e-16.
    ``RADTRANO`` gave 0 below 0 K, ``planck()`` a negative value, so
    ``rt4::radtrano`` rejects negative temperatures.
  * The ground is external to ``rt4::radtrano``: it takes ``SURF_REFLECT``
    and ``GND_RADIANCE`` (RT4's ground type ``'A'``) for every kind of
    ground, and ``GROUND_TEMP``, ``GROUND_TYPE``, ``GROUND_ALBEDO``,
    ``GROUND_INDEX`` and ``GROUND_REFLEC`` are gone.  Every ground routine
    of RT4 makes the same surface layer (no reflection from above, the
    identity as transmission, no source); only the reflection back up,
    ``REFLECT(..., 2)``, and the ground's radiance depend on the ground, and
    ``EXTERNAL_SURFACE`` makes that layer from exactly these two.
    ``rt4::ground_surface`` (``src/core/polradtran/rt4/radutil4.h``) makes them from
    an ``rt4::surface``: for the Lambertian, Fresnel and specular grounds
    with the ``*_SURFACE`` and ``*_RADIANCE`` routines, for a
    ``discrete_surface`` as given.  This is bit-identical to the ground
    types inside ``RADTRANO``.  The ``*_SURFACE`` routines are ported as
    ``rt4::*_surface_layer`` (they make the ground as a layer for the
    adding; ``polradtran::fresnel_surface`` and ``rt4::specular_surface`` are the
    types of ``rt4.h``); those RT3 shares are ``polradtran``'s
    (``radutil.h``, below).  ``LAMBERT_SURFACE`` and ``LAMBERT_RADIANCE``
    are ``polradtran::lambert_surface_layer`` and
    ``polradtran::lambert_radiance`` (``radutil.h``; RT4's
    ``LAMBERT_RADIANCE`` is RT3's in mode 0 with the thermal source alone),
    as is, the radiance with
    ``planck()`` in SI instead of ``PLANCK_FUNCTION`` per micrometre (and
    rejecting a negative ground temperature, where ``PLANCK_FUNCTION``
    gave 0); the reflection is bit-identical, the radiances of a Lambertian
    ground change by at most 3.3e-15 of I.  ``FRESNEL_SURFACE`` and
    ``FRESNEL_RADIANCE`` are ``polradtran::fresnel_surface_layer`` and
    ``polradtran::fresnel_radiance`` (RT3's in mode 0), with ARTS's ``fresnel()`` amplitudes
    (``physics_funcs.h``) and ``rtepack::fresnel_reflectance``, whose
    Mueller matrix is RT4's (``R1`` and ``R2`` in the [I, Q] block, ``R3``
    and ``R4`` in the [U, V] block), and the emission ``(1 - R) B``.
    ``fresnel()`` was made exact for this: it used the real Snell angle of
    ``Re n2``, which for an absorbing ground is off by up to 5e-3 in
    reflectivity (about 1.4 K over water); it now uses the complex
    transmitted cosine, like RT4 (for a real ``n2``, as
    ``spectral_surf_reflFlatRealFresnel`` passes, the change is rounding).
    The reflection matches RT4's to 3.4e-15, the radiances of a Fresnel
    ground change by at most 3e-15 of I.  ``SPECULAR_SURFACE`` and
    ``SPECULAR_RADIANCE`` are ``rt4::specular_surface_layer``, as is
    (bit-identical), and ``rt4::specular_radiance``, the emission
    ``(1 - R) B`` with rtepack (RT4's ``[(1 - R(I, I)) B, -R(Q, I) B]``;
    for more than two Stokes components it also gives the U and V that RT4
    set to 0); the radiances of a specular ground change by at most 3.1e-15
    of I.  ``rt4::ground_surface`` calls no Fortran.  ``EXTERNAL_SURFACE``,
    which makes the surface layer in ``RADTRANO`` from ``SURF_REFLECT``, is
    ``polradtran::external_surface_layer``, as is (bit-identical), without the
    ``RADIANCE`` argument that ``EXTERNAL_SURFACE`` does not use (the
    ground's radiance goes to ``INTERNAL_RADIANCE``).  ``THERMAL_RADIANCE``,
    the sky, is ``polradtran::thermal_radiance`` (RT3's in mode 0) with ``planck()`` in SI, so
    ``rt4::radtrano`` has no unit conversion left; its radiances change by
    at most 6e-16 of I, and an isothermal atmosphere now reproduces
    ``planck()`` to 1.7e-16 (1.6e-14 with ``PLANCK_FUNCTION`` for the sky
    and the ground).
  * ``NONSCATTER_LAYER``, ``INITIAL_SOURCE`` and ``INITIALIZE`` are
    ``polradtran::nonscatter_layer`` (``radintg.h``, RT3's in mode 0),
    ``rt4::initial_source`` and ``rt4::initialize``
    (``src/core/polradtran/rt4/radintg4.h``), as is, with ``= 0.0``
    for ``MZERO`` and ``Constant::two_pi`` for ``C``; they call no Fortran.
    All are bit-identical to the Fortran, except where gfortran on glibc
    vectorises ``NONSCATTER_LAYER``'s ``DEXP`` to libmvec's (up to 3.5 ulp
    off): its source, whose terms cancel to about the Planck function
    times the path, then differs by up to 5e-9 at a path of 1e-4.
  * ``DOUBLING_INTEGRATION`` is ``polradtran::doubling_integration``
    (``radintg.h``), RT3's with the linear (thermal) source alone, which
    has an overload with RT4's arguments; with matpack for Evans' matrix helpers: ``MCOPY`` is
    ``=``, ``MSCALARMULT`` and ``MADD`` on vectors ``*=`` and ``+=``,
    ``MINVERT`` (LINPACK ``DGEFA``/``DGEDI``) is ``inv_inplace`` (LAPACK
    ``dgetrf``/``dgetri``), and ``MMULT`` is ``mult`` (``DGEMM``), whose
    ``alpha`` and ``beta`` absorb the ``MIDENTITY`` and ``MSUB`` of
    ``1 - R R`` and the ``MADD`` that follows a product.  The row-major
    matpack matrix of a Fortran matrix is its transpose, so ``MMULT``'s
    ``C = A B`` is ``mult(C, B, A)``, the same ``DGEMM`` call, and a
    matrix-vector product ``y = A x`` is ``mult(y, transpose(A), x)``,
    ``DGEMV`` (``MMULT`` used ``DGEMM`` with one column; OpenBLAS gives the
    same result except for 1 x 1, by 1 ulp).  This mapping goes when
    ``COMBINE_LAYERS`` and ``INTERNAL_RADIANCE``, which read the same
    arrays, are ported and the matrices can be stored as the equations read.
    Two parts change the numbers, measured
    separately on the full ``rt4.solve``: with the products and additions
    unfused and LINPACK's inverse, the port is bit-identical; the fused
    ``beta`` (one rounding fewer per product) changes the radiances by at
    most 2.2e-15 of I; LAPACK's inverse instead of LINPACK's by at most
    2e-12 of I (median 4e-18), in an optically thick, strongly scattering
    layer, where about 26 doublings each invert a poorly conditioned
    ``1 - R R``.  Each doubling about squares T, doubling its relative
    error, so n doublings amplify rounding by 2^n: against the Fortran
    routine directly the difference grows from 2e-14 for 6 doublings to
    1.2e-9 for 24 (1.3e-11 on Apple arm64 with OpenBLAS), and in quad
    precision both are 1.9e-9 off.
  * ``COMBINE_LAYERS`` is ``polradtran::combine_layers`` (``radintg.h``), in the
    same way as ``polradtran::doubling_integration``.  Against the Fortran
    routine it differs by at most 7.7e-16 (one combination does not
    amplify the rounding as repeated doublings do); the radiances of
    ``rt4.solve`` change by at most 3e-15 of I.
  * ``INTERNAL_RADIANCE`` is ``polradtran::internal_radiance`` (``radintg.h``),
    in the same way: the matrix-vector products are ``DGEMV`` with
    ``beta`` absorbing the ``MADD`` after them.  Against the Fortran
    routine it differs by at most 8.4e-16; the radiances of ``rt4.solve``
    change by at most 1.3e-15 of I.  With it, the Fortran mutex of
    ``rt4::solve`` is gone.
  * All work arrays of ``rt4::radtrano`` and the routines it calls are in
    one ``polradtran::workdata`` (``src/core/polradtran/polradtran_workdata.h``,
    shared with RT3), grouped by
    lifecycle: the layers' R, T and S, ``RADTRANO``'s arrays on the
    streams, ``DOUBLING_INTEGRATION``'s linear source, and the scratch that
    ``DOUBLING_INTEGRATION``, ``COMBINE_LAYERS`` and ``INTERNAL_RADIANCE``
    use one at a time (``X``, ``Y``, ``GAMMA``, two vectors and LAPACK's
    workspace, as RT4's COMMON blocks shared them).  ``rt4::radtrano``
    sizes it, allocating only where an array grows, so one work data kept
    over repeated calls with the same streams (as over frequency) allocates
    nothing; ``rt4::solve`` makes one per call.  Reusing one over the 84
    cases of ``cpp.fast.rt4-radtrano-test``, of different sizes, gives the
    same bits as a fresh one per case.
  * What RT4 and RT3 share is one C++, in the namespace ``polradtran``
    (``src/core/polradtran``, the library ``arts_polradtran``): the
    routines that are the same Fortran in both, line by line
    (``COMBINE_LAYERS``, ``INTERNAL_RADIANCE``, ``LAMBERT_SURFACE`` and
    ``FRESNEL_SURFACE``; RT3's names have the prefix ``RT3_``), those of
    which RT4's is RT3's in the azimuth mode 0 with the thermal source
    alone (``DOUBLING_INTEGRATION``, ``NONSCATTER_LAYER``,
    ``THERMAL_RADIANCE``, ``LAMBERT_RADIANCE`` and ``FRESNEL_RADIANCE``;
    RT4 passes mode 0), ``EXTERNAL_SURFACE``, which ``rt3::radtran`` also
    uses, and the work data, of which ``rt3::rt3_workdata`` extends RT4's.  Both
    ``cpp.fast.rt4-radtrano-test`` and ``cpp.fast.rt3-radtran-test``
    compare it with their own Fortran, and sharing it left every result
    of both ports bit-identical.
    The level loop of ``RADTRANO`` and ``RADTRAN``, which adds the layers
    above and below a level and calls ``INTERNAL_RADIANCE``, is
    ``polradtran::level_radiance``, and their initial sublayer of a
    scattering layer and its number of doublings
    ``polradtran::initial_sublayer``.  The interfaces share their streams
    and grounds (``polradtran.h``: ``quadrature_type``, ``quadrature``,
    ``get_quadrature``, ``lambertian_surface`` and ``fresnel_surface``, in
    Python ``pyarts3.arts.polradtran``) and the layers of an ARTS
    propagation path (``polradtran::layers_from_path`` in
    ``polradtran_arts.h``, behind both ``problem_from_path``).

  The port was tested against the Fortran ``RADTRANO``, built beside it, on
  random inputs over every branch of ``RADTRANO``.  The port was
  bit-identical until the quadratures were replaced.  Mathematically
  equivalent evaluations (FMA contraction, vectorised libm functions, other
  BLAS kernels) are accepted: the tolerances allow 16 epsilon times what a
  computation amplifies rounding by.  Every output had to agree to 1e-11 of
  the largest value in it, for the replaced quadratures and ``planck()``,
  plus 16 epsilon times 2^n for the layer doubled n times most (the largest
  difference was 4e-13 on Apple arm64 with OpenBLAS, 1.5e-9 on AMD x86_64
  with MKL), and each ported routine was checked against Evans' routine.
  A step meant to leave the numbers alone was also checked bit for bit
  against the step before; one that changes them, like the quadratures, had
  the change measured.  When the Fortran was removed, its outputs for nine
  cases that cover every branch (the four grounds, the three quadratures,
  extra angles, gas, thin, thick and shared layers, a coarse
  ``max_delta_tau`` and a 0 K top) were kept as constants
  (``src/core/polradtran/rt4/test/rt4-radtrano-reference.h``), with which
  ``cpp.fast.rt4-radtrano-test`` compares the port to the same tolerances.

  Evans' original programs ``rt4.f`` and ``scatcnv.f``, built from the tar,
  ran his two RT4 scripts and reproduced his tables exactly, before they
  were removed with the Fortran.  ``cpp.fast.polradtran-rt4-arts`` gives
  ARTS's RT4 the optics that ``rt4.f`` reads for ``runtestc`` (cirrus of
  horizontally oriented ice columns at 340 GHz from the DDA file
  ``cl340d14.dda``, 8 Lobatto streams, a tropical atmosphere over land) and
  reproduces his table to 0.005 K.  ``runtestr`` (a 2 mm/h rain layer of
  spherical drops at 85 GHz over water, 8 Gauss streams and a Fresnel
  surface), whose RT4 scattering file Evans' ``scatcnv`` made, is solved
  from its Mie Legendre series by
  ``tests/core/disort/evans-benchmarks.rt3.rt4.py``.

Build
-----

RT4 is C++ and always built.  Its MIT licence is registered as bundled
code, see :doc:`dev.licenses`; it is allowed in LGPL builds, as polradtran
is MIT licensed (``3rdparty/polradtran/LICENSE``) and the ARTS changes fall
under the ARTS licence.

There are two tests of the solver wrapper:

* ``cpp.fast.rt4-test`` (``src/core/polradtran/rt4/test/rt4-test.cpp``);
* ``tests/core/rt4/closed-form.rt4.py``.

Every reference in both is a closed form derived in the test, never an RT4
output.  The C++ test covers:

* quadrature exactness;
* Planck against an isothermal blackbody;
* gas-only layers;
* doubling convergence;
* a layout test with lower-triangular, stream-dependent K, Stokes-asymmetric
  absorption and a forward-only phase matrix;
* Fresnel, specular, Lambertian and discrete surfaces;
* isothermal Kirchhoff with Rayleigh and a non-reciprocal Rayleigh variant,
  also with ``nstokes * nmu_total = 82`` and with 450 layers, beyond the
  Fortran's array sizes;
* 14 error paths.

The Python test repeats the quadrature, Fresnel, layout, Kirchhoff and
error-path checks through the bindings.

The comparison of VDISORT against RT4 is ``cpp.fast.vdisort-rt4-test``.  See
`Mapping to VDISORT inputs`_.  The
tests of the inputs from ARTS data are listed in `Inputs from ARTS data`_.

Interface
---------

C++ (``#include <rt4.h>``, namespace ``polradtran::rt4``, with the streams
and the grounds it shares with RT3 in ``polradtran.h``, namespace
``polradtran``):

.. code-block:: cpp

  // polradtran.h, namespace polradtran (shared with RT3)
  enum class quadrature_type { gauss, double_gauss, lobatto };   // 'G', 'D', 'L'
  struct quadrature { Vector mu; Vector weights; };
  quadrature get_quadrature(Index nmu, quadrature_type type);
  struct lambertian_surface { Numeric albedo; };
  struct fresnel_surface { Complex refractive_index; };

  // rt4.h, namespace polradtran::rt4
  inline constexpr Index down = 0, up = 1;
  struct layer_optics { Tensor4 extinction; Tensor3 absorption; Tensor6 phase; };
  struct specular_surface { Matrix reflectivity; };
  struct discrete_surface { Tensor4 reflection; Matrix emission; };
  using surface = std::variant<lambertian_surface, fresnel_surface, specular_surface, discrete_surface>;
  struct problem {
    Index nstokes{2}; Index nmu{8}; quadrature_type quad{quadrature_type::double_gauss};
    Vector extra_mu; Numeric max_delta_tau{1e-6}; Numeric frequency;
    Vector height; Vector temperature; Vector gas_extinction;
    std::vector<layer_optics> optics; ArrayOfIndex layer_optics_index;
    Numeric sky_temperature; Numeric surface_temperature; surface ground;
  };
  struct result { Vector mu; Vector weights; Tensor3 up; Tensor3 down; };
  result solve(const problem& p);

  // #include <rt4_arts.h>: inputs from ARTS data
  layer_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                 const AtmPoint& atm_point, Numeric frequency,
                                 const Vector& mu, Index nstokes, Index azimuth_count);
  struct path_settings {
    Index nstokes{2}; Index nmu{8}; quadrature_type quad{quadrature_type::double_gauss};
    Vector extra_mu; Numeric max_delta_tau{1e-6}; Index azimuth_count{64};
  };
  problem problem_from_path(const ArrayOfPropagationPathPoint& ray_path,
                            const ArrayOfAtmPoint& atm_path,
                            const ArrayOfPropmatVector& spectral_propmat_path,
                            const AscendingGrid& freq_grid, Index freq_index,
                            const ArrayOfScatteringSpecies& scattering_species,
                            const path_settings& settings, const surface& ground,
                            Numeric surface_temperature, Numeric sky_temperature);

Python (``pyarts3.arts.rt4``, with what RT3 shares in
``pyarts3.arts.polradtran``) mirrors this:

* ``polradtran.QuadratureType`` (``gauss``, ``double_gauss``, ``lobatto``);
* ``polradtran.get_quadrature(nmu, type)``, which returns a
  ``polradtran.Quadrature`` with ``mu`` and ``weights``;
* the surfaces ``polradtran.LambertianSurface(albedo)``,
  ``polradtran.FresnelSurface(refractive_index)``, and RT4's own
  ``SpecularSurface(reflectivity)`` and ``DiscreteSurface(reflection,
  emission)``;
* ``down``/``up``;
* ``LayerOptics(extinction, absorption, phase)`` and ``ArrayOfLayerOptics``;
* ``Problem(...)``, which takes every field as a keyword argument with the
  C++ default;
* ``RT4Result``, with read-only ``mu``, ``weights``, ``up`` and ``down``;
* ``solve(problem)``;
* ``scattering_optics(scattering_species, atm_point, frequency, mu,
  nstokes=2, azimuth_count=64)``, ``PathSettings(...)`` (every field as a
  keyword argument with the C++ default) and ``problem_from_path(ray_path,
  atm_path, spectral_propmat_path, freq_grid, freq_index,
  scattering_species, settings, ground, surface_temperature,
  sky_temperature)``.

Every array attribute accepts numpy arrays and lists.  ``Problem.ground``
returns a copy, so assign a new surface to change it.  ``solve`` releases the
GIL.

.. code-block:: python

  import numpy as np
  from pyarts3 import arts

  rt4 = arts.rt4
  p = rt4.Problem(nstokes=2, nmu=8, extra_mu=[1.0], frequency=89e9,
                  height=np.array([3000.0, 2000.0, 1000.0, 0.0]),       # m, top-down
                  temperature=np.array([220.0, 240.0, 265.0, 285.0]),   # K at interfaces
                  gas_extinction=np.array([1e-4, 3e-4, 5e-4]),          # 1/m per layer
                  layer_optics_index=[-1, 0, -1],
                  sky_temperature=2.725, surface_temperature=290.0,
                  ground=arts.polradtran.FresnelSurface(3.0 + 0.2j))

  n = p.nmu + len(p.extra_mu)
  q = arts.polradtran.get_quadrature(p.nmu, p.quad)
  w = np.append(q.weights, 0.0)
  sigma, kabs = 6e-4, 4e-4                       # isotropic scattering, per metre
  phase = np.zeros((2, 2, n, n, 2, 2))
  phase[..., 0, 0] = sigma / (4 * np.pi)
  ext = np.zeros((2, n, 2, 2))
  ext[..., 0, 0] = ext[..., 1, 1] = sigma + kabs
  absorption = np.zeros((2, n, 2))               # energy conservation on the streams
  absorption[..., 0] = sigma + kabs - 2 * np.pi * np.einsum("i,ohij->hj", w, phase[..., 0, 0])
  p.optics = [rt4.LayerOptics(ext, absorption, phase)]

  r = rt4.solve(p)
  top_up = np.asarray(r.up)[0]                   # [nmu_total, 2], W m-2 Hz-1 sr-1

Conventions
-----------

**Stokes basis.**

* [I, Q], with the meridional plane as reference: "vertical" polarization
  lies in the plane of the ray and the z axis.
* I = I_v + I_h and Q = I_v - I_h.
* The same basis is used in both hemispheres, so Q does not change sign
  between up and down.  A warm dielectric surface emits Q > 0 at oblique
  angles.
* U and V are not computed.

**Hemispheres.**

* Index ``down`` (0) is radiation propagating downward, toward increasing
  optical depth.  This is RT4's "+".
* ``up`` (1) is propagating upward, RT4's "-".
* ``result.up[l, i]`` is what a sensor at level ``l`` sees looking down at
  nadir angle ``acos(mu_i)``.

**Streams.**

* Streams are given per hemisphere as ``mu = |cos(zenith)|`` in (0, 1].
  Both hemispheres use the same ``mu`` values.
* The first ``nmu`` are the quadrature nodes (``get_quadrature``), in
  ascending order.  None of the rules includes ``mu = 0``.
* The zero-weight ``extra_mu`` angles follow, in the order given; each must
  be in (0, 1].  ``nmu_total = nmu + len(extra_mu)``.
* Weights are for the integral over mu in [0, 1] and sum to 1.  The 2 pi
  azimuth factor is not included.
* The three rules:

  * ``double_gauss`` ('D'): an nmu-point Gauss-Legendre rule mapped to
    [0, 1];
  * ``gauss`` ('G'): the positive half of a 2 nmu-point Gauss-Legendre rule
    on [-1, 1];
  * ``lobatto`` ('L'): the positive half of a 2 nmu-point Lobatto rule; it
    includes mu = 1.
* The extra angles receive scattering and reflection but contribute nothing
  to the angular integrals.  They are output directions only, except in
  ``DiscreteSurface``, whose columns are applied without weights (see
  Surfaces).

**Layers and levels.**

* Both are ordered top-down.  Level 0 is the top of the atmosphere and level
  ``nlay`` is just above the surface.
* ``height`` holds ``nlay + 1`` interfaces.  Only ``|height[l] - height[l+1]|``
  is used.  Any length unit works if extinctions are given in its
  reciprocal.
* ``temperature`` holds ``nlay + 1`` interface temperatures, which must be
  > 0.
* ``gas_extinction`` holds ``nlay`` values, which must be >= 0.
* ``layer_optics_index`` holds ``nlay`` values: an index into ``optics``, or
  a negative value for a gas-only layer.

**Optics shapes** (``layer_optics``), all per unit length:

* ``extinction[h, mu, row, col]``: the extinction matrix K for propagation in
  hemisphere h at mu.
* ``absorption[h, mu, s]``: the absorption vector a, which RT4 multiplies by
  the Planck function.
* ``phase[h_out, h_in, mu_out, mu_in, s_out, s_in]``.

The gas extinction is added by RT4 to every Stokes diagonal of K and to the
I component of a.  It is scalar and unpolarized.

**Phase-matrix quadrants and normalisation.**

* The four ``(h_out, h_in)`` blocks are RT4's quadrants:

  * ``(down, down)`` is q=1 (+ <- +);
  * ``(down, up)`` is q=2 (+ <- -);
  * ``(up, down)`` is q=3 (- <- +);
  * ``(up, up)`` is q=4 (- <- -).

  The Fortran buffer index is q = 2 h_out + h_in + 1.
* ``phase`` is the azimuthal mean ``(1 / 2 pi) int Z dDelta-phi`` of the
  phase matrix in the meridional basis, per unit length and per steradian.
  It includes the number density, is not normalised to 4 pi, and contains no
  quadrature weights.
* The scattering source into ``(h, mu_i)`` is
  ``2 pi sum_j sum_h' w_j Z(h <- h')(i, j) I(h', mu_j)``.
* Energy conservation on the streams reads
  ``K11(h, mu_j) = a1(h, mu_j) + 2 pi sum_i w_i [Z(up <- h) + Z(down <- h)](1, i; 1, j)``.
  Neither RT4 nor the wrapper checks or enforces it; ARTS 2.6 rescaled the
  phase matrices to satisfy it.
* For an isothermal Kirchhoff test (I = B, Q = 0 on the streams), use
  ``a_s(h, mu_j) = K_sI - 2 pi sum_{i, h'} w_i Z_sI(h <- h')(j, i)``.

**Planck source.**

* Within a layer the Planck function B, not the temperature, is linear in
  optical depth: B is linear in height, and the layer is homogeneous.
* Gas-only layers are solved analytically, with ``exp(-tau / mu)`` and the
  exact linear-in-tau source integral.
* Scattering layers are doubled from a first-order initial sublayer, so the
  source is piecewise constant per sublayer.  The error is first order in
  that sublayer's slant thickness, which is at most
  ``max_delta_tau / mu_min``.
* RT4 chooses the number of doublings from ``extinction[down, 0, 0, 0]``
  plus the gas extinction only.

**Units.**

* ``frequency`` is in Hz, and the radiances (``result.up``,
  ``result.down`` and ``DiscreteSurface.emission``) are W m-2 Hz-1 sr-1;
  ``rt4::radtrano`` works in these units (see the port above).
* The sky is an isotropic, unpolarized blackbody at ``sky_temperature``.
  RT4's Planck function for the sky and the ground is 0 for a temperature
  <= 0.
* ``result.weights`` comes from a second quadrature call, because RADTRANO
  does not return its weights.

**Surfaces.**

* ``polradtran.LambertianSurface`` (RT4 'L'): reflection ``2 A mu_j w_j`` into every
  stream, I to I only (depolarizing); emission ``[(1 - A) B_s, 0]``.
* ``polradtran.FresnelSurface`` (RT4 'F'): the medium above has index 1, and the
  reflection is specular and stream by stream.
  ``R = [[R1, R2], [R2, R1]]``, with ``R1 = (|r_v|^2 + |r_h|^2) / 2`` and
  ``R2 = (|r_v|^2 - |r_h|^2) / 2``.  Emission is ``[(1 - R1) B_s, -R2 B_s]``.
  The sign of ``Im n`` does not matter for [I, Q].
* ``SpecularSurface`` (RT4 'S'): a fixed ``R(out, in)`` of shape
  ``[nstokes, nstokes]``, applied specularly to every stream.  Emission is
  ``[(1 - R(I, I)) B_s, -R(Q, I) B_s]``.
* ``DiscreteSurface`` (RT4 'A'): ``I_up(i) = sum_j reflection[i, j] I_down(j)
  + emission[i]``.

  * ``reflection`` is ``[nmu_total out, nmu_total in, s_out, s_in]`` and
    must include all quadrature factors; a Lambertian surface is
    ``2 A mu_j w_j``.
  * The columns of extra angles are applied without weights.  Set them to
    zero to keep the extra angles as pure outputs.
  * ``emission`` is ``[nmu_total, nstokes]`` in W m-2 Hz-1 sr-1.
  * ``surface_temperature`` is not used.

**Fortran buffers** (column-major, first index fastest; for maintainers):

* ``EXTINCT_MATRIX(row, col, mu, hem, set)`` and
  ``EMIS_VECTOR(s, mu, hem, set)``;
* ``SCATTER_MATRIX(s_out, mu_out, s_in, mu_in, q, set)``;
* ``SCATLAYERS`` is the 1-based set index, or 0 for gas-only;
* ``SURF_REFLECT(s_out, mu_out, s_in, mu_in)``;
* the specular R is passed row-major as is, because RT4 reads it
  transposed;
* ``UP_RAD``/``DOWN_RAD(s, mu, level)`` is exactly the row-major
  ``[level, mu, s]`` of the result.

Inputs from ARTS data
---------------------

``src/core/polradtran/rt4/rt4_arts.h`` builds RT4 inputs from ARTS scattering species,
atmospheric points and propagation paths.  It adds no physics: the optics
are ARTS's bulk scattering properties in the laboratory frame
(``ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_gridded``),
so azimuthally randomly oriented (ARO) species work as well as totally
randomly oriented (TRO) ones.  RT4 is the only one of RT4, RT3 and VDISORT
that can take ARO particles.

**scattering_optics** returns the ``layer_optics`` of the species at one
atmospheric point on the streams ``mu`` (the quadrature nodes followed by the
extra angles), per metre and steradian:

* RT4's stream ``(down, mu)`` propagates toward the surface, which is ARTS's
  propagation zenith angle ``180 - acos(mu)`` deg; ``(up, mu)`` is
  ``acos(mu)``.  ARTS's scattering data take propagation directions, on
  zenith-angle grids that are strictly ascending (``ZenGrid``), so the
  values of ``mu`` must be distinct.
* ``extinction[h, i]`` is the ARO extinction matrix ``[[K11, K12], [K12,
  K11]]`` for propagation in ``(h, mu_i)``, ``absorption[h, i]`` is
  ``[a1, a2]`` (TRO species give K11 and a1 only).
* ``phase[ho, hi, o, i]`` is the [I, Q] block of the azimuthal mean of
  ARTS's laboratory-frame phase matrix, by the periodic midpoint rule at
  ``(k + 1/2) 360 / N`` deg, ``N = azimuth_count`` (even).  The samples above
  180 deg mirror those below, whose [I, Q] blocks are equal for the mirror
  symmetric media of ARTS's TRO and ARO formats, so ARTS is only asked for
  the azimuths in (0, 180) deg.  The midpoints avoid the principal plane,
  near which ARTS's rotation coefficients snap angles within about 1e-3 rad.
* ARTS's laboratory-frame phase matrix and RT4 use the same Q = I_v - I_h,
  and the [I, Q] block of the mean does not depend on the azimuth sense.

Accuracy:

* For a scattering matrix that is a regular Legendre series of degree L,
  the laboratory-frame matrix is a trigonometric polynomial of degree L in
  the azimuth difference, and the mean is exact for ``N > L`` (Rayleigh:
  ``N >= 4``); otherwise it converges as fast as the Fourier series of Z.
* The data are as accurate as ARTS's laboratory-frame phase matrix.
  ``GasScatterer`` and ``HenyeyGreensteinScatterer`` evaluate their
  closed-form scattering matrix at the exact scattering angle of every
  direction pair.  Particle habits interpolate linearly on their own
  scattering-angle grid, a ``ZenGrid`` (strictly ascending in [0, 180]
  deg).
* A vertical ray (``mu = 1``, e.g. the last Lobatto node or an extra
  angle) has the meridional plane of its azimuth label in ARTS's laboratory
  frame, the limit along its meridian, so vertical rays need no special
  treatment: for a pair of them the reference plane turns with the azimuth
  and the mean of Q is 0, as for any other pair.

**problem_from_path** follows the conventions of the DISORT workspace
methods (``disort_settingsOpticalThicknessFromPath``,
``disort_settingsLayerThermalEmissionLinearInTau``):

* ``ray_path``, ``atm_path`` and ``spectral_propmat_path`` have one entry per
  level, top first; the ``ray_path`` altitudes must decrease strictly and are
  the heights [m].  Only the altitudes of ``ray_path`` are used.
* The level temperatures are ``atm_path``'s.
* ``spectral_propmat_path`` is the gas propagation matrix only, per metre,
  with ``freq_grid.size()`` entries per level.  A layer's gas extinction is
  the mean of A at its two levels.  Polarized gas propagation matrices are
  rejected.
* The frequency is ``freq_grid[freq_index]``.
* Each layer gets the mean of its two levels' ``scattering_optics`` on RT4's
  quadrature for ``settings``, or is gas-only when that mean is all zero.
* ``settings``, ``ground``, ``surface_temperature`` and ``sky_temperature``
  go into the problem unchanged.

Tests (references external to the code under test):

* ``cpp.fast.rt4-arts-test`` (``src/core/polradtran/rt4/test/rt4-arts-test.cpp``):
  ARTS's Rayleigh ``GasScatterer`` through ``scattering_optics`` against
  ``sigma / (4 pi)`` times the m = 0 closed form
  (``P_II = 3/8 (3 - a - b + 3 a b)``, ``P_IQ = 3/8 (1 - 3 a)(1 - b)``,
  ``P_QI = 3/8 (1 - a)(1 - 3 b)``, ``P_QQ = 9/8 (1 - a)(1 - b)``,
  ``a = mu_out^2``, ``b = mu_in^2``) in all four quadrants, with the extra
  angles 0.35 and 1, for ``N = 4`` and 64: 5e-15 of sigma / (4 pi)
  (tolerance 1e-13); the path builder against its inputs and the level
  optics; error paths.
* ``cpp.fast.vdisort-arts-comparison``
  (``src/core/disort-cpp/test/vdisort/vdisort-arts-comparison.cpp``): RT4,
  RT3 and VDISORT through their path
  builders on an ARTS atmosphere (``AtmField`` with a temperature profile, a
  Rayleigh ``GasScatterer``, a cloud of 1.5 mm water spheres from ARTS's Mie
  code and gas absorption) at 89 GHz, thermal emission, Lambertian and
  Fresnel surfaces, at ``max_delta_tau = 1e-7``.  RT3 and RT4 run Evans'
  identical doubling, so RT4 - RT3 isolates the input routes
  (laboratory-frame azimuthal mean against RT3's Legendre series): 5.3e-10
  of max I.  RT4's layer phase matrices equal VDISORT's to 1.2e-15
  (relative), and RT4 - VDISORT is 5.4e-7 of max I, RT4's doubling error.
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

* **Scope.** Only [I] or [I, Q]: ``nstokes`` must be 1 or 2.  RT4 is
  thermal-only (no beam source), plane-parallel, and limited to azimuthally
  symmetric media; it returns only the m = 0 azimuthal mode.  With
  ``nstokes = 1``, the Q coupling of a polarizing surface or medium is
  dropped.
* **Mirror symmetry.** ``SYMMETRIC`` is hard-coded ``.TRUE.``.  The doubling
  discards the "-" reflection and transmission after the first step, so the
  medium must be mirror symmetric between the hemispheres:

  * ``extinction[down] == extinction[up]``;
  * ``absorption[down] == absorption[up]``;
  * ``phase[down, down] == phase[up, up]``;
  * ``phase[down, up] == phase[up, down]``.

  ``solve()`` rejects optics that break this by more than 1e-10 relative to
  the largest magnitude of each quantity.  Because RT4 requires
  ``phase[down, up] == phase[up, down]``, a test can never detect an
  exchange of those two quadrants.
* **Reentrant.** The C++ port keeps no state between calls, so concurrent
  ``rt4::solve`` calls run in parallel.
* **Lambertian with G or L quadrature.** The Lambertian surface conserves
  energy on the streams only with ``double_gauss``, where
  ``2 sum mu w = 1``.  With ``gauss`` or ``lobatto`` it is off by about
  3e-3 A for 8 streams, so use ``double_gauss``.
* **No array limits.** Every array is sized to the problem.  The Fortran's
  static arrays limited ``nstokes * nmu_total`` to 64, ``nlay`` to 400 and
  ``(nlay + 1) * (nstokes * nmu_total)^2`` to 301 * 4096; the port does not
  (``cpp.fast.rt4-test`` runs its isothermal Kirchhoff check with
  ``nstokes * nmu_total = 82`` and with 450 layers).  The preconditions
  above and the shapes are checked before the solve.
* **Zero pivot.** The inverse of ``1 - R R`` (LAPACK's, where RT4 had
  ``MINVERT``) throws for an exactly singular matrix, which is not expected
  for physical inputs.
* **Beam, m > 0 modes, U and V.** RT4 has none of them; RT3
  (:doc:`dev.rt3`) is the reference for those paths.

Mapping to VDISORT inputs
-------------------------

Feed both solvers the same layer optics, so that a comparison tests the
solvers and not two preprocessing chains.  VDISORT has no workspace layer
either.  Use the low-level ``pyarts3.arts.cppvdisort`` (or
``vdisort::main_data``) with ``NFourier = 1``.  VDISORT can represent only a
subset of RT4's inputs:

* K = k 1, the same for all streams;
* a = [a1, 0], the same for all streams;
* randomly oriented or Rayleigh-like scatterers whose m = 0 [I, Q] block
  does not couple to [U, V].

It cannot represent a direction-dependent or dichroic K(mu), such as an
off-diagonal K12, or a direction-dependent a(mu).  Such cases must be
rejected, not approximated.

For layer l with particle extinction ``k_e``, absorption ``a1``, scattering
``k_s = k_e - a1`` (assuming the optics conserve energy) and gas ``k_g``:

* **Quadrature.**

  * Use RT4 ``double_gauss`` with ``nmu = NQuad / 2``.  The nodes and
    weights are then the same as VDISORT's double-Gauss streams, to
    round-off, and both are ascending.
  * VDISORT streams ``0 .. N-1`` are upward (mu > 0) and ``N .. 2N-1`` are
    downward (mu = -mu[i]).  So VDISORT stream ``i`` is RT4 ``(up, i)``, and
    stream ``N + i`` is RT4 ``(down, i)``.
  * RT4's ``extra_mu`` correspond to VDISORT's user angles
    (``ungridded_u_user``): ``+mu`` for ``up`` and ``-mu`` for ``down``.
    Over a Fresnel surface an upward user angle ``+mu`` needs ``-mu`` as
    well: VDISORT reflects the downward user-angle radiance into it, and
    throws without it.
    ``user_phase[alpha, 0, l, u, j]`` is the same ``4 pi Z / k_s`` as below,
    taken from the RT4 ``phase`` row of the extra angle.  The RT4 phase rows
    of the extra angles must hold the physical phase matrix; their columns
    carry no weight.
* **Optical depth.** ``tau_arr`` holds the cumulative optical depth at the
  layer bottoms, ``sum (k_e + k_g) dz``.  RT4 level ``l`` corresponds to
  VDISORT ``tau = tau_l``, with ``tau_0 = 0``.
* **Output.** ``result.up[l, i]`` is stream ``i`` and ``result.down[l, i]``
  is stream ``N + i`` of VDISORT's ``u`` at ``tau_l``.
* **Single-scattering albedo.** ``omega = k_s / (k_e + k_g)``, and 0 for
  gas-only layers.
* **Phase matrix.**

  * ``phase_matrix[alpha=0, m=0, l, o, i][0:2, 0:2] = 4 pi Z / k_s``, with
    ``Z`` the RT4 ``phase`` quadrant that matches the hemispheres of VDISORT
    streams ``o`` (out) and ``i`` (in).
  * Do not include omega, weights or a ``2 - delta_m0`` factor.  VDISORT's
    ``(omega / 2) sum_j W_j P`` then equals RT4's ``2 pi sum_j w_j Z`` per
    unit length.
  * The [U, V] block, alpha = 1, does not affect [I, Q] at m = 0; it can be
    zero or the true [U, V] block.
* **Source.**

  * VDISORT emits ``(1 - omega) B`` per unit optical depth, which equals
    RT4's ``(a1 + k_g) B`` per unit length.
  * ``s_poly_coeffs[l] = [c0, c1]`` (Stokes ``[c, 0, 0, 0]``) in the global
    optical depth: ``c1 = (B_bottom - B_top) / dtau_l`` and
    ``c0 = B_top - c1 tau_top``.  This is the same linear-in-tau Planck
    function.
* **Top boundary.** ``b_neg[0, 0, i] = [B(T_sky), 0, 0, 0]``.
* **Lambertian surface.**

  * ``vdisort.lambertian_fourier_modes(A, 1)``, whose discrete operator
    ``(2 A / pi) pi W mu`` equals RT4's ``2 A mu w``;
  * ``b_pos[0, 0, i] = [(1 - A) B_s, 0, 0, 0]``.
* **Fresnel surface.**

  * ``vdisort.fresnel_fourier_modes(n, 1)``, which reproduces RT4's
    specular R at the nodes;
  * ``b_pos[0, 0, i] = [(1 - R1) B_s, -R2 B_s, 0, 0]``, using the same
    ``R1``/``R2`` as above, so Q > 0.
  * VDISORT drops Fresnel reflection at off-node user angles.  For an
    off-node user angle, ``fresnel_fourier_modes`` reflects nothing, and the
    emission is interpolated from the nodes.  So compare upward radiances at
    the quadrature streams only.  The comparison test measures this: with
    n = 3+0.2i and 8 streams per hemisphere, the missing ``R I_down`` is
    22% (mu = 0.35) and 13% (mu = 1) of max I.  The emission interpolation
    error is 7e-4 and 1.4e-3 of B_s.
* **Other surfaces.**

  * A ``DiscreteSurface`` with ``reflection[i, j] = pi w_j mu_j rho(mu_i,
    mu_j)`` corresponds to a custom ``vdisort::BDRF``.  Its cosine callback
    returns the BRDF ``rho`` (the [I, Q] block of R^0) at the cosines it is
    given, because VDISORT adds the ``pi W_j mu_j`` itself at m = 0.  The
    sine callback can be zero.
  * Use the same emission for both.  RT4 takes it in W m-2 Hz-1 sr-1, and
    VDISORT takes it as ``b_pos``.
* **nstokes = 1.** Put only ``4 pi Z_II / k_s`` in the M00 element and
  compare I.  With no I-Q coupling in the medium, a polarizing surface does
  not feed Q back into I.
* **Units.** Both solvers are linear in B.  Pass B_nu in W m-2 Hz-1 sr-1 to
  VDISORT to compare directly with the RT4 result.
* **Stokes basis.** Both use [I, Q] with Q = I_v - I_h and the meridional
  reference plane.

RT4 has the opposite limits: no U and V, no m > 0 modes, no beam source, and
mirror symmetry between the hemispheres is required.

``cpp.fast.vdisort-rt4-test``
(``src/core/disort-cpp/test/vdisort/vdisort-rt4-comparison.cpp``) implements
this mapping.  It compares I and Q at every level, stream and direction for
the cases listed in that directory's ``COVERAGE.md``.  Measured results:

* Gas-only atmospheres agree to 5e-15 relative to max I.
* Otherwise, RT4's first-order doubling error dominates.  At
  ``max_delta_tau = 1e-7`` the difference is 0.1 to 0.7 times the
  initial-layer thickness, at most 4e-8.  It halves exactly when
  ``max_delta_tau`` halves, and does not grow with the number of streams.
* Against RT4 Richardson-extrapolated to ``max_delta_tau = 0``, the
  difference is at most 3.4e-10.
