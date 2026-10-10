RT3 reference solver
====================

RT3 is Evans' polarized doubling-adding solver for randomly oriented
particles, from the polradtran package.  ARTS 3 keeps it only as an external
reference for validating other solvers, in particular the solar-beam, m > 0
and U, V paths of the polarized discrete-ordinate solver VDISORT, which RT4
(:doc:`dev.rt4`) cannot check.  RT3 solves the plane-parallel radiative
transfer equation with a direct (solar) beam and thermal sources, for every
Fourier azimuth mode and the Stokes components [I], [I, Q], [I, Q, U] or
[I, Q, U, V].  It has a core C++ interface, the namespace ``polradtran::rt3``
(``src/core/polradtran/rt3``; ``rt3::`` below), and a low-level Python
interface (``pyarts3.arts.rt3``).  It has no workspace methods, variables or
agendas.  Its equations are in :doc:`concept.rt3`.

Provenance
----------

* **Original code.** K. F. Evans, polradtran (RT3/RT4), distributed from
  https://nit.coloradolinux.com/polrad.html under the MIT licence
  (``3rdparty/polradtran/LICENSE``).  RT3 is described by
  :cite:t:`Evans1999`.  ARTS 2 never shipped RT3.
  Of ``PolRadTran.tar`` (sha256 ``7b0eff79...a6f7cff9d``),
  ``3rdparty/polradtran`` keeps the licence, ``README``, the four test
  scripts and ``cl340d14.dda``, unchanged.  The Fortran was removed once the
  port to C++ was complete; the history of ``README.ARTS`` there lists the
  ARTS changes to it.
* **ARTS 3 changes** to the Fortran, each marked ``c ARTS3:`` in the source,
  before it was ported:

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
    link; ``librt3`` and ``librt4`` exported no common symbol.
  * Planck constants computed from the exact SI h, c and k, as for RT4.  The
    original 5-digit constants give a bias of -3e-5 in the Rayleigh-Jeans
    limit, growing to -2e-4 and -3e-4 at 3 um and 300 and 200 K.
  * A new ``rt3_c_interface.f90`` with ``ISO_C_BINDING`` entry points.

* **Port to C++.** RT3 is ported to C++ one routine at a time, as RT4 was
  (:doc:`dev.rt4`), on matpack and rtepack types with the method unchanged,
  sharing no code with the RT4 port while it lasts (so that every step
  compares with RT3's own Fortran).  The port is not bitwise: the C++ is
  written in its natural order and the compiler may contract multiply-adds
  into FMAs, so results differ from the Fortran by rounding.  Every routine
  that ``RADTRAN`` calls is ported, and ``rt3::solve`` calls no Fortran;
  over the whole port, the 148 capture problems changed by at most 4.7e-13
  of the m = 0 I (in optically thick layers, through LAPACK's inverse).
  The Fortran, and with it the lock that serialised ``rt3::solve``, was
  then removed.  The steps:

  * ``RADTRAN`` is ``rt3::radtran`` (``src/core/polradtran/rt3/radtran3.h``).  It
    follows the Fortran step by step and calls the same subroutines, all
    now C++.  The Fortran ones were called through entry points in
    ``rt3_c_interface.f90``, as the tests called the Fortran until it was
    removed.
    Evans' matrix helpers are matpack (``MZERO`` is ``= 0.0``,
    ``MIDENTITY`` ``matpack::identity``, ``MCOPY`` ``=``, ``MSCALARMULT``
    ``*=``).  It takes no counts (``NSTOKES``, ``NUMMU``, ``AZIORDER``,
    ``NUM_LAYERS``, ``NSL``, ``LDCOEF`` and ``NOUTLEVELS`` are the extents
    of its arrays), the extra angles as an input of their own instead of
    ``QUAD_TYPE 'E'`` and non-zero ``MU_VALUES``, the ground as data (see
    below), and ``polradtran::quadrature_type``.  Its work arrays are sized to the
    problem, including the 210 MB static scattering-matrix buffer, and are
    those of an ``rt3::rt3_workdata`` (see below); its ``STOP`` checks
    throw.  It does not print RADTRAN's message when it truncates a
    Legendre series (``rt3::solve`` rejects a truncation that drops a
    non-zero coefficient).
  * The quadratures are ARTS's: ``polradtran::get_quadrature`` (used by
    ``rt3::radtran`` and ``rt3::solve``) is the positive half of
    ``scattering::DoubleGaussQuadrature``, ``GaussLegendreQuadrature`` or
    ``LobattoQuadrature`` of degree ``2 nmu``.  They are RT3's rules.
    Against 50-digit references for nmu up to 64, RT3's Gauss and
    double-Gauss weights are off by up to 2.3e-12 relative and ARTS's by
    3e-16; the Lobatto rules are equally accurate (4.5e-14).  The radiances
    of the 148 capture problems changed by at most 4.3e-14 of the m = 0 I
    at the same level and stream (median 2.7e-15).  In
    ``cpp.fast.rt3-radtran-test`` one optically thick, strongly scattering
    case changes by 1.7e-11 of the largest radiance: its smallest node
    moves by 3 ulp, which flips a pivot of LINPACK's ``DGEFA`` in the
    doubling, where the inverse is ill-conditioned.
  * All of ``radscat3.f`` is C++ (``src/core/polradtran/rt3/radscat3.h``, with the FFT
    in ``rt3_fft.h``); it calls no Fortran and keeps no static state.  The counts are the extents of their
    arrays.  A Legendre series is a ``CompactPlanarMuelmatVector``, each
    coefficient the six elements of a scattering-plane phase matrix by name
    (not in the column order of RT3's files), the phase matrix in the
    scattering plane is a ``CompactPlanarMuelmat``, and the phase matrices in
    the meridional planes are ``Muelmat`` (element (r, c) of the Fortran
    matrix is ``[r - 1, c - 1]``), of which only the leading
    ``nstokes x nstokes`` is transformed (the rest is 0):

    * ``rt3::get_scat_set`` (``GET_SCAT_SET``) throws where the delta-M
      scaling divides by zero (an extinction that is not positive,
      ``1 - f = 0`` or ``1 - albedo f = 0``), where the Fortran returned NaN
      or infinity.  Its delta-M scaling of a coefficient is
      ``rt3::delta_m_scaled``, which ``rt3::solve`` also uses to check the
      scaled series.
    * ``rt3::scattering`` (``SCATTERING``) writes the part of ``SCATBUF`` of
      one set, a ``MuelmatTensor4 [aziorder + 1, 2, nummu, nummu]``, each
      mode straight into it.  The Fortran's limits (``FFT1DR``
      took at most 512 azimuths, ``FOURIER_MATRIX`` 1024) were the sizes of
      its buffers, which the port does not have.
    * ``rt3::direct_scattering`` (``DIRECT_SCATTERING``) writes the part of
      ``DIRECTBUF`` of one set, a ``StokvecTensor3 [aziorder + 1, 2, nummu]``: the
      first column of each mode of the phase matrix from the sun's
      direction (the cosine modes of I and Q, the sine modes of U and V).
      The Fortran's limit of 512 azimuths and modes is not the port's.
    * ``rt3::get_scattering`` (``GET_SCATTERING``) copies the leading
      ``nstokes x nstokes`` of one mode of a set's ``SCATBUF`` part into
      ``SCATTER_MATRIX``
      (``[4, nummu, nstokes, nummu, nstokes]``), and
      ``rt3::scatter_symmetry`` (``SCATTER_SYMMETRY``) makes P-- and P-+
      from it, copying the diagonal 2 x 2 Stokes blocks and negating the
      others.
    * ``rt3::get_direct`` (``GET_DIRECT``) copies the leading ``nstokes`` of
      one mode of a set's ``DIRECTBUF`` part, ``[2, nummu, nstokes]``.
    * ``rt3::check_norm`` (``RT3_CHECK_NORM``) throws where the Fortran
      stopped (the I-I term not integrating to 1 within 1e-7), and also for
      NaN, which the Fortran let pass.
    * ``rt3::sum_legendre`` (``SUM_LEGENDRE``) sums the series of compact
      matrices, or that of F11 alone for nstokes 1.  That is ``NUMBER_SUMS``'s
      choice where it matters (it also skipped F34 for nstokes 2 and 3, which
      the rotation does not mix into the leading 3 x 3, and took F22 and F44
      from F11 and F33 when they were equal), so ``NUMBER_SUMS`` is gone.  It
      sums with ARTS's Legendre
      polynomials (``Legendre::legendre_polynomials``, Boost's recurrence),
      the generator of every Legendre series in ARTS, made once for all
      six series where RT3 ran its own recurrence for each.  Against
      50-digit references, for series (2l + 1) g^l up to degree 1023 and
      cosines clustered at +-1, both recurrences err by at most 2e-13 of
      :math:`\sum_l|c_l|`, ARTS's by half as much on average (7e-16 against
      1.3e-15).  The cosine is clamped to [-1, 1], which it can leave by
      rounding (one ulp at 1 moves P_1023 by 1e-10).  Against the Fortran
      it agrees to 1.1e-16 of the largest value; the 148 capture problems
      changed by at most 1.1e-15 of the m = 0 I.
    * ``rt3::rotate_phase_matrix`` (``ROTATE_PHASE_MATRIX``) finds RT3's
      rotation angles and rotates the compact matrix with
      ``rtepack::rotated``, the closed form of the two Stokes rotations
      around it.
    * ``MATRIX_SYMMETRY``, which negates the off-diagonal 2 x 2 blocks, is
      ``rtepack::mirror``.
    * ``rt3::fourier_matrix`` and ``rt3::fourier_basis`` (``FOURIER_MATRIX``,
      ``FOURIER_BASIS``, with ``rt3::fourier_direction`` in place of the
      sign of ``DIRECTION``); the basis order is passed, as the basis has
      ``order + 1`` or ``2 order + 1`` elements, and so is nstokes.
    * ``rt3::combine_phase_modes`` (``COMBINE_PHASE_MODES``): its
      ``SINFLAG`` table is the block structure of ``rtepack::mirror``.

  * ``RT3_THERMAL_RADIANCE`` is ``polradtran::thermal_radiance``
    (``src/core/polradtran/radutil.h``, shared with RT4), with ARTS's
    ``planck()`` in SI in place of RT3's ``PLANCK_FUNCTION``; it throws for
    a negative temperature, where RT3 gave 0.  ``planck()`` evaluates
    ``expm1``: against 50-digit references it is within 1.7e-15, RT3's
    ``exp(x) - 1`` within 2.8e-13 (at small h nu / k T).  With it
    ``rt3::radtran`` works in SI (W m-2 Hz-1 sr-1) at the frequency, and
    ``rt3::solve`` no longer converts (the Fortran ground radiances, per
    micrometre, were converted until they were ported).  The 148 capture
    problems changed by at most 7.4e-15 of the m = 0 I.
  * The ground is an input of ``rt3::radtran``, as RT4's is of
    ``rt4::radtrano``, in place of ``GROUND_TEMP``, ``GROUND_TYPE``,
    ``GROUND_ALBEDO`` and ``GROUND_INDEX``.  Both grounds of RT3 make the
    same surface layer (no reflection from above, the identity as
    transmission, no source); only the reflection back up depends on the
    ground and the azimuth mode.  ``radtran`` takes, for every mode,
    that reflection (``surf_reflect``,
    ``[aziorder + 1, nummu, nstokes, nummu, nstokes]``), the ground's own
    radiance (``gnd_radiance``, ``[aziorder + 1, nummu, nstokes]``) and,
    unlike RT4, which has no beam, the radiance reflected from the direct
    beam per unit of direct flux (``direct_reflect``, sr-1, same shape).
    The direct flux that reaches the ground is computed inside
    ``radtran``, so the ground cannot add that part itself: ``radtran``
    makes ``GND_RADIANCE = gnd_radiance + F direct_reflect`` with the
    solar source.  ``polradtran::external_surface_layer`` makes the surface layer
    (RT4's ``EXTERNAL_SURFACE``; RT3 has none).
    ``rt3::ground_surface`` (``src/core/polradtran/rt3/radutil3.h``) makes the three
    arrays from an ``rt3::surface``: the Lambertian ground reflects and
    emits in mode 0 only, emits only with the thermal source, and reflects
    the beam as ``A / pi``.  The Fresnel ground reflects the same in
    every mode, always emits (in mode 0, as in RT3), and throws with the
    solar source, because RT3 cannot reflect the beam specularly.  It
    calls no Fortran (the ground routines are ported, below).  Making the
    ground an
    input changed the 148 capture problems by at most 4.3e-16 of the
    m = 0 I (the beam's reflection is added in SI).
  * ``RT3_LAMBERT_SURFACE`` is ``polradtran::lambert_surface_layer`` (named, as
    in the RT4 port, for the layer it makes, since
    ``polradtran::lambertian_surface`` is the ground's type): ``2 A mu_j w_j``
    from stream j into every stream, I to I only, in mode 0, as joker
    slices of a ``[2, nummu, nstokes, nummu, nstokes]`` view.
    Bit-identical.
  * ``RT3_LAMBERT_RADIANCE`` is ``polradtran::lambert_radiance``, in SI with
    ARTS's ``planck()``: in mode 0 the emission ``(1 - A) B`` with the
    thermal source and the reflected beam ``F A / pi`` with the solar
    source.  It throws for a negative temperature where it uses it (RT3
    gave 0).  Against the Fortran it differs as the Planck functions do
    (2.8e-13); the capture changed by at most 3.3e-15 of the m = 0 I.
  * ``RT3_FRESNEL_SURFACE`` and ``RT3_FRESNEL_RADIANCE`` are
    ``polradtran::fresnel_surface_layer`` and ``polradtran::fresnel_radiance``, as in
    the RT4 port: ARTS's ``fresnel()`` amplitudes at ``acos(mu)`` and
    ``rtepack::fresnel_reflectance`` for the Mueller matrix (RT3's
    ``R``, with ``R(U, V) = -R4`` and ``R(V, U) = R4``), and the emission
    ``(1 - R) B`` with ``planck()``.  The layer agrees with the Fortran to
    1.4e-15 of its largest value, the radiance as the Planck functions do
    (6.1e-13 at 1 GHz and 150 K); the capture changed by at most 3.2e-15
    of the m = 0 I.  With these, ``rt3::radtran`` uses nothing of
    ``radutil3.f`` (its Planck function and quadratures are ARTS's).
  * ``RT3_NONSCATTER_LAYER`` is ``polradtran::nonscatter_layer``
    (``src/core/polradtran/radintg.h``, shared with RT4): the
    reflection, transmission and source of a purely absorbing layer, the
    source in mode 0 only.  ``radtran`` passes it
    ``[2, nummu, nstokes, nummu, nstokes]`` and ``[2, nummu, nstokes]``
    views of its layer arrays.
  * ``RT3_INITIALIZE`` is ``rt3::initialize``: the thin starting layer's
    reflection and transmission from the phase function, extinction and
    albedo, as row slices of ``[2, nummu, nstokes, nummu, nstokes]`` views.
    The diagonal of the transmission is kept in the Fortran's form,
    ``1 - f (1 - albedo P)`` rounded once: the layer's extinction ``f`` is
    as small as ``max_delta_tau`` (1e-6), so the diagonal's last bit is a
    relative 1e-10 of it, which the doubling carries to the radiances
    (computed as ``(1 - f) + f albedo P``, rounded twice, some
    ``cpp.fast.rt3-radtran-test`` cases moved by up to 9e-10).
  * ``RT3_DOUBLING_INTEGRATION`` is ``polradtran::doubling_integration``,
    which RT4 shares (with the linear source alone): the products are BLAS
    ``mult`` (DGEMM, and DGEMV for the source vectors) whose alpha and beta
    absorb the ``MIDENTITY``, ``MSUB`` and ``MADD`` around them, and
    ``MINVERT`` is LAPACK's ``inv_inplace``; it doubles RT3's exponential
    (solar) source as well as the linear (thermal) one.  Its ``T_EXP``,
    the doubled solar source of a step, is held in the output ``t_source``,
    which is written only after the doubling, so the work data needs no
    array that RT4 would not use.  Each doubling about squares T, doubling
    its relative error, and inverts ``1 - R R``, whose condition number
    kappa amplifies the inverse's rounding, so n doublings amplify rounding
    by 2^n kappa.  Against the Fortran it agrees to 0.8 epsilon times
    2^n kappa: 7.6e-13 of the largest value on Apple arm64 with OpenBLAS,
    1.8e-10 (20 doublings) on AMD x86_64 with MKL.  The 148 capture
    problems changed by at most 4.5e-13 of the m = 0 I (in optically thick
    layers, through the inverse).
  * ``RT3_COMBINE_LAYERS`` is ``polradtran::combine_layers``, which puts one
    layer on top of another in the adding, written as
    ``polradtran::doubling_integration`` (BLAS ``mult`` with alpha and beta,
    LAPACK's ``inv_inplace``).  It is RT4's ``COMBINE_LAYERS``, line by
    line, so both ports share it (below).
    On the layers ``RADTRAN`` combines (thin and thick scattering layers,
    gas, the Lambertian and Fresnel grounds) it agrees with the Fortran to
    6.8e-16 of the largest value; the 148 capture problems changed by at
    most 2.9e-15 of the m = 0 I.
  * ``RT3_INTERNAL_RADIANCE`` is ``polradtran::internal_radiance``, the
    radiances at a level from the atmosphere above and below it, written
    as ``polradtran::combine_layers``.  It agrees with the Fortran to 3.5e-16 of
    the largest value; the 148 capture problems changed by at most
    6.2e-16 of the m = 0 I.  With it ``rt3::radtran`` calls no Fortran.
  * Every work array of ``rt3::radtran`` and the routines it calls is in
    ``rt3::rt3_workdata`` (``src/core/polradtran/rt3/rt3_workdata.h``): the
    arrays it shares with RT4 (``polradtran::workdata``, below) and its own,
    the scattering sets and the layers' optics, the
    scratch of ``SCATTERING``, ``DIRECT_SCATTERING`` and
    ``FOURIER_MATRIX`` with ``FFT1DR``'s table, the layers and the ground,
    ``RADTRAN``'s arrays on the streams, ``DOUBLING_INTEGRATION``'s
    sources, and one ``X``, ``Y``, ``GAMMA``, two vectors and LAPACK's
    workspace for the inverse, shared by the doubling and the adding.
    ``radtran`` takes it last and sizes it, so a repeated call with the
    same sizes (as over frequency) allocates nothing; the scattering
    routines size their own scratch, as only they know how many azimuths
    they sample, and the doubling and adding throw if it is not sized for
    their streams.  Its values between calls are unspecified, except the
    FFT's table.  ``cpp.fast.rt3-radtran-test`` runs every case again with
    one work data shared over all cases and with one sized for the case
    whose every array is NaN: neither changes a bit.  ``rt3::solve`` keeps
    one per call.
  * What RT3 and RT4 share is one C++, in the namespace ``polradtran``
    (``src/core/polradtran``, the library ``arts_polradtran``): the
    routines that are the same Fortran in both, line by line
    (``RT3_COMBINE_LAYERS``, ``RT3_INTERNAL_RADIANCE``,
    ``RT3_LAMBERT_SURFACE`` and ``RT3_FRESNEL_SURFACE`` are RT4's
    routines without the prefix), those of which RT4's is RT3's in the
    azimuth mode 0 with the thermal source alone
    (``RT3_DOUBLING_INTEGRATION``, ``RT3_NONSCATTER_LAYER``,
    ``RT3_THERMAL_RADIANCE``, ``RT3_LAMBERT_RADIANCE`` and
    ``RT3_FRESNEL_RADIANCE``), RT4's ``EXTERNAL_SURFACE``, and
    ``polradtran::workdata``, RT4's work data, which ``rt3::rt3_workdata``
    extends.  Both ``cpp.fast.rt3-radtran-test`` and
    ``cpp.fast.rt4-radtrano-test`` compare it with their own Fortran, and
    sharing it left every result of both ports bit-identical.
    The level loop of ``RADTRAN`` and ``RADTRANO``, which adds the layers
    above and below a level and calls ``INTERNAL_RADIANCE``, is
    ``polradtran::level_radiance``, and their initial sublayer of a
    scattering layer and its number of doublings
    ``polradtran::initial_sublayer``.  The interfaces share their streams
    and grounds (``polradtran.h``: ``quadrature_type``, ``quadrature``,
    ``get_quadrature``, ``lambertian_surface`` and ``fresnel_surface``, in
    Python ``pyarts3.arts.polradtran``) and the layers of an ARTS
    propagation path (``polradtran::layers_from_path`` in
    ``polradtran_arts.h``, behind both ``problem_from_path``).
  * ``RT3_INITIAL_SOURCE`` is ``rt3::initial_source``: the thin starting
    layer's source, delta_z / mu times the extinction times a source vector
    (the solar pseudo source or the thermal one), per angle.
  * ``FFT1DR``, Evans' real FFT with ``FFTC``, ``FIXREAL`` and
    ``MAKEPHASE``, is ``rt3::fft1dr`` in a file pair of its own
    (``src/core/polradtran/rt3/rt3_fft.h``).  It is RT3's FFT and the default.  A
    build may use FFTW instead, but only as a compile-time option: FFTW's
    license keeps it out of the default build.  The rest of RT3 uses only
    ``fft1dr``, ``fft_direction`` and ``fft_workdata``.  The header gives
    the packed format and scaling that an FFTW ``fft1dr`` must also give
    (FFTW's r2c followed by a complex conjugate, and a conjugate followed
    by c2r), and ``cpp.fast.rt3-radtran-test`` checks it against direct
    sums.  ``FFT1DR``'s SAVEd phase table is ``fft_workdata``, which the
    caller owns and passes down (it is part of ``rt3_workdata``; FFTW would
    keep its plans there).  Its limit of 512 values is kept (``FFT1DR``'s
    ``STOP`` throws).  ``MAKEPHASE``'s table fits in ``4 nmax`` values only
    for ``nmax`` a power of two, as ``FFT1DR`` uses it; for another
    ``nmax`` the Fortran writes past it, so ``rt3::makephase`` throws.

  The port was tested against the Fortran ``RADTRAN``, built beside it, on
  random inputs over every branch of ``RADTRAN`` (178 cases).
  Mathematically equivalent evaluations (FMA contraction, vectorised libm
  functions, other BLAS kernels) are accepted: the tolerances allow 16
  epsilon times what a computation amplifies rounding by.  Every output had
  to agree to 1e-10 of its largest value, for the replaced quadratures and
  Planck function, plus 16 epsilon times 2^n for the layer doubled n times
  most (the largest difference was 1.2e-9 on AMD x86_64 with MKL), and
  ``polradtran::doubling_integration`` to 16 epsilon times 2^n kappa.
  ``polradtran::get_quadrature``, ``rt3::ground_surface`` and each ported
  routine were checked against the Fortran: to rounding (1e-14 of the
  largest value; 1e-13 for ``rt3::scattering`` and
  ``rt3::direct_scattering``, whose Legendre sums amplify the rounding of
  the scattering angle), or exactly for the integer and copying ones
  (``number_sums``, ``matrix_symmetry``, ``get_scattering``,
  ``scatter_symmetry``, ``get_direct``).  Each step was also measured
  against the step before.  Since the FFT and the natural order of
  operations, the 148 capture problems differed from the all-Fortran
  ``SCATTERING`` by at most 5.4e-16 of the m = 0 I.  When the Fortran was
  removed, its outputs for twelve cases that cover every branch (1 to 4
  Stokes components, the three quadratures and an extra angle, azimuth
  orders, delta-M, every source code, both grounds, gas, thin, thick and
  shared layers, the three summation cases and a coarse ``max_delta_tau``)
  were kept as constants
  (``src/core/polradtran/rt3/test/rt3-radtran-reference.h``), with which
  ``cpp.fast.rt3-radtran-test`` compares the port to the same tolerances.

  ``rt3.f``, Evans' original main program, built from the tar, ran his two
  RT3 scripts and reproduced his tables, every value to one unit in its
  last printed digit, before it was removed with the Fortran.
  ``OUTPUT_FILE`` also showed how RT3 sums its Fourier series.

Build
-----

RT3 is C++ and always built, as RT4 is (:doc:`dev.rt4`).  Its arrays are
sized to the problem: the 230 MB of static arrays of the Fortran (its
scattering-matrix buffer ``SCATBUF`` alone was 210 MB) are gone.

Tests
-----

There are two tests, ``cpp.fast.rt3-test``
(``src/core/polradtran/rt3/test/rt3-test.cpp``) and
``tests/core/rt3/closed-form.rt3.py``.  Every reference is either Evans' own
benchmark output of the original program or a closed form derived in the
test; none is an output of this build.  The C++ test covers:

* **(a) Quadrature exactness** of the three rules (moment errors 1e-14).
* **(b) runmietest**, the Mie case of :cite:t:`Evans1999`: tau = 1,
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
* **(h) 20 error paths.**
* **(i) Extra angles** with Lobatto and double-Gauss leave the radiances on
  the quadrature nodes unchanged, and an extra angle on a node reproduces
  it (to rounding times 2^n).
* **(j) Isothermal Kirchhoff**, I = B and Q = U = V = 0 in an isothermal
  Rayleigh medium over a Fresnel ground, with N = 68 and aziorder 1 and
  with 420 layers and 210 scattering sets, beyond the Fortran's array
  sizes.

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

The comparison of VDISORT against RT3 is ``cpp.fast.vdisort-rt3-test``.  See
`Mapping to VDISORT inputs`_.  Its
part E solves the problems of Evans' two scripts, read from the scripts, with
VDISORT; see ``src/core/disort-cpp/test/vdisort/COVERAGE.md``.  The
tests of the inputs from ARTS data are listed in `Inputs from ARTS data`_.

Interface
---------

C++ (``#include <rt3.h>``, namespace ``polradtran::rt3``, with the streams
and the grounds it shares with RT4 in ``polradtran.h``, namespace
``polradtran``):

.. code-block:: cpp

  // polradtran.h, namespace polradtran (shared with RT4)
  enum class quadrature_type { gauss, double_gauss, lobatto };   // 'G' (RT3's 'E'), 'D', 'L'
  struct quadrature { Vector mu; Vector weights; };
  quadrature get_quadrature(Index nmu, quadrature_type type);
  struct lambertian_surface { Numeric albedo; };
  struct fresnel_surface { Complex refractive_index; };

  // rt3.h, namespace polradtran::rt3
  Index max_legendre_degree(Index nmu, quadrature_type type);    // RT3's NLEGLIM
  struct scattering_set { Numeric extinction; Numeric scattering; CompactPlanarMuelmatVector legendre; };  // [nleg + 1]
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

Python (``pyarts3.arts.rt3``, with what RT4 shares in
``pyarts3.arts.polradtran``: ``QuadratureType``, ``get_quadrature(nmu,
type)``, ``Quadrature``, ``LambertianSurface(albedo)`` and
``FresnelSurface(refractive_index)``) mirrors this:
``max_legendre_degree(nmu, type)``, ``ScatteringSet(extinction, scattering,
legendre)`` and ``ArrayOfScatteringSet``, ``Problem(...)`` (every field as a
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
  rayleigh = np.array([[1.0, -0.5, 1.0, 0.0, 0.0, 0.0],    # F11 F12 F22 F33 F34 F44, l = 0
                       [0.0, 0.0, 0.0, 1.5, 0.0, 1.5],     # l = 1
                       [0.5, 0.5, 0.5, 0.0, 0.0, 0.0]])    # l = 2
  p = rt3.Problem(nstokes=4, nmu=8, aziorder=2,
                  direct_flux=1.0, direct_mu=0.6, thermal=False,
                  frequency=6e14,
                  height=[1.0, 0.0], temperature=[0.0, 0.0], gas_extinction=[0.0],
                  scattering_sets=[rt3.ScatteringSet(0.1, 0.1, rayleigh)],
                  layer_scattering_index=[0],
                  ground=arts.polradtran.LambertianSurface(0.0))
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
  convention of :cite:t:`Evans1999`, in which RT3 reproduced the V of
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
* ``extra_mu`` is allowed with every rule (RT3's 'E' type had them with
  ``gauss`` only).  The extra angles receive scattering and reflection but
  carry no weight, so they do not change the radiances on the quadrature
  nodes.
* ``down`` is radiation propagating downward, toward increasing optical
  depth (RT3's "+", printed with mu > 0 by rt3.f); ``up`` propagates upward
  (RT3's "-", printed with mu < 0).
* :doc:`concept.rt3` labels the hemispheres as DISORT does, "+" upward, and
  names the operators of a slab by their outgoing and incident hemispheres.
  The code's (Evans') ``R``, ``T`` and ``S`` of the "+" part (``[0]``) are
  its :math:`R^{-+}`, :math:`T^{--}` and :math:`S^-`, and those of the "-"
  part (``[1]``) its :math:`R^{+-}`, :math:`T^{++}` and :math:`S^+`.  The
  same holds for RT4.

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
* ``legendre`` is ``[nleg + 1]``: the scattering-plane phase matrix

  ::

    [[F11, F12,   0,   0],
     [F12, F22,   0,   0],
     [  0,   0, F33, F34],
     [  0,   0, -F34, F44]]

  as a plain Legendre series in cos(Theta),
  ``F = sum_l legendre[l] P_l(cos(Theta))``, each coefficient a
  ``CompactPlanarMuelmat`` with its elements by name.  As an array (in
  Python) it is ``[nleg + 1, 6]`` in rtepack's order F11, F12, F22, F33,
  F34, F44; RT3's scattering files have the columns F11, F12, F33, F34,
  F22, F44, and must be reordered.
* The basis is that of the scattering plane with Q = I_par - I_perp, so
  Rayleigh scattering has F12 = -3/4 sin^2(Theta) (the ``rayleigh`` array
  in the example).  The coefficients include the factor 2 l + 1
  (Henyey-Greenstein is ``(2 l + 1) g^l``), and the phase function is
  normalised to 1 over 4 pi: ``legendre[0].F11()`` must be 1.
* RT3 sums only F11 for nstokes 1.  Trailing all-zero coefficients are
  dropped before the call.
* With ``delta_m``, every set is scaled with M = 2 nmu_total (including the
  extra angles, because RT3 passes its NUMMU): f = legendre[M].F11() /
  (2 M + 1), extinction (1 - omega f) k, albedo (1 - f) omega / (1 - omega f),
  series ``(2 l + 1) (c_l / (2 l + 1) - f id) / (1 - f)`` with ``id`` the
  identity (so F12 and F34 become ``c_l / (1 - f)``), degree M - 1.  The beam is
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

* ``polradtran.LambertianSurface`` (RT3 'L'): reflection ``2 A mu_j w_j`` of the m = 0
  mode into every stream, I to I only; direct-beam reflection
  ``A F_direct(surface) / pi``.  The diffuse reflection conserves energy on
  the streams only with ``double_gauss`` (``2 sum mu w = 1``).
* ``polradtran.FresnelSurface`` (RT3 'F'): specular reflection for every mode under a
  medium of index 1, ``[[R1, R2, 0, 0], [R2, R1, 0, 0], [0, 0, R3, -R4],
  [0, 0, R4, R3]]`` with ``R1 = (|r_v|^2 + |r_h|^2) / 2``,
  ``R2 = (|r_v|^2 - |r_h|^2) / 2``, ``R3 = Re(r_v r_h*)``,
  ``R4 = Im(r_v r_h*)``; emission ``[(1 - R1) B, -R2 B, 0, 0]``.  Without a
  beam the field is azimuthally symmetric, so R3 and R4 never act.

**Fortran buffers** (for maintainers): ``SCAT_COEF(6, LDCOEF, set)`` is a
``CompactPlanarMuelmatMatrix [set, LDCOEF]`` that the ``legendre`` of each
set leads; ``SCATLAYERS`` is the 1-based set or 0;
``OUTLEVELS`` lists every level; ``UP_RAD``/``DOWN_RAD(s, mu, m, level)`` is
exactly the row-major ``[level, m, mu, s]`` of the result, and
``UP_FLUX``/``DOWN_FLUX(s, level)`` the row-major ``[level, s]``.  For the
'E' type the first ``nmu`` entries of ``MU_VALUES`` were passed as 0.

Inputs from ARTS data
---------------------

``src/core/polradtran/rt3/rt3_arts.h`` builds RT3 inputs from ARTS scattering species,
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
* The elements are ARTS's, by name, without any sign change, and the series
  is normalised by ``c_0.F11()`` so that ``legendre[0].F11() = 1``.
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
case of :cite:t:`Evans1999` and of Garcia and Siewert (1989): spheres
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

* ``cpp.fast.rt3-arts-test`` (``src/core/polradtran/rt3/test/rt3-arts-test.cpp``):
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
* **Normalisation.** ``legendre[0].F11()`` must be 1 to 1e-9 after delta-M
  scaling.  RT3's ``CHECK_NORM`` stops the process when the discrete
  normalisation is off by more than 1e-7; for a series within NLEGLIM that
  discrete normalisation equals ``legendre[0].F11() - 1`` up to round-off.
* **No array limits.** Every array is sized to the problem.  The Fortran's
  static arrays limited N = nstokes * nmu_total to 64, the layers and the
  scattering sets to 200 each, the scattering-matrix, direct-beam and
  azimuth-basis buffers (2 A + 1 to 512 with a beam and 1024 without), the
  Legendre degree to 1023, and with A > 0 the degree it sums to 251 (its
  FFT held 512 azimuth samples); the port has none of them.
  ``cpp.fast.rt3-test`` runs its isothermal Kirchhoff check with N = 68 and
  with 420 layers and 210 sets, and ``cpp.fast.rt3-radtran-test`` checks
  ``fft1dr`` up to 4096 values and the m = 0 mode of ``SCATTERING`` and
  ``DIRECT_SCATTERING`` with 1024 and 4096 azimuths (degree 300 and 1100).
  Every scattering set is precomputed, also an unused one, so memory grows
  as sets * (A + 1) * 2 N^2.
* **Accuracy.** The initial doubling sublayer is first order in its slant
  thickness ``max_delta_tau / mu_min``.  Its transmission is stored as
  1 - O(max_delta_tau), so results carry a round-off floor of about
  eps / max_delta_tau (2e-10 at the default 1e-6).  The thermal source of a
  scattering layer is linear in optical depth.
* **Cost.** For each output level, RT3 adds all layers above and below anew,
  so the cost grows as nlay^2 per azimuth mode; ``solve()`` returns every
  level.
* **Reentrant.** The C++ port keeps no state between calls, so concurrent
  ``rt3::solve`` calls run in parallel.
* **Zero pivot.** The inverse of ``1 - R R`` (LAPACK's, where RT3 had
  ``MINVERT``) throws for an exactly singular matrix, which is not expected
  for physical inputs.
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
