DISORT implementation and validation
=====================================

The scalar and vector solvers are implemented in ``src/core/disort-cpp``.
See :doc:`concept.disort` for the equations and notation used below.

Implementation array layouts
============================

The similarly named arrays in the two cores do not always store the same
factorization:

.. list-table::
   :header-rows: 1
   :widths: 24 35 41

   * - Mathematical object
     - CPP-DISORT
     - VDISORT
   * - Layer transport matrix :math:`\boldsymbol A`
     - Implicit in ``D_pos``, ``D_neg``, ``apb``, and ``amb``; only the
       :math:`N\times N` reduced product ``sqr`` is diagonalized
     - The :math:`4N_q\times4N_q` real ``transport`` is diagonalized in real
       arithmetic, and each conjugate pair is expanded to two complex columns
   * - Propagation constants :math:`\boldsymbol K`
     - ``K_collect[m, layer, eigen]`` with negative half followed by its
       positive partners
     - ``K_collect[alpha, m, layer, eigen]`` after complex lexicographic sort
   * - Eigenvectors :math:`\boldsymbol G`
     - ``G_collect[m, layer, state, eigen]``
     - ``G_collect[alpha, m, layer, state, eigen]``
   * - Modal constants :math:`\boldsymbol c_\ell`
     - The band-solve ``RHS`` is overwritten by the constants and copied to
       ``C_collect[m, layer, eigen]``
     - The real band-solve ``rhs`` holds the real and imaginary parts of one
       member of each conjugate pair, from which ``GC_collect`` is rebuilt
   * - ``GC_collect``
     - Caches :math:`G_{qe}c_e` for ordinary exponential evaluation; the
       conservative pair is evaluated from ``C_collect`` and its coupled
       two-column basis
     - Despite the name, stores only :math:`c_e`; multiplication by
       :math:`G_{qe}` occurs during field reconstruction
   * - Beam particular solution :math:`\boldsymbol B`
     - ``B_collect[m, layer, stream]``
     - ``B_collect[alpha, m, layer, stream]`` of Stokes vectors
   * - Polynomial particular solution
     - Coefficients are cached in
       ``source_collect[layer, stream, power]``; ``SRC0``, ``SRC1``, and
       ``SRCB`` retain boundary values used by the ordinary fast assembly
     - Coefficients are Stokes vectors in
       ``source_collect[alpha, m, layer, stream, power]``
   * - Layer-bottom total field
     - ``um[layer, m, stream]``
     - ``um[layer, alpha, m, stream]`` of Stokes vectors

In particular, ``um`` is a cached total quadrature field, not an eigenvector or
a modal-coefficient array.

Numerical behavior and limitations
**********************************

The following details are intentional and should be considered when changing
or comparing the cores:

* **Conservative scattering.**  Exact :math:`\omega=1` gives a defective zero
  pair in the zeroth cosine mode.  Both cores keep the input albedo unchanged
  and represent the physical energy pair by the constant/linear centered basis
  above.  The same basis is selected through the near-conservative interval
  :math:`[1-10^{-8},1]`; ordinary modes retain the fast anchored-exponential
  path.  VDISORT stabilizes the one physically required energy-conservation
  pair.  A contrived Mueller operator with additional independently conserved
  polarization quantities can contain more zero pairs and is not covered by
  this single-pair representation.  Generalizing this is intentionally deferred
  until a physically occurring scattering model requiring it is identified.
* **Complex VDISORT modes.**  A real physical solution can use complex
  eigenpairs.  The reconstructed imaginary part is required to cancel to a
  relative tolerance of about :math:`2\,10^{-8}` before the real part is
  returned.  Failure indicates an ill-conditioned or incorrectly assembled
  system, not a physical complex radiance.
* **VDISORT spectral split.**  Anchoring requires half the non-neutral modes to
  propagate toward each boundary.  After sorting, the implementation checks
  this sign split.  A relative neutral band of :math:`10^{-10}` times the
  eigenspectrum scale permits mathematically zero or nearly imaginary modes to
  occupy either half; a non-neutral sign imbalance or a mode on the wrong side
  is rejected before boundary assembly.

* **No absorption cutoff.**  An absorption-optical-depth shortcut is not
  implemented.  The complete physical solution is evaluated instead.
* **No special-boundary shortcut.**  ``IBCND=1`` albedo/transmission is not a
  separate algorithm here.  The same quantities are obtained from the regular
  boundary-value solution.


* **Update cost.**  Constructing or updating either solver performs the
  eigendecompositions and global boundary solves.  Evaluation at cached
  quadrature points is much cheaper.  VDISORT's update allocates nothing:
  its scratch (the transport operator, the eigensolution, the right-hand
  sides, the banded boundary system and the surface terms) is sized by the
  constructor and reused.  Only what is accumulated or written
  conditionally is reset: the transport operator and the band system before
  each mode, the conservative-pair flags, the right-hand side of the exact
  conservative pair, the source recursion, the surface terms of modes
  without a surface, and the beam and source terms of modes and layers
  without them.  For many frequencies, give each thread one solver from the
  size-only constructor, set its inputs and update it.
* **What an update saves.**  A layer and Fourier mode without scattering
  (single-scattering albedo or phase matrix zero) is not diagonalized
  numerically, since its eigenvectors are the states themselves.  A
  mirror-symmetric one is diagonalized through the half-size product problem
  (``mirror_symmetric`` and ``half_size_eigensolution``), except in the cosine
  system of order zero, and ``conjugate_partner`` records the conjugate of
  each complex mode.  The boundary system is real (``dgbsv`` instead of
  ``zgbsv``).  At 16 streams, 16 Fourier modes and 50 layers an update takes
  0.30 s instead of 0.99 s with every layer scattering, and 0.15 s instead
  of 0.43 s with 40 clear layers.  A Fourier
  mode that nothing drives (no source, beam, boundary radiance or surface
  term) has the zero solution and is neither assembled nor solved.  Each
  mode's exponential at a layer's top and bottom is computed once per
  assembly, not once per state.  The source polynomial of a layer shares
  one LU decomposition among all its degrees.  The gridded radiances and
  fluxes are then combinations of the cached layer-bottom fields.  The
  skipped work, the cached exponentials and the shared LU decomposition change
  no arithmetic (with the compiler's default contraction of multiplications
  and additions they move results by about one ulp).  The half-size problem
  and the real boundary system are exact reformulations that change results
  at the rounding level: by at most 1.3e-14 of the largest radiance in the
  test problems.  VDISORT intentionally contains no
  OpenMP parallel region in these core algorithms; ARTS normally parallelizes
  independent frequencies outside the solver.  The linear-algebra library may
  itself be threaded.  With an OpenMP-threaded BLAS and no outer parallel
  region, e.g. a single solve from Python or a C++ test, every small BLAS call
  inside the eigensolutions wakes and parks the library's threads.  In the
  polarized delta-M test this made the solves up to three times slower than
  with one thread.  Inside an outer parallel region OpenMP does not nest by
  default, so the library calls there run on one thread.
* **Thread use.**  A fully constructed solver may be shared for read-only
  evaluation when each caller owns its scratch objects.  Updating the solver
  or sharing scratch storage concurrently is unsafe.

Validation scope
****************

The scalar core is checked against published numerical reference values for
the supported plane-parallel test problems, including arbitrary angles, fluxes,
DFDT, delta-M corrections, physical BRDFs, conservative scattering, multilayer
continuity, and delta-M-plus.  The pseudo-spherical test case is intentionally
unsupported.

The strict scalar embedding of VDISORT is checked against the same supported
reference problems, with :math:`Q=U=V=0` asserted.  These tests establish the
normalization and scalar limit of the vector equations.  Analytic polarized
two-stream tests additionally cover :math:`I/Q` coupling, complex :math:`U/V`
eigenpairs, a polarized direct beam, polarized absorption, vector internal
sources, and a polarized reflecting boundary.  The comparisons with Evans'
doubling-adding solvers RT4 and RT3 (:doc:`dev.rt4`, :doc:`dev.rt3`) cover
reference-plane rotations, many-stream Mueller problems, all Fourier modes,
U and V, the direct beam, and the formal solution at user angles between the
streams.  Polarized delta-M with the IMS and TMS corrections is checked
against VDISORT itself without delta-M, on enough streams to carry the whole
phase-matrix series.  That reference shares the transport solver, which the
RT3 and RT4 comparisons validate, but not the scaling, the removed peak or the
corrections.  A polarized beam, and the reflection of a beam by the polarized
Cox-Munk surface, are checked against closed forms built from electric fields.
The coverage document in the VDISORT test directory lists the measured
agreement of every comparison.

VDISORT inputs from ARTS data
*****************************

``src/core/disort-cpp/vdisort_arts.h`` builds VDISORT inputs from ARTS
scattering species, atmospheric points and propagation paths.  Its Python
counterparts are in ``pyarts3.arts.vdisort``.  The header and the docstrings
are the reference for names, fields and defaults.  The builders follow the
same paradigms as those of RT3 and RT4 (:doc:`dev.rt3`, :doc:`dev.rt4`).

**The phase-matrix Fourier coefficients at an atmospheric point** come from
the azimuthal Fourier modes that the scattering species give of ARTS's
laboratory-frame phase matrix:

* They are taken at the propagation zenith angles of signed direction
  cosines (> 0 upward): VDISORT's streams, and the beam's cosine for the
  beam column.
* They are VDISORT's ordinary coefficients, without :math:`2-\delta_{m0}`,
  of the phase matrix ``P = 4 pi Z / sigma`` normalised to 1 over 4 pi.
* The phase integral ``sigma = int Z11 dOmega`` comes from the species'
  modes too, so nothing is integrated numerically.
* They are as accurate as the species' modes:

  * exact for the closed-form scatterers (Rayleigh gas scattering,
    Henyey-Greenstein) and for Legendre series;
  * gridded ARO data give them only at their own zenith angles;
  * gridded TRO data must first be converted to a Legendre series.
* ARTS's propagation directions (za, aa), with aa clockwise seen from above,
  are VDISORT's with ``mu = cos(za)`` and ``phi = -aa``.  ARTS's Stokes vector
  is VDISORT's meridional one (``h = k x z / |k x z|``, ``v = h x k``,
  Q = I_v - I_h, U = 2 Re(E_v E_h*)).  See :doc:`dev.rt3` for the test of
  this and for the sign of F34.
* VDISORT's beam propagates toward its azimuth ``phi0``, and its radiance at
  ``phi0 + psi`` is RT3's at ``psi``.
* VDISORT's optical depth and single-scattering albedo are scalars, so the
  particles' optics must be scalar:

  * K11, a1 and sigma must not depend on the incidence direction;
  * K12, K34 and a2 must vanish;
  * sigma must equal the scattering coefficient K11 - a1.

  These are checked to a tolerance relative to K11, and an infinite
  tolerance disables the checks.  Totally randomly oriented (TRO) particles
  always meet the first two conditions.  Azimuthally randomly oriented
  (ARO) particles generally do not, but ARO species whose optics are scalar
  can be used.

**A solved problem from a propagation path** uses the path conventions of
the RT3 and RT4 builders and of the DISORT workspace methods:

* The ray path, the atmospheric path and the gas propagation matrices have
  one entry per level, top first, with strictly decreasing altitudes in
  metres.
* The gas propagation matrices are unpolarized, per metre, and a layer's gas
  extinction is the mean of A at its two levels.
* The frequency is selected from the frequency grid by index.
* A layer has the mean extinction and scattering of its two levels'
  coefficients on VDISORT's streams, and the scattering-weighted mean of
  their normalised coefficients.
* The optical depth is the cumulative gas plus particle optical depth, and
  it must increase in every layer.  The single-scattering albedo is the
  scattering over the total extinction.
* With thermal emission, the source is ARTS's Planck function at the level
  temperatures, linear in optical depth within each layer, and the ground
  emits.  The sky is a blackbody at the sky temperature.
* A Lambertian ground uses VDISORT's Lambertian modes and emits
  ``[(1 - A) B, 0, 0, 0]``.
* A Fresnel ground uses VDISORT's Fresnel modes, whose reflection is their
  specular part R(mu), and emits ``B ([1, 0, 0, 0] - R[:, 0])``.  At an
  upward user angle, VDISORT reflects the downward user-angle radiance at
  the same cosine.  User angles must therefore come in up and down pairs.
* A beam has the Stokes irradiance ``[F / mu0, 0, 0, 0]`` normal to it, with
  F its flux on the horizontal plane.  Its cosine must not be a stream.

**Tests:**

* VDISORT's ARTS-input test in ``src/core/disort-cpp/test/vdisort`` checks:

  * the coefficients of ARTS's Rayleigh gas scattering against closed forms
    derived from the dipole Jones matrix in the meridional basis, for the
    diffuse and beam columns, including mu = +-1;
  * ARTS's laboratory-frame phase matrix of a polarizing Mie particle
    against the phase matrix of the polarization geometry
    (:ref:`Sec Polarization scattering`), with ``mu = cos(za)`` and
    ``phi = -aa``.  Both are built from that geometry, which its own test
    checks against electric fields, so this pins the direction convention
    and the sign of F34;
  * the coefficients of that particle against the Fourier modes of that
    phase matrix;
  * the coefficients of Henyey-Greenstein scattering against the addition
    theorem;
  * the modes of gridded ARO data (Rayleigh scattering on a grid) against
    the exact ones, to the interpolation error of the grid;
  * the path builder against its inputs.
* The same directory runs VDISORT, RT3 and RT4 on one ARTS atmosphere
  (see :doc:`dev.rt3` and :doc:`dev.rt4`).
* The Kirchhoff test there checks both solvers in an isothermal cavity over
  depolarizing and polarizing surfaces: the radiance must be the Planck
  radiance, unpolarized, on every stream and at user angles, to rounding.
  Giving the surface's emission as the upward boundary instead (the earlier
  convention) must miss it.

The VDISORT workspace layer
***************************

VDISORT has DISORT's workspace layer, method for method.

* What the two solvers share has one solver-neutral name: the quadrature,
  Legendre and Fourier dimensions, the quadrature and the profiles of the
  flux components with their group.  The named setup options are shared
  too, and a solver refuses an option that it does not support.
* The settings, the radiance field, the agendas and the methods are each
  solver's own, with ``vdisort`` in place of ``disort`` in their names and
  ``VDisort`` in place of ``Disort`` in the groups.
* In ``src/workspace_methods.cpp`` and ``src/workspace_meta_methods.cpp``,
  a lambda makes the VDISORT record from a DISORT record, and each DISORT
  record that VDISORT shares is followed by that call.  A new DISORT method
  gets its VDISORT counterpart by adding the call after it.
* VDISORT also has Fresnel surface methods, which DISORT lacks.

The settings and radiance types are in
``src/core/disort-cpp/vdisort_settings.h``.  Their conventions:

* Where DISORT has a number, VDISORT has a Stokes vector: the sun, the
  boundary conditions, the source polynomial and the radiance field.  The
  fluxes are those of I.
* The settings hold each layer's scattering as the Legendre moments of the
  six scattering-plane elements F11, F12, F22, F33, F34 and F44, normalised
  so that the first moment of F11 is one.  The Fourier coefficients of the
  phase matrices on the streams, and from the beam, are built from them per
  frequency when the solver is set up.
* The boundary conditions are in VDISORT's combined Fourier representation.
  An azimuth-independent Stokes vector [I, Q, U, V] is [I, Q, 0, 0] in the
  cosine system and [0, 0, U, V] in the sine system of mode 0.
* The surface is held per frequency as data, as in DISORT: the Fourier
  reflection between the streams, the reflection of the sun's beam with the
  cosine it was evaluated at, and, for VDISORT, a specular part (Fresnel)
  or none.  The setter methods evaluate a surface model once for constant
  parameters, again where the sun's zenith angle changes, and per frequency
  for spectral ones; the solver's surface is allocated once, from the number
  of modes and whether there is a specular part, and the settings copy each
  frequency into it in place.  A solver with a beam refuses a surface
  evaluated for another beam cosine.
* The sources follow DISORT's methods.  The layers' source is linear in
  optical depth, from the Planck function or the non-LTE source vector.
  Surface emission by temperature sets the upward boundary to the Planck
  radiance B(T), and the solver transmits it with 1 - R, Kirchhoff's
  emission of whatever surface the reflection methods set.  The
  sky is the cosmic background, nothing, or the Stokes vector that the
  observer agenda gives at the top of the path.
* The extinction is scalar, as in DISORT, so a polarized gas propagation
  matrix is refused.
* The radiance field's zenith and azimuth grids are those of the propagation
  direction, as is the line of sight of a ray path point.  ARTS's azimuth,
  clockwise seen from above, is VDISORT's ``phi = -aa``.  The settings' solar
  zenith and azimuth are those of the line of sight toward the sun, so the
  beam propagates at the solar azimuth plus 180 degrees.
* DISORT methods without a VDISORT core counterpart are not mirrored:
  delta-M-plus scaling, the subsurface and coupled calculations, and cDISORT.

**Tests:** ``tests/core/disort/vdisort-workspace.rt3.py`` checks the layer
against references derived outside VDISORT:

* the exact solution of the transfer equation through clear-sky layers whose
  Planck function is linear in optical depth, below a polarized sky and with
  an empty layer;
* the first-order closed form of single scattering of the sun by an
  optically thin Henyey-Greenstein layer, which pins the field's grids to
  propagation directions and the beam to the sun's position;
* RT3 on the same path, through its path builder, for thermal emission and
  for the sun, with polarizing Rayleigh and Mie scattering in all four Stokes
  components.

The sign of U and V relative to ARTS's sensor polarization convention is not
tested.

Caller conventions are described in :doc:`user.disort`.
