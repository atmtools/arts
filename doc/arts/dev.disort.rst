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
     - Local :math:`4N_q\times4N_q` real ``A_real`` is copied to a complex
       matrix and diagonalized
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
     - The complex band-solve ``rhs`` is copied into ``GC_collect``
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
  quadrature points is much cheaper.  VDISORT intentionally contains no
  OpenMP parallel region in these core algorithms; ARTS normally parallelizes
  independent frequencies outside the solver.
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
U and V and the direct beam.  Genuinely polarized IMS/TMS corrections have no
external reference.

VDISORT inputs from ARTS data
*****************************

``src/core/disort-cpp/vdisort_arts.h`` (namespace ``vdisort``, part of
``vdisort-cpp``) builds VDISORT inputs from ARTS scattering species,
atmospheric points and propagation paths.  VDISORT has a scalar extinction,
so only species with totally randomly oriented (TRO) data can be used.

.. code-block:: cpp

  struct fourier_optics { Numeric extinction; Numeric scattering;
                          rtepack::muelmat_tensor3 cosine; rtepack::muelmat_tensor3 sine; };
  fourier_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                   const AtmPoint& atm_point, Numeric frequency,
                                   const Vector& mu_out, const Vector& mu_in, Index nfourier,
                                   Index azimuth_count, Index scattering_angle_count,
                                   Numeric normalisation_tolerance);
  struct lambertian_surface { Numeric albedo; };
  struct fresnel_surface { Complex refractive_index; };
  using surface = std::variant<lambertian_surface, fresnel_surface>;
  struct path_settings {
    Index nquad{16}; Index nfourier{1}; Index azimuth_count{64};
    Index scattering_angle_count{512}; Numeric normalisation_tolerance{1e-3};
    bool thermal{true}; Numeric beam_flux{0.0}; Numeric beam_mu{0.5}; Numeric beam_azimuth{0.0};
  };
  main_data main_data_from_path(const ArrayOfPropagationPathPoint& ray_path,
                                const ArrayOfAtmPoint& atm_path,
                                const ArrayOfPropmatVector& spectral_propmat_path,
                                const AscendingGrid& freq_grid, Index freq_index,
                                const ArrayOfScatteringSpecies& scattering_species,
                                const path_settings& settings, const surface& ground,
                                Numeric surface_temperature, Numeric sky_temperature);

The same names are in ``pyarts3.arts.vdisort`` (``FourierOptics``,
``scattering_optics``, ``LambertianSurface``, ``FresnelSurface``,
``PathSettings``, ``main_data_from_path``, which returns a
:class:`~pyarts3.arts.cppvdisort`).  Its read-only ``tau``, ``omega``,
``mu`` and ``weights`` give the layers and streams, ``u(tau, phi)`` the
radiance, and :func:`pyarts3.plots.cppvdisort.plot` draws it on the streams
at a layer boundary.

**scattering_optics** returns, for signed direction cosines ``mu_out`` and
``mu_in`` (> 0 upward; VDISORT's streams, and ``[-mu0]`` for the beam
column), the ordinary Fourier coefficients without :math:`2-\delta_{m0}`,
``C^m, S^m = (1 / 2 pi) int P(mu_o, 0; mu_i, phi) {cos, sin}(m phi) dphi``,
``m = 0 .. nfourier - 1``, of the laboratory-frame phase matrix
``P = 4 pi Z / sigma``:

* Z is ARTS's laboratory-frame phase matrix
  (``get_bulk_scattering_properties_aro_gridded``) from the incident
  ``za = acos(mu_i)`` to the scattered ``za = acos(mu_o)`` at
  ``delta_aa = aa_scat - aa_inc = phi``.  ARTS's propagation directions
  (za, aa), aa clockwise from above, are VDISORT's with ``mu = cos(za)`` and
  ``phi = -aa``, and its Stokes vector is VDISORT's meridional one
  (``h = k x z / |k x z|``, ``v = h x k``, Q = I_v - I_h,
  U = 2 Re(E_v E_h*); z up, ``k = (sin cos phi, sin sin phi, mu)``); see
  :doc:`dev.rt3` for the test and the F34 sign.
* The phi integral is the periodic midpoint rule at
  ``(k + 1/2) 2 pi / azimuth_count``.  For a regular Legendre series of
  degree L, Z is a trigonometric polynomial of degree L in phi, and the rule
  is exact for ``azimuth_count > L + nfourier - 1``.
* ``sigma = 2 pi int F11 dcos(Theta)`` by a ``scattering_angle_count``-point
  Gauss-Legendre rule normalises P to 1 over 4 pi.  ``extinction`` is K11 and
  ``scattering`` is K11 - a1; sigma must equal K11 - a1 to
  ``normalisation_tolerance`` times K11.
* The results go to ``combine_phase_matrices`` (diffuse) and
  ``combine_beam_phase_matrices`` (beam column).  VDISORT's beam propagates
  toward its azimuth ``phi0``; its radiance at ``phi0 + psi`` is RT3's at
  ``psi``.

**main_data_from_path** uses the path conventions of the RT3 and RT4
builders and of the DISORT workspace methods (one entry per level, top
first, strictly decreasing altitudes in metres, unpolarized gas propagation
matrix per metre, gas extinction the mean of A at the two levels, frequency
``freq_grid[freq_index]``).  A layer has the mean extinction and scattering
of its two levels' ``scattering_optics`` on VDISORT's streams and the
scattering-weighted mean of their coefficients; ``tau`` is the cumulative
(gas + particle) optical depth, which must increase in every layer, and
``omega`` the scattering over the total extinction.  With ``thermal`` the
source is ARTS's ``planck`` at the level temperatures, linear in optical
depth within each layer (``[c0, c1]`` in the global optical depth), and the
surface emits; the sky is a blackbody at ``sky_temperature``.  A Lambertian
surface uses ``brdf::lambertian_fourier_modes`` and emits
``[(1 - A) B, 0, 0, 0]``; a Fresnel surface uses
``brdf::fresnel_fourier_modes`` (the specular part R(mu); at an upward user
angle VDISORT reflects the downward user-angle radiance at the same angle),
and emits ``B ([1, 0, 0, 0] - R[:, 0])``.  A beam has the Stokes irradiance
``[beam_flux / beam_mu, 0, 0, 0]`` normal to it, and ``beam_mu`` must not be
a stream.

Tests: ``cpp.fast.vdisort-arts-test``
(``src/core/disort-cpp/test/vdisort/vdisort-arts-test.cpp``) compares
``scattering_optics`` for ARTS's Rayleigh ``GasScatterer`` with closed forms
of C^m and S^m for m = 0 .. 3 derived from the dipole Jones matrix in the
meridional basis (diffuse and beam column, including mu = +-1; 6.1e-15),
ARTS's laboratory-frame phase matrix with the vector geometry for a
polarizing Mie particle (6.7e-15), the coefficients with those of the
vector-geometry phase matrix (1.3e-15), and the path builder with its
inputs.
``cpp.fast.vdisort-arts-comparison`` runs VDISORT, RT3 and RT4 on one ARTS
atmosphere (see :doc:`dev.rt3` and :doc:`dev.rt4`).

Caller conventions are described in :doc:`user.disort`.
