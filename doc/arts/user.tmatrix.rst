T-matrix interface
==================

The Python interface is ``pyarts3.arts.tmatrix``.  ``available()`` reports
whether the solver was built (``ENABLE_TMATRIX``); the solver is double
precision.  ``fixed`` computes a particle and evaluates its amplitude and phase
matrices at one illumination/scattering geometry.  ``fixed_batch`` computes
one particle and evaluates a list of geometries, returning a list of fixed
results in the same order.  A geometry is a ``tmatrix.Geometry``, with the
fields ``theta_incident``, ``theta_scattered``, ``phi_incident``,
``phi_scattered``, ``alpha`` and ``beta``, all in degrees and 0 by default.
``random`` computes a
randomly oriented size distribution, defaulting to an effectively monodisperse
particle.  ``azimuthally_random`` computes a particle whose axis has a
distribution of tilts from the vertical and every azimuth alike, on ARTS's
grids of scattering data.  These functions calculate individual optical
results.  Use ``ParticleHabit.tmatrix`` and ``ParticleHabit.tmatrix_aro``
below to prepare a habit for scattering species.

Units and polarization
----------------------

Angles are in degrees and are used as given, also along the particle axis
and at the laboratory poles, where the meridional basis is the limit along the
given azimuth.  For fixed particles and the default and power-law
size distributions, use a common length unit for radius and wavelength.
Metres give SI outputs.  ``equivalent_radius`` is an
:class:`~pyarts3.arts.EquivalentRadius`: the radius is that of the sphere of
the same volume (``"volume"``, the default) or of the same surface area
(``"surface_area"``).
``shape`` is a :class:`~pyarts3.arts.ParticleShape`.  ``"spheroid"`` has
``aspect_ratio`` its horizontal/rotational axis ratio; ``"cylinder"`` its
diameter/length ratio.  ``"chebyshev"`` is a Chebyshev particle
:math:`r_0(1+\varepsilon\cos n\theta)` of ``degree`` :math:`n\ge1`, with
``aspect_ratio`` the deformation :math:`\varepsilon`, :math:`|\varepsilon|<1`.
``degree`` is 0 for the other shapes.  Odd degrees are not mirror symmetric
about the equatorial plane; even degrees are.  ``"generalized_chebyshev"`` is
the particle :math:`r_0(1+c_0+\sum_{n\ge1}c_n\cos n\theta)` of the
``coefficients`` :math:`c_0,\ldots,c_N`, such as of a raindrop's shape, with
:math:`1+c_0>\sum_{n\ge1}|c_n|`; its ``degree`` is 0 and ``aspect_ratio`` is
unused.  It is mirror symmetric without odd coefficients.  The other shapes
take no coefficients.

The fixed result contains the amplitude matrix as a :class:`~pyarts3.arts.Jones`
matrix between ARTS's meridional bases of the two rays (length), its Mueller
matrix as a ``Muelmat`` phase matrix (length squared), and orientation-averaged
extinction and scattering cross sections (length squared).  Both are in ARTS's
Stokes convention (:doc:`concept.polarization`).

Random results contain the scattering matrix as a
:class:`~pyarts3.arts.CompactPlanarMuelmatVector`, its six independent elements
F11, F12, F22, F33, F34 and F44 in the scattering-plane basis: ``legendre``
holds its exact Legendre series, and ``phase`` its values at ``angles`` equally
spaced scattering angles from 0 to 180 degrees, none unless asked for.  They are dimensionless, with F11
integrating to :math:`4\pi` over solid angle.  Multiply by
``scattering / (4*pi)`` to obtain differential scattering cross sections.  They
are in ARTS's Stokes convention: F34 is the negative of the original solver's,
whose :math:`U` has the other sign.  ``distribution`` is a
:class:`~pyarts3.arts.TMatrixSizeDistribution`.

``azimuthally_random`` takes the grids of ARTS's azimuthally random scattering
data: ``za_inc`` and ``za_scat``, the zenith angles of the incident and
scattered propagation directions, and ``delta_aa``, the scattered minus the
incident azimuth, all in degrees.  ``tilt_angles`` are the tilts of the
particle's axis from the vertical, in degrees, with ``tilt_weights``, which
need not be normalized; a single tilt with weight 1 is that tilt alone.  The
result holds the phase matrix on the three grids as a ``MuelmatTensor3``
between ARTS's meridional bases, and the extinction matrix (K11, K12, K34) and
absorption vector (I, Q) on ``za_inc`` as compact planar types, in length
squared.  The average over the azimuth and the absorption vector's integral
over the incident directions are exact for the truncated T-matrix.

Results remain valid after subsequent calls, and calls may run concurrently.
Random orientation uses OpenMP threads, unless called inside a parallel region,
with results independent of the number of threads.

Random orientation computes Clebsch-Gordan coefficients from the Wigner tables,
which must be initialized first, for example with *WignerInit*.  The tables
bound the particle size: random orientation fails, naming the tables it needs,
once its truncation would exceed them.  Neither orientation has other limits
on size, accuracy or quadrature.  It is an error when the requested
``accuracy`` cannot be resolved in double precision, when the T-matrix does not
converge before it is no longer finite or before an added quadrature node stops
reducing its change, when the area of a (generalized) Chebyshev particle does
not converge to the accuracy, and when the power-law radii have no root that
double precision resolves.  Random orientation of a generalized Chebyshev
particle has no reference in the original program; it is checked against the
uniform tilt distribution of ``azimuthally_random``.  A particle too large to allocate fails at once.

See :doc:`concept.tmatrix` for the physical definitions, normalization and the
equations solved.
Build configuration and implementation notes are in :doc:`dev.tmatrix`.

Using particles as scattering species
-------------------------------------

``ParticleHabit.tmatrix`` generates a totally randomly oriented habit in memory.
Supply temperature and frequency grids, volume-equivalent diameters in metres,
a complex refractive-index matrix with axes (temperature, frequency), and
material density in kg/m³.  The refractive index uses a nonnegative imaginary
part.  Shape and aspect-ratio conventions match the direct T-matrix interface.
``aspect_ratio`` has one element for all diameters or one per diameter, and so
do the ``coefficients`` of a generalized Chebyshev particle, an
:class:`~pyarts3.arts.ArrayOfVector` that the other shapes leave empty.
Each particle's maximum dimension is the largest distance between two points of
its surface: closed forms for spheroids and cylinders, and a numerical maximum
for (generalized) Chebyshev particles.

For example, a single-size ice population can be set up as follows::

    import numpy as np
    import pyarts3 as pa

    A = pa.arts
    ws = pa.Workspace()
    ws.WignerInit()
    density = 917.0
    habit = A.ParticleHabit.tmatrix(
        t_grid=[250.0, 280.0],
        f_grid=[229e9, 231e9],
        diameters=[200e-6],
        refractive_index=np.full((2, 2), 1.78 + 0.01j),
        density=density,
        aspect_ratio=[1.5],
    )
    number = A.ScatteringSpeciesProperty(
        "ice", A.ParticulateProperty.NumberDensity)
    psd = A.MonodispersePSD(number, 250.0, 280.0)
    ws.scat_species = [A.ScatteringHabit(
        habit, psd, density * np.pi / 6, 3.0)]

The constant refractive index above is illustrative; supply the material model
appropriate to the calculation.  Set the ``number`` property in the atmosphere
to the particle number density in m⁻³.  The explicit mass-size relation uses
volume-equivalent diameter: mass equals ``density * pi / 6 * diameter**3``.
``MonodispersePSD`` requires exactly one diameter; habits containing several
sizes can be combined with the other ARTS PSDs.

The factory computes each particle size separately.  ``ScatteringHabit`` then
applies the atmospheric PSD and interpolates the stored optical properties;
it does not rerun T-matrix during bulk-property evaluation.  Choose the
temperature, frequency and diameter grids to resolve changes in the optical
properties.  The particles hold spectral TRO data, the exact Legendre series of
their scattering matrices, so there is no scattering-angle grid to choose:
``to_tro_gridded`` and gridded bulk properties evaluate the series at the
angles asked for.

The generated data already have physical cross-section units, including
forward and backscatter matrices.  No additional phase normalization is needed
before passing the habit to ``ScatteringHabit``.

``ParticleHabit.tmatrix_aro`` generates an azimuthally random habit, with the
arguments of ``ParticleHabit.tmatrix`` and those of ``azimuthally_random``:
``za_inc_grid``, ``delta_aa_grid``, a ``ZenithAngleGrid`` ``za_scat_grid``,
``tilt_angles`` and ``tilt_weights``, one tilt distribution for all diameters.
The angle grids must be strictly increasing.  Its particles hold gridded ARO
data, exact at the grid nodes.  Bulk properties
interpolate between them, so the grids must resolve the phase matrix.  The
forward and backscatter matrices are taken from the grids, exactly when
``delta_aa_grid`` holds 0 and 180 and ``za_scat_grid`` holds every incident
zenith angle and 180 degrees minus it, and interpolated otherwise.  The
factories do not generate particles of a fixed orientation.
