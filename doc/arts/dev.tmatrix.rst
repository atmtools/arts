T-matrix implementation and validation
======================================

The T-matrix solver computes scattering by homogeneous spheroids, finite
circular cylinders, and Chebyshev and generalized Chebyshev particles.  It is a C++ port of the
double-precision codes of Mishchenko (``3rdparty/tmatrix``), whose licence
permits use in not-for-profit research.  ``-DENABLE_TMATRIX=ON`` builds it and
registers that licence; LGPL builds refuse it.  No Fortran compiler is needed.
The port has no shared state, and its calls may run concurrently.

Structure
---------

``tmatrix.h`` is the interface: it validates the inputs, converts the radius
definition and assembles the results.  ``tmatrix_solver.h`` holds the solver,
one function per section of :ref:`Sec T-matrix equations`, which states every
equation solved.  The code names the section it implements, and keeps the
symbols of the equations.  A shape is the ``ParticleShape`` option with, for
Chebyshev particles, a degree, and for generalized Chebyshev particles the
user's coefficients, which the solver holds normalized to the volume.

The amplitude matrix is assembled from three parts: the particle frame of an
orientation, the angular functions of a direction, and the contraction of the
T-matrix with the angular functions of one fixed direction.  A contraction
costs as much as one amplitude matrix, and every further direction paired
with that fixed one costs one order of :math:`n_{\max}` less.  Azimuthally
random orientation contracts the incident direction of each incident zenith
angle and azimuth node for its phase matrices, and the outgoing direction for
the incident directions of its absorption integral.  It runs its incident
zenith angles in OpenMP threads, each computed by one thread, so the results
do not depend on the number of threads.  The method is ARTS's own, the
averages of fixed orientations by the rules of :ref:`Sec T-matrix aro`, and
not a port.

The port follows the original's numerical procedure: its convergence of the
truncation and of the quadrature, and its starting rules for the truncation
and the Bessel recurrences, with the cube root where the original writes
0.333333 and 0.33333.  A start one degree apart only shifts where the
convergence begins.  The double- and extended-precision originals differ only in such
precision-tuned constants, in which stages run in extended precision, and in an
equivalent form of one parity factor; the port is double precision.

The user chooses the accuracy and the quadrature density; the port puts no
limits of its own on them, on the truncation or on the quadrature.  Where a
method has an inherent limit, the port either uses the limit the mathematics
gives, at the point where the general expression stops being defined, or fails
where it cannot estimate one: the convergence fails on the problem itself, a
surface too large to allocate fails first on its largest arrays, and a root or
an integral that double precision cannot resolve fails with it.

It changes six things:

* The generalized Chebyshev particle takes the user's coefficients, where the
  original fixes those of one drop in its code.  Its volume is exact by a
  Gauss-Legendre rule of the polynomial's degree and its area is integrated to
  the requested accuracy, where the original uses a 60-point rule for both.
* Every array is double precision and sized to the truncation; the original
  stores the T-matrix and the intermediate arrays of the orientation average in
  single precision, in arrays of a fixed maximal size.
* The expansion series of random orientation is complete; the original ends it
  at the first coefficient below a fixed threshold, whatever the requested
  accuracy.
* The convergence has no limits from array sizes; where the original reaches
  the limit of its quadrature unconverged and uses that matrix, the port fails
  once an added node no longer reduces the change.
* The amplitude matrix is evaluated at the given angles.  The original moves
  them off the singular directions, because its particle-frame azimuth divides
  by :math:`\sin\theta_p`, and uses the pole limits of the angular functions
  near the axis.  The port computes both particle-frame angles with ``atan2``
  from the components of the direction, which needs no division and is accurate
  to rounding next to the axis, and uses the limits only on the axis itself
  (:ref:`Sec T-matrix amplitude`).
* Shared ARTS routines replace the original's own: ``Legendre`` for the
  Gauss-Legendre quadratures, the Schmidt functions behind :math:`d^n_{0m}` and
  the Legendre polynomials, the LU solve of ``right_solve_inplace`` for
  :math:`T=-\mathrm{Rg}Q\,Q^{-1}` in place of the original's inverse and
  product, and wigxjpf for the Clebsch-Gordan coefficients.  Boost's TOMS 748
  root finder solves for the power-law radii on their whole interval, where the
  original uses Brent's ZEROIN in a bracket that excludes radii near its ends.
  Boost's adaptive Gauss-Kronrod rule integrates the area of a Chebyshev
  particle to the requested accuracy, where the original uses a fixed 60-point
  rule; a generalized Chebyshev particle's likewise.  The spherical Bessel recurrences serve the complex arguments
  inside the particle, which Boost does not support, and the real ones outside;
  they and the generalized spherical functions are ported.

The Clebsch-Gordan coefficients follow the line-mixing code: the Wigner tables
must have been initialized (*WignerInit*), which random orientation therefore
requires.  The expansion checks that the tables cover its symbols before it
starts, and allocates the per-thread wigxjpf scratch once around its loops.
wigxjpf evaluates each symbol exactly, at a cost that grows with the angular
momenta; the original's recursion is cheaper but loses accuracy as
:math:`n_1` grows.

The expansion runs its independent outer loops, over :math:`n` and over blocks
of :math:`l`, in OpenMP threads unless it is called inside a parallel region.
Each term is computed by one thread in the order of its sums, so the results do
not depend on the number of threads.

Storage is allocated around the loops that use it: a surface's special
functions write into its own arrays, the orders of a surface share one set of
work arrays and one :math:`Q`, :math:`T` is solved in its own storage, the
amplitude matrices of a batch of geometries share their angular functions, and
each thread of the expansion its sums.

Not ported: the extended precision; the original's
diagnostic van der Mee and Hovenier test, which only prints; the inversions it
does not select (Gaussian elimination, and Waterman's method through NAG).

Validation
----------

``src/core/tmatrix/test/tmatrix-sphere-test.cpp`` compares spheres with the
independent Mie solver of the scattering module: cross sections, the
random-orientation scattering matrix, and all sixteen elements of the
fixed-orientation phase matrix against the laboratory-frame phase matrix ARTS
builds from Mie's, also at the poles, in exact backscattering and at
:math:`\beta=90` degrees.  The spheres are spheroids and undeformed Chebyshev
particles of even and odd degree, so they also cover the integration over the
whole sphere.  The tolerances are the requested accuracy and, for the
scattering matrix, the end of the expansion series.

``tmatrix-singular-test.cpp`` approaches the particle axis, the laboratory
poles and :math:`\beta=90` degrees for spheroids, cylinders and Chebyshev
particles.  Along a smooth path the amplitude matrix is linear plus a second
order term, and the test bounds the residual of the linear term by the second
order, estimated from the path, plus the rounding of the evaluations.  A limit
taken away from the singular direction, or an angle that loses accuracy next to
it, fails it.

``tmatrix-shape-test.cpp`` compares the numerical maximum dimension with the
closed forms of spheroids and cylinders, and the radius ratio of Chebyshev
particles with fixed-order Gauss-Legendre rules, which give their volume
exactly and their area to rounding.  A generalized Chebyshev particle of one
coefficient is the Chebyshev particle of that degree: its surface, volume,
maximum dimension and radius ratio agree, and only odd coefficients break its
mirror symmetry.

``tmatrix-aro-test.cpp`` checks azimuthally random orientation.  Spheres of
several tilts against ARTS's conversion of Mie's scattering matrix to the
laboratory frame in all sixteen elements, and against Mie's cross sections.
A uniform tilt distribution, Gauss-Legendre nodes in :math:`\cos\beta`,
against the laboratory frame of random orientation's series, for a spheroid,
an odd Chebyshev particle and a generalized Chebyshev particle, whose random
orientation the original cannot compute.  The averages over the azimuth and
the absorption integral against ``fixed_batch`` with twice their nodes, to
rounding; fewer nodes than the rules give fail it.  The extinction against the
coherent transmission of the forward amplitude, the mirror symmetry in
:math:`\Delta\varphi`, and the vanishing absorption of a lossless particle.
``tests/core/tmatrix/aro.tmatrix.py`` checks the Python interface: a vertical
axis against ``fixed_batch`` at the angles of :ref:`Sec T-matrix aro`, the
normalization of the tilt weights, and the rules of the arguments.

``tmatrix-reference-test.cpp`` compares the port with outputs of the original
Fortran, which ARTS no longer builds.  ``tmatrix_fortran_cases.h`` holds the
inputs and outputs of the double-precision programs, recorded from them as C++
data, over spheroids and cylinders, prolate, spherical and oblate, Chebyshev
particles of even and odd degree and either sign of deformation, generalized
Chebyshev particles in fixed orientation with the original's drop coefficients
as test data, both radius definitions, absorbing and weakly absorbing particles,
polar and edge-of-range geometries, all size distributions and several
quadrature densities.  The Fortran moves the angles of the amplitude matrix off
the singular directions, so the port is given the moved angles.  The tolerances
follow the differences of the two: rounding times the truncation, for the
replaced quadrature nodes, times the componentwise condition number of
:math:`Q` for the cross sections, and the single-precision storage of the
Fortran times the magnitudes of the terms of the amplitude series and of the
expansion coefficients for the amplitude and scattering matrices.  The
comparison and its tolerances are in ``tmatrix_comparison.h``.

The test also compares the port with every number the original programs print
in ``3rdparty/tmatrix/*.ref``, of the double- and extended-precision programs
for fixed and random orientation, run with the inputs of their main routines:
amplitude and phase matrices, truncation, radius ratio, power-law radii,
effective radius and variance, cross sections, albedo, asymmetry, expansion
coefficients and scattering matrices.  The tolerance of a printed number is
half a unit of its last digit plus the bound of the Fortran comparison; the
extended-precision programs also store the T-matrix in single precision.  Each
expansion coefficient has its own bound: the relative single-precision error
that the bound of the scattering matrices sums over the coefficients, plus the
terms the original drops.

``tests/core/tmatrix/reference.tmatrix.py`` compares the Python interface with
the unchanged printed outputs of the original programs (``*.ref``) to their
printed precision, at the angles the original moves to, and checks unit
scaling, invalid inputs, the failures of the convergence and of the power-law
radii, ownership and concurrent calls.

See :doc:`user.tmatrix` for the interface and :doc:`concept.tmatrix` for
physical definitions and the equations.

Native particle-habit integration
---------------------------------

``ParticleHabit::tmatrix`` in ``particle_habit_tmatrix.cc`` generates native
TRO spectral single-scattering data, and ``ParticleHabit::tmatrix_aro`` native
ARO gridded data, without legacy-data adapters or files.
The scattering library links to the T-matrix interface; when the solver is
disabled the factory reports its unavailability through the same exception
path as the direct interface.

For each diameter, temperature, and frequency, the factory uses the direct
random-orientation interface with one size quadrature point at the specified
volume-equivalent radius.  It does not perform the atmospheric PSD integration.
The particle holds the solver's exact Legendre series of the scattering
matrix, multiplied by ``Csca/(4*pi)``, to the degree of the longest series of
its temperatures and frequencies.  Extinction is Cext and absorption is Cext
minus Csca.  Both endpoint scattering matrices are the series at 180 and 0
degrees.  The laboratory frame of the series, as azimuthally random bulk
properties take it, is the series at the exact scattering angle of every
direction pair.

The ARO factory calls the direct azimuthally random interface for each
diameter, temperature and frequency on its grids, and extracts the endpoint
matrices from them.  Both factories take the aspect ratio, and the
coefficients of generalized Chebyshev particles, as one entry for every
diameter or one entry each, in the ARTS types that XML reads and writes.

Particle mass is density times the equivalent-sphere volume.  Maximum
dimension is the largest point-to-point distance of the surface: the larger
axial diameter of a spheroid, the diagonal of a cylinder, and a numerical
maximum for a (generalized) Chebyshev particle (:doc:`concept.tmatrix`).
Temperature/frequency grids and all material inputs are validated before
starting the solver calls, and bulk evaluation subsequently uses the ordinary
habit machinery.

``tests/core/tmatrix/habit.tmatrix.py`` checks the direct-to-native normalization,
the series against the solver's own scattering matrix, forward/backscatter
extraction, metadata, and the analytic small-sphere Rayleigh limit, evaluating
the series independently with numpy.  It also exercises temperature/frequency
interpolation and repeated number-density scaling through
``ArrayOfScatteringSpecies``, verifying that these evaluations leave the
original particle data intact.  It checks per-diameter shapes against
habits of single diameters and the generalized Chebyshev maximum dimension
against a grid search, every ARO particle against the direct interface and its
endpoint matrices against its grids, and the ARO bulk properties of a sphere
against those of its TRO habit.

MC azimuthally random particle reproduction
-------------------------------------------

``tests/core/scat/mc_general_arts2.tmatrix.py`` generates its oblate ice particle
in memory, replacing the large scattering XML fixture.  This test requires
``ENABLE_TMATRIX``.  ``tools/compare_mc_tmatrix.py`` reuses the generator to
compare against an explicitly supplied original ARTS2 XML.  It uses
``tmatrix.fixed_batch`` to compute one T-matrix per frequency/temperature and
evaluate all geometries.  The single geometry interface delegates to the same
implementation.  Each returned result owns its amplitude and Mueller matrix.
The C++ batch interface takes a span of ``Geometry`` and writes into a
caller-owned ``std::span<FixedResult>`` with one element per geometry.  The
scalar wrapper uses stack storage.  A vector-returning C++ convenience overload
allocates the output and delegates to the span overload; the Python binding
uses this convenience overload.

The test records the four refractive indices obtained independently from the
historical PyARTS REFICE material model.  With vertical symmetry axis, an oblate
spheroid needs no orientation averaging.  Directional extinction comes from
forward amplitudes, and absorption is extinction minus the angular scattering
integral.  The original XML is read only after generation to report differences.
The generated data use the legacy ARO layout and existing native habit adapter.
Like the ARTS2 data they reproduce, they hold Mishchenko's Stokes convention:
the test converts the solver's phase matrices from ARTS's convention for them
(:doc:`concept.polarization`).

Run with a T-matrix-enabled Python environment::

    python tools/compare_mc_tmatrix.py --reference /path/to/original.xml

The default output and JSON error report go to ``tmp/mc-tmatrix-comparison``.
``--legacy-numerics`` tests historical float32 wavelength rounding and rounded
quadrature coefficients; ``--quadrature`` changes the absorption integration
order.  The original reference path is protected from overwrite.
