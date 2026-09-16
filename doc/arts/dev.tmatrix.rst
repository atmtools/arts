T-matrix implementation and validation
======================================

The optional Mishchenko T-matrix backend computes scattering by homogeneous
spheroids and finite circular cylinders.  Enable it with
``-DENABLE_TMATRIX=ON``.  A Fortran compiler is required only when enabled.
``-DENABLE_TMATRIX_QUAD=ON`` selects the retained extended-precision solver;
inputs and outputs remain double precision, and both solvers retain some
single-precision internal storage.  The two backends are alternative builds.

On macOS with GNU Fortran, the backend is built as ``libtmatrix.dylib``
instead of a static archive.  The Fortran driver links this library with its
own Darwin unwind settings, avoiding compact-unwind conversion warnings for
GNU stack frames.  C++ consumers retain their normal compact-unwind tables
and exception handling.  Applying ``-no_compact_unwind`` to the mixed C++
link instead breaks exception handling with the tested Clang/ARM64 toolchain, even with
``-keep_dwarf_unwind``.  The shared backend must accompany the ARTS binaries;
CMake supplies the build-tree runtime search path.  Other platforms retain
the static backend.  No Fortran source, numerical compiler options, or
warning-suppression flags are changed.

Each call owns its output and holds one common mutex throughout computation
and amplitude evaluation.  Results survive subsequent calls.  Concurrent
calls are safe but do not run the Fortran solver concurrently.

Port and reference tests
------------------------

The complete ARTS2 ``3rdparty/tmatrix`` directory was retained, including its
original programs, parameter files, README, license and four reference files.
Array limits and numerical algorithms were preserved.  The build file was
adapted for optional Fortran compilation, position-independent code, compiler
symbol mangling, and GNU equivalents of the old quad-precision intrinsics.
The reference programs remain optional targets ``tmatrix_ampld`` and
``tmatrix_tmd``.

The following defects in the ARTS-specific sources were corrected:

* Initialize the Fortran error string on entry, and consistently use default
  Fortran integers for the quiet flag in callers and callees.
* Return immediately when the random-orientation stored expansion exceeds
  NPN4; the old wrapper reported failure but continued writing the arrays.
* Return distribution-averaged CEXT and CSCA, matching the returned albedo and
  phase matrix, rather than cross sections from the last size quadrature point.
* In the quad fixed wrapper, honor the integer shape argument, call its own
  ACJB routine, and stop redirecting standard output to a file named ``test``.
* Declare all six output arrays in the quad random wrapper.  Complete the
  ARTS LAPACK symbol prefix on its ``tmzswap`` implementation.

``tests/core/tmatrix/reference.tmatrix.py`` exercises the selected backend through
nanobind and the C++ interface, using the unchanged ``.ref`` files.  It compares
complex amplitudes, Mueller elements, size-distribution cross sections,
albedo, asymmetry and effective size statistics.  Printed precision determines
the comparison tolerances.  Timing and printed expansion-coefficient tables
are not part of the exposed result or these comparisons.

The fixed double reference uses volume-equivalent radius; the fixed quad
reference uses surface-area-equivalent radius.  The random double reference
uses refractive index 1.53+0.008i; the quad reference uses 1.33+0.001i.
These differences are inputs, not relaxed numerical comparisons.

The random reference's two distributions are evaluated separately, preserving
its seven and four size quadrature points.  The final extra CEXT/CSCA line in
``tmatrix_tmd.ref`` records the old return-value bug; the test compares the
unchanged distribution-average summaries printed earlier in that same file.
No reference numbers were replaced.  Radius normalization at the C++ boundary
for the power-law and gamma distributions avoids the original POWER routine's
absolute root-bracket scale when using SI-sized particles.

Tests also cover unit scaling, invalid inputs, convergence failure, ownership,
and interleaved Python-thread calls.  The original backend retains its
convergence envelope and some internal fatal error paths; the port is not a
claim that arbitrary particle parameters converge.  Input angle validation
prevents the directly exposed AMPL angle-error STOP paths.

The reference test reads each selected ``.ref`` file once.  The two original
random size distributions are each solved once for comparison; the smaller
result is reused as the baseline for the SI scaling check.  Ownership and
thread-safety checks use small absorbing nonspherical particles, with a barrier
to exercise concurrent entry.  They do not repeat the expensive reference
size distributions.  This preserves the original reference inputs and
tolerances while avoiding redundant serialized Fortran work.

See :doc:`user.tmatrix` for the interface and :doc:`concept.tmatrix` for
physical definitions and normalization.

Native particle-habit integration
---------------------------------

``ParticleHabit::tmatrix`` in ``particle_habit_tmatrix.cc`` generates native
TRO gridded single-scattering data, without legacy-data adapters or files.
The scattering library links to the optional-backend adapter; when the backend
is disabled the factory reports its unavailability through the same exception
path as the direct interface.

For each diameter, temperature, and frequency, the factory uses the direct
random-orientation interface with one size quadrature point at the specified
volume-equivalent radius.  It does not perform the atmospheric PSD integration.
The Mueller matrices are multiplied by ``Csca/(4*pi)`` and packed in native
TRO order F11, F12, F22, F33, F34, F44.  Extinction is Cext and absorption is
Cext minus Csca.  Both endpoint scattering matrices are extracted from the
populated phase data.

Particle mass is density times the equivalent-sphere volume.  For spheroids,
maximum dimension is the larger axial diameter.  For cylinders it is the
largest point-to-point distance, including both cylinder length and diameter.
Temperature/frequency grids and all material inputs are validated before
starting the solver calls.  The existing T-matrix mutex serializes Fortran
access, and bulk evaluation subsequently uses the ordinary habit machinery.

``tests/core/tmatrix/habit.tmatrix.py`` checks the direct-to-native normalization and
component ordering, forward/backscatter extraction, metadata, and the analytic
small-sphere Rayleigh limit.  It also exercises temperature/frequency
interpolation and repeated number-density scaling through
``ArrayOfScatteringSpecies``, verifying that these evaluations leave the
original particle data intact.

MC azimuthally random particle reproduction
-------------------------------------------

``tests/core/scat/mc_general_arts2.tmatrix.py`` generates its oblate ice particle
in memory, replacing the large scattering XML fixture.  This test requires
``ENABLE_TMATRIX``.  ``tools/compare_mc_tmatrix.py`` reuses the generator to
compare against an explicitly supplied original ARTS2 XML.  It uses
``tmatrix.fixed_batch`` to compute one T-matrix per frequency/temperature and
evaluate all geometries while holding the Fortran solver lock.  The single
geometry interface delegates to the same implementation.  Each returned
result owns its amplitude and Mueller matrix.  The C++ batch interface accepts
``ConstMatrixView`` input and writes into a caller-owned
``std::span<FixedResult>`` with one element per geometry row.  The scalar wrapper
uses stack storage.  A vector-returning C++ convenience overload allocates the
output and delegates to the span overload; the Python binding uses this
convenience overload.

The test records the four refractive indices obtained independently from the
historical PyARTS REFICE material model.  With vertical symmetry axis, an oblate
spheroid needs no orientation averaging.  Directional extinction comes from
forward amplitudes, and absorption is extinction minus the angular scattering
integral.  The original XML is read only after generation to report differences.
The generated data use the legacy ARO layout and existing native habit adapter.
MC acceptance values and native ``ParticleHabit.tmatrix`` orientation support
remain unchanged.

Run with a T-matrix-enabled Python environment::

    python tools/compare_mc_tmatrix.py --reference /path/to/original.xml

The default output and JSON error report go to ``tmp/mc-tmatrix-comparison``.
``--legacy-numerics`` tests historical float32 wavelength rounding and rounded
quadrature coefficients; ``--quadrature`` changes the absorption integration
order.  The original reference path is protected from overwrite.
