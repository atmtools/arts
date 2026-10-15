.. _dev-arts2notintended:

ARTS 2 Features Not Intended for ARTS 3
=======================================

This page lists ARTS 2 interfaces and implementation patterns that ARTS 3
deliberately does not port, with the reason for each decision.  The science
behind an item may still be needed: if ARTS 3 cannot produce the result, the
capability belongs in :doc:`dev.arts2missing` and must be implemented with the
ARTS 3 design.  Anything that exists in ARTS 3 belongs on neither page.

Rules for editing this page:

* List only ARTS 2 features that are absent from ARTS 3 by decision.  Name
  the ARTS 3 replacement in a sentence at most; do not document ARTS 3 here.
* If a listed feature is added to ARTS 3 after all, delete its entry in the
  same change.  Do not annotate it; the git history records that.

ARTS 2 scripting language and executable
----------------------------------------

The ARTS 2 controlfile language is not an ARTS 3 compatibility target.  This
includes its parser, the ``arts`` controlfile executable, syntax-only checks,
controlfile-to-Python conversion, and command-line introspection tied to that
runtime.  Python is the orchestration language, so the language's
conveniences are not ported either:

* ``ForLoop``: use Python loops.
* ``Copy``, ``AgendaCopy``, and ``ArrayOfAgenda``: use Python assignment,
  copies, lists, and callables where no native ARTS object semantics are
  required.
* ``ybatch``, ``DOBatch``, and ``YCalcAppend``: use Python loops and array
  assembly.  If a large workload needs a parallel C++ kernel, that kernel is a
  missing performance capability, not a reason to restore the construct.

Do not port these by adding Python implementations that pretend to be native
workspace methods.

Global model dimensionality and master grids
--------------------------------------------

ARTS 3 geometry is always three-dimensional on a reference ellipsoid, so the
following are not ported:

* ``AtmosphereSet1D`` and ``AtmosphereSet2D`` as global model modes;
* general or master ``z_grid``, ``lat_grid``, and ``lon_grid`` workspace
  variables that determine the shape of every field; and
* separate 1-D, 2-D, and 3-D propagation-path implementations.

Individual atmospheric quantities can still be constant or low-dimensional.

Global Stokes and Jacobian switches
-----------------------------------

The ``stokes_dim`` workspace variable is not ported.  ARTS 3 data types carry
the full Stokes vector, and sensors select polarization through Stokes
projection weights, so no global dimension switch is needed.

The global ``jacobian_do`` flag is not ported either.  Requested derivatives
are the Jacobian targets, so an empty target set, or a method's explicit
inputs, decides whether derivative work is done.

Formal cloudbox
---------------

``cloudbox_on``, ``cloudbox_limits``, the cloudbox checked flags, the cloudbox
field, and the cloudbox interpolation and setup machinery are not ported.
Scattering solvers take explicit atmospheric fields, scattering species,
paths or domains, and solver settings instead.  This does not remove the need
for the scattering solvers themselves; see DOIT in :doc:`dev.arts2missing`.

Grid-and-tensor field collections
---------------------------------

ARTS 2 describes the atmosphere, surface, and scattering state as global grids
plus parallel tensors.  ARTS 3 uses composite types whose dimensions and
invariants are fixed by the type.  Do not split ``AtmField``,
``SurfaceField``, ``SubsurfaceField``, the model state, or the scattering
state back into that collection to reproduce an ARTS 2 method signature.

Sensor workspace variables and instrument frameworks
----------------------------------------------------

The ARTS 2 global sensor position, line-of-sight, frequency-grid, response,
polarization, and normalization workspace variables are not ported; ARTS 3
describes a measurement with ``measurement_sensor``.  The
instrument-specific controlfile frameworks built on those variables, such as
MetMM with its dedicated grid and response methods, are not ported either.
Port verified channel definitions as predefined sensors or
``pyarts3.arts.sensor.Builder`` setups instead.

Scattering-data layout
----------------------

The ARTS 2 ``ScatSpeciesInit``, ``pnd_field``, and ``scat_data`` object layout
and its agenda setup are not ported; ARTS 3 has typed scattering species.
Port a PSD, habit, optical-property model, or solver result against the ARTS 3
scattering-species interface, and do not assume that a working ARTS 2 setup
sequence defines that interface.

Checked flags and silent switches
---------------------------------

The eager ``*_checkedCalc`` methods and their checked-state flags
(``atmfields_checkedCalc``, ``atmgeom_checkedCalc``, ``cloudbox_checkedCalc``,
``sensor_checkedCalc``, and others) are not ported.  ARTS 3 validates values
where they are used and checks the dimensions that several inputs share.

Silent atmosphere-class switches, class-ignore flags, and solver-side zeroing
are not ported.  A calculation that needs a quantity set to zero gets an
explicitly constructed field; a solver may ignore quantities outside its
contract, but must not mutate or zero them.

Retired ARTS 2 tests
--------------------

These ARTS 2 tests exercise the interfaces above.  They are not ported; their
useful behavior is tested through the ARTS 3 interfaces:

* ``cmdline.*``, ``converted.*``, ``TestWSMCalls``, and the
  controlfile-conversion example: the controlfile language.
* ``TestForloop``, ``TestAgendaCopy``, ``TestArrayOfAgenda``,
  ``TestDOBatch``, and ``TestYCalcAppend``: Python orchestration.
* ``TestPpath1D`` and ``TestPpath2D``: the common 3-D path implementation.
* ``TestCloudboxAuto`` and ``Testcloudbox_fieldInterp2Azimuth``: the
  scattering solvers' explicit inputs.
* ``TestHSE`` and the pressure-grid regridding family
  (``GriddedFieldPRegrid``, ``AtmFieldPRegrid``, master-grid refinement):
  ``atm_fieldHydrostaticPressure`` and independent field interpolation.
* ``jacobianAdjustAndTransform``: transformations of individual Jacobian
  targets and model-state mappings.
* the MetMM instrument controlfiles: predefined sensors and observation
  elements.

Deciding where an ARTS 2 feature belongs
----------------------------------------

Before porting an ARTS 2 test or method, ask:

#. Is a physical, numerical, data, or performance capability actually absent
   from ARTS 3?
#. Can ordinary Python or a maintained recipe express it clearly?
#. Which ARTS 3 field, sensor, scattering, path, or Jacobian type owns it?
#. Would the port restore global grids, cloudbox state, ``stokes_dim``,
   checked flags, or controlfile-language conveniences?
#. Can a compact regression compare a meaningful result without depending on
   the ARTS 2 source tree?

If only the ARTS 2 spelling or orchestration is absent, do not port it.  If
the result itself cannot be produced, add it to :doc:`dev.arts2missing` with
the evidence and an ARTS 3-native acceptance test.
