.. _dev-arts2missing:

ARTS 2 Capabilities Missing from ARTS 3
=======================================

This page lists scientific and numerical capabilities of ARTS 2 that ARTS 3
does not have.  Anything that exists in ARTS 3 belongs on neither this page
nor :doc:`dev.arts2notintended`, however it was ported: as a workspace method,
a core library, Python bindings, a recipe, or a composition of existing
methods.  ARTS 2 interfaces that are deliberately not ported are listed in
:doc:`dev.arts2notintended`.

Rules for editing this page:

* Delete an entry in the same change that makes the capability available in
  ARTS 3.  Do not mark it as ported, reconnected, or partly done; the git
  history records that.
* Keep only the missing part of a capability that ARTS 3 partly provides.
* Track capabilities, not ARTS 2 spellings.  An ARTS 2 helper is not missing
  because ARTS 3 has no workspace method of that name if plain Python, a
  recipe, or a composition of ARTS 3 methods produces the same result.
* Missing test coverage of something ARTS 3 implements is not a missing ARTS 2
  capability and does not belong here.
* Each entry states the ARTS 2 capability and its evidence, the ARTS 3 design
  constraint, and the remaining action.  If a decision is made not to port an
  item, move it to :doc:`dev.arts2notintended` with the reason.

The list comes from comparing the ARTS 2.6 and ARTS 3 test inventories, with a
source audit where a test alone could not settle the question.  The complete
method and group inventories, optional builds, data formats, and ARTS 2
workflows without tests have not been compared, so a capability that is absent
from this page is not necessarily present in ARTS 3.

DOIT scattering solver
----------------------

ARTS 2 tests DOIT in ``TestDOIT``, ``TestDOITaccelerated``,
``TestDOITprecalcInit``, ``TestDOITsensorInsideCloudbox``, and
``TestDOITpressureoptimization``.  ARTS 3 has no DOIT implementation, and its
DISORT and Monte Carlo solvers are not a deterministic polarized solver of the
same kind.  A port must use the ARTS 3 atmospheric-field and
scattering-species models.  It must not restore the ARTS 2 cloudbox, global
spatial grids, or a global ``stokes_dim`` to preserve the old setup sequence.

FASTEM surface model
--------------------

ARTS 3 has no FASTEM.  FASTEM provides separate emissivity and reflectivity
and an atmospheric-transmittance correction, which the current surface
interface cannot represent: its closed-surface reflectance path applies
Kirchhoff consistency and has no atmospheric dependency.  A port needs a
surface-boundary contract that carries emission, reflection, and the
atmospheric coupling explicitly, without discarding FASTEM outputs or
weakening the consistency of the existing reflectance agenda.  The
``ENABLE_FASTEM`` hooks in ``config.h.cmake`` and ``src/CMakeLists.txt`` have
no implementation behind them.

MetMM instrument presets
------------------------

The ARTS 2 MetMM framework defines ATMS, DEIMOS, HATPRO, ISMAR up and down,
MARSS, MHS, MWHS-2, and SAPHIR.  The ARTS 3 predefined-sensor enumeration has
none of them.  Decide which channel definitions to port, and port them
individually as predefined sensors or ``pyarts3.arts.sensor.Builder`` setups
with an end-to-end reference case each.  The MetMM framework itself is not
intended for ARTS 3.

Raw calibration and time-series corrections
-------------------------------------------

ARTS 2 registers ``raw/calib.py`` and ``raw/corr.py`` for hot/cold
calibration, timestamp sorting, a simple tropospheric correction, and time
averaging.  ARTS 3 has no equivalent.  These are array and time-series
operations: if users need them, provide maintained Python recipes with the
ARTS 2 semantics rather than workspace methods.

NetCDF optical-depth export
---------------------------

ARTS 2 ``WriteMolTau`` exports molecular optical depths to NetCDF for
libRadtran.  ARTS 3 has no equivalent.  If it is still needed, a Python/xarray
exporter is the natural implementation.
