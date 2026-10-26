"""
Adding Jacobian targets must invalidate an earlier jac_targetsFinalize.
"""

import pyarts3 as pyarts

ws = pyarts.workspace.Workspace()

ws.abs_bands = pyarts.arts.AbsorptionBands()
ws.measurement_sensor = []
ws.subsurf_field = pyarts.arts.SubsurfaceField()
ws.surf_fieldPlanet(option="Earth")
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_field["wind_u"] = 1.0


def assert_not_finalized(msg):
    try:
        n = ws.jac_targets.x_size()
    except Exception:
        return
    raise AssertionError(f"{msg}: x_size() returned {n} instead of raising")


ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
assert_not_finalized("Added without finalizing")

ws.jac_targetsFinalize()
n = ws.jac_targets.x_size()
assert n > 0

ws.jac_targetsAddPressure()
assert_not_finalized("Added after finalizing")

ws.jac_targetsFinalize()
assert ws.jac_targets.x_size() == 2 * n

ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
assert_not_finalized("Added after re-initializing")

ws.jac_targetsInit()
ws.jac_targetsAddWindField(component="u")
ws.jac_targetsFinalize()
ws.jac_targetsAddOverlappingWindField()
assert_not_finalized("Added overlapping wind after finalizing")

ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
ws.jac_targetsFinalize()
ws.jac_targets.atm.append(pyarts.arts.JacobianAtmTarget())
assert_not_finalized("Appended from Python after finalizing")
