"""
Overlapping targets retrieve the magnitude of a vector field.  They must be
magnitude fields and cannot be transformed, and the shared model state must
map back to all components.
"""

import pyarts3 as pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

ws.abs_bands = pyarts.arts.AbsorptionBands()
ws.measurement_sensor = []
ws.subsurf_field = pyarts.arts.SubsurfaceField()
ws.surf_fieldPlanet(option="Earth")
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

alts = [0.0, 50e3, 100e3]
grids = (alts, [0.0], [0.0])


def set_wind(absolute):
    for key, scl in (("wind_u", 1.0), ("wind_v", -2.0), ("wind_w", 0.5)):
        data = pyarts.arts.Tensor3(np.array([1.0, 2.0, 3.0]).reshape(3, 1, 1) * scl)
        ws.atm_field[key] = pyarts.arts.GeodeticField3(data=data, grids=grids)
        ws.atm_field[key].set_extrapolation(extrapolation="Nearest")
    if absolute:
        ws.atm_fieldAbsoluteWindField()


def add_overlapping_wind():
    ws.jac_targetsInit()
    ws.jac_targetsAddWindField(component="u")
    ws.jac_targetsAddOverlappingWindField()


def assert_raises(f, msg):
    try:
        f()
    except Exception:
        return
    raise AssertionError(f"{msg}: did not raise")


def winds():
    return np.array(
        [[ws.atm_field[k](a, 0.0, 0.0) for a in alts] for k in ("wind_u", "wind_v", "wind_w")]
    )


# %% The shared model state maps back to all three components

set_wind(True)
add_overlapping_wind()
ws.jac_targetsFinalize()
assert ws.jac_targets.x_size() == len(alts)

orig = winds()
ws.model_state_vecInit()
ws.model_state_vecFromAtmosphere()
ws.model_state_vec *= 2.0
ws.atm_fieldFromModelState(model_state_targets=ws.jac_targets)
assert np.allclose(winds(), 2 * orig), "Overlapping components not scaled with the magnitude"

# %% Overlapping targets cannot be transformed, whatever the order of the calls

set_wind(True)
add_overlapping_wind()
ws.jac_targetsFinalize()
assert_raises(
    lambda: ws.jac_targetsToggleRelativeAtmTarget(key="wind_u"),
    "Transformed original after adding overlaps",
)
assert_raises(
    lambda: ws.jac_targetsToggleLogarithmicAtmTarget(key="wind_v"),
    "Transformed an overlapping target",
)

ws.jac_targetsInit()
ws.jac_targetsAddWindField(component="u")
ws.jac_targetsToggleRelativeAtmTarget(key="wind_u")
assert_raises(ws.jac_targetsAddOverlappingWindField, "Added overlaps to a transformed original")

# %% Overlapping targets must be magnitude fields

set_wind(False)
add_overlapping_wind()
assert_raises(ws.jac_targetsFinalize, "Finalized overlaps of plain component fields")
