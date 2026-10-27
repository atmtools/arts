"""
With nonnegative relative humidity, the field does not depend on negative
model state elements, so their Jacobian columns must be zero.
"""

import pyarts3 as pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

ws.abs_bands = pyarts.arts.AbsorptionBands()
ws.measurement_sensor = []
ws.subsurf_field = pyarts.arts.SubsurfaceField()
ws.surf_fieldPlanet(option="Earth")
ws.abs_speciesSet(species=["H2O"])
ws.water_equivalent_pressure_operatorMK05()
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

ws.jac_targetsInit()
ws.jac_targetsAddSpeciesVMR(species="H2O")
ws.jac_targetsFinalize()
ws.jac_targetsToggleRelativeHumidityAtmTarget(key="H2O", nonnegative=1)

ws.model_state_vecInit()
ws.model_state_vecFromAtmosphere()
x = np.array(ws.model_state_vec)
assert np.all(x > 0)

negative = np.arange(x.size) % 2 == 0
x[negative] = -0.1
ws.model_state_vec = x

ws.measurement_jac = np.ones((3, x.size))
ws.measurement_jacAtmosphereTransformation()
jac = np.array(ws.measurement_jac)

assert np.all(jac[:, negative] == 0.0), "Clamped elements have a Jacobian"
assert np.all(jac[:, ~negative] > 0.0), "Unclamped elements lost their Jacobian"
