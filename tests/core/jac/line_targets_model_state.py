"""
Line parameter and line shape targets must map between abs_bands and the
model state vector.
"""

import pyarts3 as pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)

ws.measurement_sensor = []
ws.subsurf_field = pyarts.arts.SubsurfaceField()
ws.surf_fieldPlanet(option="Earth")
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

band = list(ws.abs_bands.keys())[0]

ws.jac_targetsInit()
ws.jac_targetsAddLineParameter(band=band, line=0, parameter="f0")
ws.jac_targetsAddLineShapeParameter(
    band=band, line=0, species="O2", parameter="G0", coefficient="X0"
)
ws.jac_targetsFinalize()

ws.model_state_vecInit()
ws.model_state_vecFromBands()
x = np.array(ws.model_state_vec)

assert np.isclose(x[0], ws.abs_bands[band].lines[0].f0)
assert x[1] != 0.0

ws.model_state_vec = x * 1.5
ws.abs_bandsFromModelState(model_state_targets=ws.jac_targets)

assert np.isclose(ws.abs_bands[band].lines[0].f0, 1.5 * x[0])

ws.model_state_vecInit()
ws.model_state_vecFromBands()
assert np.allclose(ws.model_state_vec, 1.5 * x)
