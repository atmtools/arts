"""
jac_targetsFinalize must reject duplicate error targets and line targets
that are not in abs_bands.
"""

import pyarts3 as pyarts
import numpy as np

f = np.linspace(-5e6, 5e6, 11) + 118750348044.712

ws = pyarts.workspace.Workspace()

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)

ws.subsurf_field = pyarts.arts.SubsurfaceField()
ws.surf_fieldPlanet(option="Earth")
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.measurement_sensorSimpleGaussian(
    freq_grid=f, std=1, pos=[100e3, 0, 0], los=[180.0, 0.0], pol="I"
)

band = list(ws.abs_bands.keys())[0]
nlines = len(ws.abs_bands[band].lines)


def assert_finalize_raises(msg):
    try:
        ws.jac_targetsFinalize()
    except Exception:
        return
    raise AssertionError(f"{msg}: jac_targetsFinalize did not raise")


# %% Duplicate error targets

ws.jac_targetsInit()
ws.jac_targetsAddErrorPolyFit(t=f - f[0], sensor_elem=0, polyorder=1)
ws.jac_targetsFinalize()

ws.jac_targetsAddErrorPolyFit(t=f - f[0], sensor_elem=0, polyorder=1)
assert_finalize_raises("Duplicate error target")

# %% Line targets must exist in abs_bands

ws.jac_targetsInit()
ws.jac_targetsAddLineParameter(band=band, line=0, parameter="f0")
ws.jac_targetsAddLineShapeParameter(
    band=band, line=0, species="O2", parameter="G0", coefficient="X0"
)
ws.jac_targetsFinalize()

ws.jac_targetsInit()
ws.jac_targetsAddLineParameter(band=band, line=nlines, parameter="f0")
assert_finalize_raises("Line index out of range")

ws.jac_targetsInit()
ws.jac_targetsAddLineShapeParameter(
    band=band, line=0, species="H2O", parameter="G0", coefficient="X0"
)
assert_finalize_raises("Line shape species not in the line")

ws.jac_targetsInit()
ws.jac_targetsAddLineParameter(band=band, line=0, parameter="f0")
ws.abs_bands = {}
assert_finalize_raises("Band not in abs_bands")
