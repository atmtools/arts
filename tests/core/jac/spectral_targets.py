"""Spectral Jacobian targets: their model state, its round trip through the fields, XML, and OEM covariance blocks.

Each spectral target is one property of atm_spectral_field, surf_spectral_field
or subsurf_spectral_field.  Its model state is the flat data of the property,
frequency slowest, and the spectral targets follow the targets of their field
in the model state.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
ws = pyarts.Workspace()
ws.atm_field = arts.AtmField()
ws.surf_field = arts.SurfaceField()
ws.surf_field.ellipsoid = [1, 1]
ws.subsurf_field = arts.SubsurfaceField()
ws.abs_bands = arts.AbsorptionBands()
ws.measurement_sensors = arts.ArrayOfSensor()

freq = [100e9, 200e9]
ws.atm_field["t"] = arts.GriddedField3(
    data=np.array([250.0, 240.0, 230.0]).reshape(3, 1, 1),
    grid_names=["altitude", "latitude", "longitude"],
    grids=[[0.0, 1e3, 2e3], [0.0], [0.0]],
)
atm = np.arange(6.0).reshape(2, 3, 1, 1) + 1.0
surf = np.arange(4.0).reshape(2, 1, 2) + 10.0
subsurf = np.arange(4.0).reshape(2, 2, 1, 1) + 20.0
ws.atm_spectral_field["extinction"] = arts.SpectralGeodeticField4(
    data=atm, grid_names=["frequency", "altitude", "latitude", "longitude"], grids=[freq, [0.0, 1e3, 2e3], [0.0], [0.0]]
)
ws.surf_spectral_field["depolarization factor"] = arts.SpectralGeodeticField3(
    data=surf, grid_names=["frequency", "latitude", "longitude"], grids=[freq, [0.0], [-10.0, 10.0]]
)
ws.subsurf_spectral_field["loss tangent"] = arts.SpectralGeodeticField4(
    data=subsurf, grid_names=["frequency", "altitude", "latitude", "longitude"], grids=[freq, [-1.0, 0.0], [0.0], [0.0]]
)

# The targets in the order of the model state, though added in another order
ws.jac_targetsInit()
ws.jac_targetsAddSubsurfaceSpectral(target="loss tangent")
ws.jac_targetsAddSurfaceSpectral(target="depolarization factor")
ws.jac_targetsAddAtmosphereSpectral(target="extinction")
ws.jac_targetsAddTemperature()
ws.jac_targetsFinalize()

jt = ws.jac_targets
assert (jt.atm[0].x_start, jt.atm[0].x_size) == (0, 3)
assert (jt.atm_spectral[0].x_start, jt.atm_spectral[0].x_size) == (3, 6)
assert (jt.surf_spectral[0].x_start, jt.surf_spectral[0].x_size) == (9, 4)
assert (jt.subsurf_spectral[0].x_start, jt.subsurf_spectral[0].x_size) == (13, 4)
assert jt.x_size() == 17

# The model state is the flat data, frequency slowest
ws.model_state_vecInit()
ws.model_state_vecFromAtmosphere()
ws.model_state_vecFromAtmosphereSpectral()
ws.model_state_vecFromSurfaceSpectral()
ws.model_state_vecFromSubsurfaceSpectral()
x = np.asarray(ws.model_state_vec).copy()
assert np.array_equal(x, np.concatenate([[250.0, 240.0, 230.0], atm.ravel(), surf.ravel(), subsurf.ravel()]))

# The fields from a changed model state, and back
ws.model_state_targets = ws.jac_targets
ws.model_state_vec = 2.0 * x
ws.atm_spectral_fieldFromModelState()
ws.surf_spectral_fieldFromModelState()
ws.subsurf_spectral_fieldFromModelState()
assert np.array_equal(np.asarray(ws.atm_spectral_field["extinction"].data.data), 2.0 * atm)
assert np.array_equal(np.asarray(ws.surf_spectral_field["depolarization factor"].data.data), 2.0 * surf)
assert np.array_equal(np.asarray(ws.subsurf_spectral_field["loss tangent"].data.data), 2.0 * subsurf)

# A target whose property is not in its field cannot be finalized
ws.jac_targetsAddAtmosphereSpectral(target="ssa")
try:
    ws.jac_targetsFinalize()
    raise AssertionError("A spectral target must be a property of its field")
except RuntimeError as e:
    assert "ssa" in str(e), str(e)

# The targets survive XML
ws.jac_targetsInit()
ws.jac_targetsAddAtmosphereSpectral(target="extinction")
ws.jac_targetsAddSurfaceSpectral(target="depolarization factor")
ws.jac_targetsAddSubsurfaceSpectral(target="loss tangent")
ws.jac_targetsFinalize()
with tempfile.TemporaryDirectory() as directory:
    path = os.path.join(directory, "targets.xml")
    pyarts.xml.save(ws.jac_targets, path)
    loaded = pyarts.xml.load(path)
assert str(loaded) == str(ws.jac_targets)

# OEM puts the covariance block of each spectral target on its part of the diagonal
ws.oemInit()
ws.oemAddAtmosphereSpectral(target="extinction", matrix=np.eye(6) * 1.0)
ws.oemAddSurfaceSpectral(target="depolarization factor", matrix=np.eye(4) * 2.0)
ws.oemAddSubsurfaceSpectral(target="loss tangent", matrix=np.eye(4) * 3.0)
ws.measurement_vec = np.zeros(1)
ws.measurement_vec_error_covmat = arts.CovarianceMatrix()
try:
    ws.oemFinalizeDiagonal()
except RuntimeError:
    pass  # The check of the incomplete OEM data fails after the covariance is assembled
blocks = {(b.row_range.offset, b.row_range.extent): np.asarray(b.matrix) for b in ws.oem.model_state_covmat.blocks}
assert sorted(blocks) == [(0, 6), (6, 4), (10, 4)], sorted(blocks)
for (start, size), value in [((0, 6), 1.0), ((6, 4), 2.0), ((10, 4), 3.0)]:
    assert np.array_equal(blocks[(start, size)], value * np.eye(size))
