"""Outputs that a Python callback assigns or returns must reach the caller."""

import numpy as np
import pyarts3 as pyarts
from pyarts3.workspace import arts_agenda

arts = pyarts.arts

# %% A callback replacing an existing output, and creating a missing one

ws = pyarts.Workspace(False)
ws.measurement_jac = arts.Matrix(np.zeros((1, 1)))
assert not ws.has("measurement_vec_fit")


def replace(w):
    w.measurement_jac = arts.Matrix(np.ones((2, 2)))
    w.measurement_vec_fit = arts.Vector([1.0, 2.0])


arts.CallbackOperator(replace, [], ["measurement_jac", "measurement_vec_fit"])(ws)
np.testing.assert_array_equal(ws.measurement_jac, np.ones((2, 2)))
np.testing.assert_array_equal(ws.measurement_vec_fit, [1.0, 2.0])

# %% An agenda calling a function that returns its outputs

for name, group in (
    ("atm_field", arts.AtmField),
    ("abs_bands", arts.AbsorptionBands),
    ("measurement_sensors", arts.ArrayOfSensor),
    ("surf_field", arts.SurfaceField),
    ("subsurf_field", arts.SubsurfaceField),
    ("jac_targets", arts.JacobianTargets),
):
    setattr(ws, name, group())
ws.surf_field.ellipsoid = [2, 1]
ws.model_state_vec = [0.5, -0.25]
ws.atm_field[arts.AtmKey.t] = arts.GriddedField3(
    data=np.ones((2, 1, 1)),
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0, 1], [0], [0]],
)
ws.jac_targetsAddTemperature()
ws.jac_targetsFinalize()
ws.model_state_targets = ws.jac_targets

J = np.array([[1.0, 2.0], [2.0, -1.0], [1.0, 1.0]])


def forward(model_state_vec, jac_targets):
    measurement_vec_fit = arts.Vector(J @ np.asarray(model_state_vec))
    measurement_jac = arts.Matrix(J if jac_targets.x_size() else np.empty((0, 0)))
    return measurement_vec_fit, measurement_jac


@arts_agenda(ws=ws, fix=True)
def inversion_iterate_agenda(ws):
    forward()


ws.inversion_iterate_agendaExecute()
np.testing.assert_array_equal(ws.measurement_vec_fit, J @ [0.5, -0.25])
np.testing.assert_array_equal(ws.measurement_jac, J)
