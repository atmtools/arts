"""
Sensor frequency offset Jacobians must land in their own model state columns,
also when other targets come before them, and every polynomial coefficient
must get its own column.
"""

import pyarts3 as pyarts
import numpy as np

NF = 21
D = 1e3
line_f0 = 118750348044.712
f = np.linspace(-5e6, 5e6, NF) + line_f0
pos = [100e3, 0, 0]
los = [180.0, 0.0]

ws = pyarts.workspace.Workspace()
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)
ws.WignerInit()
ws.spectral_propmat_agendaAuto()
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.spectral_rad_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()


def set_sensor(freq_grid):
    ws.measurement_sensorSimpleGaussian(
        freq_grid=freq_grid, std=1, pos=pos, los=los, pol="I"
    )


def jacobian(with_t, polyorder=None):
    set_sensor(f)
    ws.jac_targetsInit()
    if with_t:
        ws.jac_targetsAddTemperature()
    if polyorder is not None:
        ws.jac_targetsAddSensorFrequencyPolyOffset(
            sensor_elem=0, d=D, polyorder=polyorder
        )
    ws.jac_targetsFinalize()
    ws.measurement_vecFromSensor()
    return np.array(ws.measurement_vec), np.array(ws.measurement_jac)


def assert_close(a, b, rtol, msg):
    scale = np.abs(b).max()
    assert scale > 0, f"{msg}: reference is all zeros"
    err = np.abs(a - b).max() / scale
    assert err < rtol, f"{msg}: relative error {err}"


# %% Placement: a temperature target before the sensor target

_, J_sensor = jacobian(False, 0)
_, J_t = jacobian(True)
_, J_both = jacobian(True, 0)

nt = J_t.shape[1]
assert J_both.shape[1] == nt + 1

assert_close(J_both[:, nt], J_sensor[:, 0], 1e-6, "Sensor column misplaced")
assert_close(J_both[:, :nt], J_t, 1e-6, "Sensor derivative leaked into temperature")

# %% Polynomial order 1: f_i -> f_i + c0 + c1 * t_i, with t the grid mapped to [-1, 1]

_, J_poly = jacobian(False, 1)
assert J_poly.shape[1] == 2

assert_close(J_poly[:, 0], J_sensor[:, 0], 1e-6, "Constant term changed")

t = (f - 0.5 * (f.max() + f.min())) / (0.5 * (f.max() - f.min()))
set_sensor(f)
ws.jac_targetsInit()
ws.measurement_vecFromSensor()
y0 = np.array(ws.measurement_vec)
set_sensor(f + D * t)
ws.measurement_vecFromSensor()
y1 = np.array(ws.measurement_vec)

# The shifts here have both signs while the Jacobian uses forward differences,
# so the two differ at second order
assert_close(J_poly[:, 1], (y1 - y0) / D, 1e-2, "Linear term wrong")

# %% The polynomial order must be less than the number of grid points

set_sensor(f)
ws.jac_targetsInit()
try:
    ws.jac_targetsAddSensorFrequencyPolyOffset(sensor_elem=0, d=D, polyorder=NF)
except Exception:
    pass
else:
    raise AssertionError("Accepted a polynomial order the grid cannot determine")
