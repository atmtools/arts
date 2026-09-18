import numpy as np
import pyarts3 as pyarts


# A general invertible affine transformation x = A (t - b).  A is
# intentionally neither diagonal nor orthogonal.
A = np.array([[2.0, 0.5], [-0.25, 1.5]])
A_inv = np.linalg.inv(A)
b = np.array([10.0, -3.0])


def transform_state(native_state, atm_field):
    assert pyarts.arts.AtmKey.t in atm_field
    return A @ (np.asarray(native_state) - b)


def inverse_state(model_state, atm_field):
    assert pyarts.arts.AtmKey.t in atm_field
    return A_inv @ np.asarray(model_state) + b


def inverse_jacobian(jacobian, model_state, atm_field):
    assert pyarts.arts.AtmKey.t in atm_field
    assert len(model_state) == 2
    return np.asarray(jacobian) @ A_inv


ws = pyarts.Workspace()
ws.measurement_sensor = []
ws.abs_bands = {}
ws.atm_fieldInit(toa=100e3)
ws.surf_fieldEarth()

native_state = np.array([12.0, 4.0])
ws.atm_field[pyarts.arts.AtmKey.t] = pyarts.arts.GriddedField3(
    data=native_state.reshape(2, 1, 1),
    grids=([0.0, 100e3], [0.0], [0.0]),
    grid_names=["Altitude", "Latitude", "Longitude"],
    name="Temperature",
)

ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
target = ws.jac_targets.atm[-1]
target.transform_state = transform_state
target.inverse_state = inverse_state
target.inverse_jacobian = inverse_jacobian
assert ws.jac_targets.atm[-1].transform_state is not None
ws.jac_targetsFinalize()

# Native data -> transformed model state.
ws.model_state_vecFromData()
expected_model_state = A @ (native_state - b)
assert np.allclose(ws.model_state_vec, expected_model_state), (
    f"model_state_vec={ws.model_state_vec}, expected={expected_model_state}"
)

# Transformed model state -> native data.
new_model_state = np.array([3.0, -2.0])
ws.model_state_vec = new_model_state
# The full state mapping is explicit; Jacobian targets independently select derivatives.
ws.atm_fieldFromModelState(model_state_targets=ws.jac_targets)
expected_native_state = A_inv @ new_model_state + b
assert np.allclose(np.asarray(ws.atm_field[pyarts.arts.AtmKey.t].data).flat, expected_native_state)

# Native Jacobian -> transformed-state Jacobian.
native_jacobian = np.array([[1.0, 2.0], [-4.0, 3.0], [0.5, -1.0]])
ws.measurement_jac = native_jacobian
ws.measurement_jacAtmosphereTransformation()
assert np.allclose(ws.measurement_jac, native_jacobian @ A_inv)

# The bound operators are readable as well as writable.
assert np.allclose(target.transform_state(expected_native_state, ws.atm_field), new_model_state)
assert np.allclose(target.inverse_state(new_model_state, ws.atm_field), expected_native_state)
assert np.allclose(
    target.inverse_jacobian(native_jacobian, new_model_state, ws.atm_field),
    native_jacobian @ A_inv,
)
