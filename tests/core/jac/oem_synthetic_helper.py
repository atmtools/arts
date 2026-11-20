"""The synthetic OEM problem shared by the OEM tests next to this file.

The affine problem has two states and three measurements with the forward
model K x + OFFSET; the nonlinear problem has one state and the forward model
x**2.  Only the target dimensions of the workspace are used by the forward
model.
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

K = np.array([[1, 2], [2, -1], [1, 1]], dtype=float)
OFFSET = np.array([0.25, -0.5, 1])
SA = np.array([[4, 1], [1, 2]], dtype=float)
SE = np.array([[1, 0.2, 0], [0.2, 2, 0.3], [0, 0.3, 0.5]])
XA = np.array([0.5, -0.25])
Y = np.array([2, -1, 1.5])


def workspace(nonlinear=False, *, calls=None, fail_at=None, evaluator=None):
    """A workspace with the OEM data and inversion_iterate_agenda of the problem

    calls collects (state, with Jacobian) of every evaluation, fail_at(state)
    makes an evaluation fail, and evaluator(state) -> (fit, jac) replaces the
    forward model.
    """
    ws = pyarts.Workspace()
    ws.oem = arts.OptimalEstimationData()
    for variable, group in (
        ("model_state_vec", arts.Vector),
        ("measurement_vec_fit", arts.Vector),
        ("measurement_jac", arts.Matrix),
        ("atm_field", arts.AtmField),
        ("abs_bands", arts.AbsorptionBands),
        ("measurement_sensors", arts.ArrayOfSensor),
        ("surf_field", arts.SurfaceField),
        ("subsurf_field", arts.SubsurfaceField),
        ("jac_targets", arts.JacobianTargets),
    ):
        setattr(ws, variable, group())

    # Only target dimensions are used by the synthetic forward model. Creating
    # tiny temperature grids satisfies the public agenda's target-size checks.
    size = 1 if nonlinear else 2
    ws.surf_field.ellipsoid = [1, 1]
    ws.atm_field[arts.AtmKey.t] = arts.GriddedField3(
        name="Synthetic state metadata",
        data=np.ones((size, 1, 1)),
        grid_names=["Altitude", "Latitude", "Longitude"],
        grids=[np.arange(size, dtype=float), [0], [0]],
    )
    ws.jac_targetsAddTemperature()
    ws.jac_targetsFinalize()

    if nonlinear:
        ws.oem.model_state_vec = [0.1]
        ws.oem.model_state_vec_apriori = [1]
        ws.oem.measurement_vec = [4]
        ws.oem.model_state_covmat = arts.CovarianceMatrix([[4]])
        ws.oem.measurement_vec_error_covmat = arts.CovarianceMatrix([[0.25]])
    else:
        ws.oem.model_state_vec_apriori = XA
        ws.oem.measurement_vec = Y
        ws.oem.model_state_covmat = arts.CovarianceMatrix(SA)
        ws.oem.measurement_vec_error_covmat = arts.CovarianceMatrix(SE)

    def forward(atm, bands, sensors, surf, subsurf, state_targets, targets, state):
        state = np.asarray(state)
        if calls is not None:
            calls.append((state.copy(), targets.x_size() != 0))
        if fail_at is not None and fail_at(state):
            raise RuntimeError("deliberate reduced agenda failure")
        if nonlinear:
            fit = np.array([state[0] ** 2])
            derivative = np.array([[2 * state[0]]])
        else:
            fit = K @ state + OFFSET
            derivative = K
        if evaluator is not None:
            fit, derivative = evaluator(state)
        if not targets.x_size():
            derivative = np.empty((0, 0))
        return atm, bands, sensors, surf, subsurf, arts.Vector(fit), arts.Matrix(derivative)

    ws.inversion_iterate_agendaSetOperator(f=forward)
    return ws
