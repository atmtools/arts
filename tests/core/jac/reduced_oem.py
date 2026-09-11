"""ReducedOEM: full-basis equivalence and constrained-state retrieval."""
import numpy as np
import pyarts3 as pyarts
from pyarts3.retrieval import information, information_from_workspace
arts = pyarts.arts


def covariance(values):
    matrix = arts.Matrix(values)
    size = np.asarray(matrix).shape[0]
    result = arts.CovarianceMatrix()
    result.blocks = [
        arts.Block(arts.Range(0, size), arts.Range(0, size), (0, 0), matrix)
    ]
    return result


def workspace(nonlinear=False, *, calls=None, fail_at=None, evaluator=None):
    ws = pyarts.Workspace()
    for variable, group in (
        ("model_state_vec", arts.Vector),
        ("measurement_vec_fit", arts.Vector),
        ("measurement_jac", arts.Matrix),
        ("atm_field", arts.AtmField),
        ("abs_bands", arts.AbsorptionBands),
        ("measurement_sensor", arts.ArrayOfSensorObsel),
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
        ws.model_state_vec = [0.1]
        ws.model_state_vec_apriori = [1]
        ws.measurement_vec = [4]
        ws.model_state_covmat = covariance([[4]])
        ws.measurement_vec_error_covmat = covariance([[0.25]])
    else:
        ws.model_state_vec_apriori = [0.5, -0.25]
        ws.measurement_vec = [2, -1, 1.5]
        ws.model_state_covmat = covariance([[4, 1], [1, 2]])
        ws.measurement_vec_error_covmat = covariance(
            [[1, 0.2, 0], [0.2, 2, 0.3], [0, 0.3, 0.5]]
        )
    jacobian = np.array([[1, 2], [2, -1], [1, 1]], dtype=float)

    def forward(local):
        state = np.asarray(local.get("model_state_vec"))
        if calls is not None:
            calls.append((state.copy(), local.get("jac_targets").x_size() != 0))
        if fail_at is not None and fail_at(state):
            raise RuntimeError("deliberate reduced agenda failure")
        if nonlinear:
            fit = np.array([state[0] ** 2])
            derivative = np.array([[2 * state[0]]])
        else:
            fit = jacobian @ state + [0.25, -0.5, 1]
            derivative = jacobian
        if evaluator is not None:
            fit, derivative = evaluator(state)
        if not local.get("jac_targets").x_size():
            derivative = np.empty((0, 0))
        # Preserve the existing output objects shared with the outer workspace.
        local.get("measurement_vec_fit").value = arts.Vector(fit)
        local.get("measurement_jac").value = arts.Matrix(derivative)

    callback = arts.CallbackOperator(
        forward,
        ["model_state_vec", "jac_targets"],
        ["measurement_vec_fit", "measurement_jac"],
    )
    agenda = arts.Agenda("inversion_iterate_agenda")
    agenda.add(arts.Method("synthetic_forward_model", callback))
    agenda.finalize(True)
    ws.inversion_iterate_agenda = agenda
    return ws


sa = np.array([[4., 1.], [1., 2.]])
basis = np.linalg.cholesky(sa)
methods = ("li", "li_m", "li_cg", "li_cg_m", "gn",
           "gn_m", "gn_cg", "gn_cg_m", "lm", "lm_cg", "ml", "ml_cg")
for method in methods:
    full, reduced = workspace(), workspace()
    options = dict(method=method, max_iter=80, stop_dx=1e-10)
    full.OEM(**options)
    reduced.ReducedOEM(measurement_red_mat=np.eye(
        3), model_state_red_mat=basis, **options)
    assert reduced.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, reduced.oem_diagnostics
    np.testing.assert_allclose(reduced.model_state_vec, full.model_state_vec, atol=2e-5)
    np.testing.assert_allclose(reduced.measurement_gain_mat,
                               full.measurement_gain_mat, atol=1e-8)

# Independent scalar posterior for a restricted affine state.
ws = workspace()
b = basis[:, :1]
j = np.array([[1., 2.], [2., -1.], [1., 1.]])
se = np.array([[1., .2, 0], [.2, 2., .3], [0, .3, .5]])
xa = np.asarray(ws.model_state_vec_apriori).copy()
y = np.asarray(ws.measurement_vec).copy()
a = j @ b
z = np.linalg.solve(np.eye(1) + a.T @ np.linalg.solve(se, a), a.T @
                    np.linalg.solve(se, y - j @ xa - [.25, -.5, 1]))
ws.ReducedOEM(measurement_red_mat=np.eye(3), model_state_red_mat=b, method="li")
np.testing.assert_allclose(ws.model_state_vec, xa + b @ z, atol=1e-10)
for invalid in (np.zeros((2, 1)), np.full((2, 1), np.nan), np.ones((2, 2)), np.ones((3, 1))):
    try:
        workspace().ReducedOEM(measurement_red_mat=np.eye(
            3), model_state_red_mat=invalid, method="li")
    except RuntimeError:
        pass
    else:
        raise AssertionError("Invalid basis accepted")

# Nonlinear LM with value-only trial evaluations.
ws = workspace(True)
ws.ReducedOEM(measurement_red_mat=[[1.]], model_state_red_mat=[
              [2.]], method="lm", max_iter=80, stop_dx=1e-10)
reference = workspace(True)
reference.OEM(method="lm", max_iter=80, stop_dx=1e-10)
np.testing.assert_allclose(ws.model_state_vec, reference.model_state_vec, atol=2e-5)

# A single damped step verifies the projected LM damping, not just its limit.
for method in ("lm", "lm_cg"):
    full, reduced = workspace(evaluator=lambda state: (
        j @ state, j)), workspace(evaluator=lambda state: (j @ state, j))
    B = np.array([[1.7, -.2], [.3, .8]])
    C = np.array([[1., .2, 0], [0, 2., .3], [.1, 0, 1.]])
    full.OEM(method=method, max_iter=1)
    reduced.ReducedOEM(model_state_red_mat=B, measurement_red_mat=C,
                       method=method, max_iter=1)
    np.testing.assert_allclose(reduced.model_state_vec,
                               full.model_state_vec, atol=1e-10)
    assert np.linalg.norm(np.asarray(full.model_state_vec) - xa) > 1e-3
    assert reduced.oem_diagnostics.iterations == full.oem_diagnostics.iterations == 1

ws = workspace()
ws.model_state_vec = xa + np.array([0., 1.])
try:
    ws.ReducedOEM(measurement_red_mat=np.eye(3), model_state_red_mat=b, method="li")
except RuntimeError as error:
    assert "affine subspace" in str(error)
else:
    raise AssertionError("Out-of-subspace starting state accepted")
ws = workspace()
ws.ReducedOEM(measurement_red_mat=np.eye(
    3), model_state_red_mat=b, method="gn", clear_matrices=1)
assert np.asarray(ws.measurement_jac).size == 0
assert np.asarray(ws.measurement_gain_mat).size == 0
np.testing.assert_allclose(ws.model_state_vec, xa + b @ z, atol=1e-8)

# An actual null-space retrieval, using a basis chosen by the loss threshold.
null_j = np.array([[2., 0.], [0., 0.], [0., 0.]])
null_basis = information(null_j, np.eye(2), np.eye(
    3)).reduction(max_lost_dofs=0).model_state_red_mat
assert null_basis.shape == (2, 1)
for method in ("li", "gn", "lm", "lm_cg"):
    pair = [workspace(evaluator=lambda state: (null_j @ state, null_j))
            for _ in range(2)]
    for item in pair:
        item.model_state_covmat = covariance(np.eye(2))
        item.measurement_vec_error_covmat = covariance(np.eye(3))
    pair[0].OEM(method=method, max_iter=80, stop_dx=1e-10)
    pair[1].ReducedOEM(measurement_red_mat=np.eye(
        3), model_state_red_mat=null_basis, method=method, max_iter=80, stop_dx=1e-10)
    np.testing.assert_allclose(pair[1].model_state_vec,
                               pair[0].model_state_vec, atol=2e-5)
    np.testing.assert_allclose(
        pair[1].measurement_gain_mat, pair[0].measurement_gain_mat, atol=1e-10)
    assert np.asarray(pair[1].model_state_vec)[1] == xa[1]

# Nonlinear truncation to one coefficient. The independent oracle solves the
# resulting quartic cost's cubic derivative and chooses its global minimum.


def curved(state):
    fit = j @ state + [.25 + .1 * state[0]**2, -.5 + .1 * state[1]**2, 1]
    derivative = j + [[.2 * state[0], 0], [0, .2 * state[1]], [0, 0]]
    return fit, derivative


fit_at_prior, jac_at_prior = curved(xa)
curved_basis = information(jac_at_prior, sa, se).reduction(rank=1).model_state_red_mat
direction = curved_basis[:, 0]
d0, d1 = fit_at_prior - y, jac_at_prior @ direction
d2 = np.array([.1 * direction[0]**2, .1 * direction[1]**2, 0])
precision = np.linalg.inv(se)
cost = np.array([d0 @ precision @ d0,
                 2 * d1 @ precision @ d0,
                 1 + d1 @ precision @ d1 + 2 * d2 @ precision @ d0,
                 2 * d2 @ precision @ d1,
                 d2 @ precision @ d2])
roots = np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(cost))
real_roots = roots.real[np.abs(roots.imag) < 1e-10]
optimum = min(real_roots, key=lambda value: np.polynomial.polynomial.polyval(value, cost))
for method in ("gn", "gn_m", "gn_cg", "gn_cg_m", "lm", "lm_cg"):
    calls = []
    ws = workspace(evaluator=curved, calls=calls)
    ws.ReducedOEM(measurement_red_mat=np.eye(
        3), model_state_red_mat=curved_basis, method=method, max_iter=80, stop_dx=1e-12)
    assert ws.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, ws.oem_diagnostics
    state = np.asarray(ws.model_state_vec)
    np.testing.assert_allclose(state, xa + direction * optimum, atol=2e-5)
    np.testing.assert_allclose(ws.measurement_vec_fit, curved(state)[0], atol=1e-10)
    np.testing.assert_allclose(ws.measurement_jac, curved(state)[1], atol=1e-10)
    # The last successful Jacobian already represents the returned state.
    assert not (len(calls) > 1 and calls[-1][1] and calls[-2][1]
                and np.array_equal(calls[-1][0], calls[-2][0]))

# Use the actual information modes through the public workspace report.
probe = workspace()
probe.model_state_vec = probe.model_state_vec_apriori
probe.model_state_targets = probe.jac_targets
probe.inversion_iterate_agendaExecute()
report = information_from_workspace(probe)
reduction = report.reduction(rank=1)
b = reduction.model_state_red_mat
a = j @ b
reduced_precision = np.eye(1) + a.T @ np.linalg.solve(se, a)
z = np.linalg.solve(reduced_precision, a.T @
                    np.linalg.solve(se, y - j @ xa - [.25, -.5, 1]))
gain = b @ np.linalg.solve(reduced_precision, a.T @ np.linalg.inv(se))
for method in methods:
    ws = workspace()
    ws.ReducedOEM(measurement_red_mat=reduction.measurement_red_mat,
                  model_state_red_mat=b, method=method, max_iter=80, stop_dx=1e-10)
    assert ws.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, ws.oem_diagnostics
    np.testing.assert_allclose(ws.model_state_vec, xa + b @ z, atol=2e-5)
    np.testing.assert_allclose(ws.measurement_gain_mat, gain, atol=1e-10)
    np.testing.assert_allclose(ws.measurement_jac, j, atol=1e-12)
    np.testing.assert_allclose(ws.measurement_vec_fit, j @
                               np.asarray(ws.model_state_vec) + [.25, -.5, 1])
    # Existing workspace error-covariance methods preserve omitted uncertainty.
    ws.measurement_averaging_kernelCalc()
    ws.model_state_covmat_smoothing_errorCalc()
    ws.measurement_vec_error_covmat_observation_systemCalc()
    uncertainty = (np.asarray(ws.model_state_covmat_smoothing_error)
                   + np.asarray(ws.measurement_vec_error_covmat_observation_system))
    np.testing.assert_allclose(
        uncertainty, reduction.posterior_covariance(), atol=1e-10)

# Dropping an exactly unobservable mode preserves the solution and uncertainty.
null_report = information([[2., 0.]], np.eye(2), np.eye(1))
null_reduction = null_report.reduction(max_lost_dofs=0)
assert null_reduction.rank == 1
np.testing.assert_allclose(null_reduction.posterior_covariance(), np.diag([.2, 1.]))

# A start-cost exit at the prior should use its already evaluated full Jacobian.
calls = []
ws = workspace(calls=calls)
ws.ReducedOEM(measurement_red_mat=np.eye(3), model_state_red_mat=b,
              method="lm", max_start_cost=1e-100)
assert ws.oem_diagnostics.status == arts.OptimalEstimationStatus.StartCostLimit
assert len(calls) == 1
np.testing.assert_allclose(ws.model_state_vec, xa)
np.testing.assert_allclose(ws.measurement_jac, j)

# LI's full-state restoration can fail after the solver itself succeeds.
ws = workspace(fail_at=lambda state: not np.array_equal(state, xa))
ws.ReducedOEM(measurement_red_mat=np.eye(3), model_state_red_mat=b, method="li")
assert ws.oem_diagnostics.status == arts.OptimalEstimationStatus.Error
assert "deliberate reduced agenda failure" in str(ws.oem_diagnostics.errors)
assert np.all(np.isnan(ws.model_state_vec))
assert np.asarray(ws.measurement_gain_mat).size == 0


def affine_reference(B, C):
    """Independent restricted Gaussian posterior in arbitrary coordinates."""
    jr = C @ j @ B
    sr = C @ se @ C.T
    prior_precision = B.T @ np.linalg.solve(sa, B)
    gain = B @ np.linalg.solve(prior_precision + jr.T @ np.linalg.solve(sr, jr),
                               np.linalg.solve(sr, jr).T) @ C
    state = xa + gain @ (y - j @ xa - [.25, -.5, 1])
    return state, gain


# Exercise both projections, without relying on prior/noise normalization.
# This includes channel selection and correlated noise after channel mixing.
for B, C in (
    (np.array([[2., .3], [-.4, .7]]),
     np.array([[1., .3, 0], [0, 2., -.2], [.2, 0, 1.]])),
    (np.eye(2), np.array([[0., 1., 0.]])),
    (np.array([[.7], [-.3]]), np.array([[1., 0, .5], [0, 2., 1.]])),
    (np.array([[.7], [-.3]]), np.eye(3)),
):
    expected_state, expected_gain = affine_reference(B, C)
    for method in methods:
        calls = []
        ws = workspace(calls=calls)
        ws.ReducedOEM(model_state_red_mat=B, measurement_red_mat=C,
                      method=method, max_iter=80, stop_dx=1e-10)
        assert ws.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, ws.oem_diagnostics
        np.testing.assert_allclose(ws.model_state_vec, expected_state, atol=2e-5)
        np.testing.assert_allclose(ws.measurement_gain_mat, expected_gain, atol=1e-10)
        # The physical agenda must always receive the original two-element state.
        assert all(state.shape == (2,) for state, _ in calls)
        np.testing.assert_allclose(ws.measurement_jac, j, atol=1e-12)
        state = np.asarray(ws.model_state_vec)
        residual = y - j @ state - [.25, -.5, 1]
        measurement_cost = residual @ np.linalg.solve(se, residual) / len(y)
        delta = state - xa
        np.testing.assert_allclose(
            ws.oem_diagnostics.measurement_cost, measurement_cost, atol=1e-10)
        np.testing.assert_allclose(ws.oem_diagnostics.final_cost,
                                   measurement_cost + delta @ np.linalg.solve(sa, delta) / len(y), atol=1e-10)

# Reduced-size normalization vectors are optional, independently of B and C.
for method, normalization in (
    ("li_cg", dict(model_state_covmat_normalization=[3.])),
    ("gn_m", dict(measurement_vec_normalization=[2., 3.])),
    ("gn_cg_m", dict(measurement_vec_normalization=[2., 3.])),
    ("lm", dict(model_state_covmat_normalization=[3.])),
):
    B, C = np.array([[.7], [-.3]]), np.array([[1., 0, .5], [0, 2., 1.]])
    ws = workspace()
    ws.ReducedOEM(model_state_red_mat=B, measurement_red_mat=C,
                  method=method, max_iter=80, stop_dx=1e-10, **normalization)
    np.testing.assert_allclose(ws.model_state_vec, affine_reference(B, C)[0], atol=2e-5)

# Full-state information survives measurement compression for this affine model.
complete = report.reduction(rank=2)
assert complete.measurement_red_mat.shape == (2, 3)
for method in methods:
    full, reduced = workspace(), workspace()
    full.OEM(method=method, max_iter=80, stop_dx=1e-10)
    reduced.ReducedOEM(model_state_red_mat=complete.model_state_red_mat,
                       measurement_red_mat=complete.measurement_red_mat,
                       method=method, max_iter=80, stop_dx=1e-10)
    np.testing.assert_allclose(reduced.model_state_vec, full.model_state_vec, atol=2e-5)
    np.testing.assert_allclose(reduced.measurement_gain_mat,
                               full.measurement_gain_mat, atol=1e-10)

# A compressed nonlinear problem has its own objective: compare against its
# quartic scalar cost, including a non-unit prior precision and correlated noise.
B = 1.7 * curved_basis
C = np.array([[1., .3, 0], [0, .5, 1.]])
direction = B[:, 0]
d0, d1 = C @ (fit_at_prior - y), C @ jac_at_prior @ direction
d2 = C @ np.array([.1 * direction[0]**2, .1 * direction[1]**2, 0])
precision = np.linalg.inv(C @ se @ C.T)
prior_precision = (B.T @ np.linalg.solve(sa, B))[0, 0]
cost = np.array([d0 @ precision @ d0, 2 * d1 @ precision @ d0,
                 prior_precision + d1 @ precision @ d1 + 2 * d2 @ precision @ d0,
                 2 * d2 @ precision @ d1, d2 @ precision @ d2])
roots = np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(cost))
real_roots = roots.real[np.abs(roots.imag) < 1e-10]
optimum = min(real_roots, key=lambda value: np.polynomial.polynomial.polyval(value, cost))
for method in ("gn", "gn_m", "gn_cg", "gn_cg_m", "lm", "lm_cg"):
    calls = []
    ws = workspace(evaluator=curved, calls=calls)
    ws.ReducedOEM(model_state_red_mat=B, measurement_red_mat=C,
                  method=method, max_iter=80, stop_dx=1e-12)
    np.testing.assert_allclose(ws.model_state_vec, xa + direction * optimum, atol=2e-5)
    state = np.asarray(ws.model_state_vec)
    np.testing.assert_allclose(ws.measurement_vec_fit, curved(state)[0], atol=1e-10)
    np.testing.assert_allclose(ws.measurement_jac, curved(state)[1], atol=1e-10)
    if method.startswith("lm"):
        assert any(not derivatives for _, derivatives in calls)

# Validate projection shape, finiteness, and rank before executing an agenda.
for C in (np.zeros((1, 3)), np.ones((2, 3)), np.eye(2),
          np.full((1, 3), np.nan), np.empty((0, 3))):
    calls = []
    try:
        workspace(calls=calls).ReducedOEM(
            model_state_red_mat=np.eye(2), measurement_red_mat=C)
    except RuntimeError:
        assert not calls
    else:
        raise AssertionError("Invalid measurement reduction accepted")

# Both matrices are required, even when one side uses identity.
for arguments in (dict(model_state_red_mat=np.eye(2)), dict(measurement_red_mat=np.eye(3))):
    try:
        workspace().ReducedOEM(**arguments)
    except (RuntimeError, TypeError):
        pass
    else:
        raise AssertionError("Missing reduction matrix accepted")
