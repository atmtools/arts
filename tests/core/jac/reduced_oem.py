"""ReducedOEM: full-basis equivalence and constrained-state retrieval.

Print results and timings, with plots unless ARTS_HEADLESS is set by CTest.
"""
import scipy.sparse
import os
from time import perf_counter

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
full_results = {}
basis_results = {None: {}, 1: {}, 2: {}}
for method in methods:
    full, reduced = workspace(), workspace()
    options = dict(method=method, max_iter=80, stop_dx=1e-10)
    full.OEM(**options)
    reduced.ReducedOEM(measurement_basis_mat=np.eye(
        3), model_state_basis_mat=basis, **options)
    assert reduced.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, reduced.oem_diagnostics
    np.testing.assert_allclose(reduced.model_state_vec, full.model_state_vec, atol=2e-5)
    np.testing.assert_allclose(reduced.measurement_gain_mat,
                               full.measurement_gain_mat, atol=1e-8)
    full_results[method] = (np.array(full.model_state_vec),
                            np.array(full.measurement_vec_fit))

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
ws.ReducedOEM(measurement_basis_mat=np.eye(3), model_state_basis_mat=b, method="li")
np.testing.assert_allclose(ws.model_state_vec, xa + b @ z, atol=1e-10)
for invalid in (np.zeros((2, 1)), np.full((2, 1), np.nan), np.ones((2, 2)), np.ones((3, 1))):
    try:
        workspace().ReducedOEM(measurement_basis_mat=np.eye(
            3), model_state_basis_mat=invalid, method="li")
    except RuntimeError:
        pass
    else:
        raise AssertionError("Invalid basis accepted")

# Nonlinear LM with value-only trial evaluations.
ws = workspace(True)
ws.ReducedOEM(measurement_basis_mat=[[1.]], model_state_basis_mat=[
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
    reduced.ReducedOEM(model_state_basis_mat=B, measurement_basis_mat=C,
                       method=method, max_iter=1)
    np.testing.assert_allclose(reduced.model_state_vec,
                               full.model_state_vec, atol=1e-10)
    assert np.linalg.norm(np.asarray(full.model_state_vec) - xa) > 1e-3
    assert reduced.oem_diagnostics.iterations == full.oem_diagnostics.iterations == 1

ws = workspace()
ws.model_state_vec = xa + np.array([0., 1.])
try:
    ws.ReducedOEM(measurement_basis_mat=np.eye(3), model_state_basis_mat=b, method="li")
except RuntimeError as error:
    assert "affine subspace" in str(error)
else:
    raise AssertionError("Out-of-subspace starting state accepted")
ws = workspace()
ws.ReducedOEM(measurement_basis_mat=np.eye(
    3), model_state_basis_mat=b, method="gn", clear_matrices=1)
assert np.asarray(ws.measurement_jac).size == 0
assert np.asarray(ws.measurement_gain_mat).size == 0
np.testing.assert_allclose(ws.model_state_vec, xa + b @ z, atol=1e-8)

# An actual null-space retrieval, using a basis chosen by the loss threshold.
null_j = np.array([[2., 0.], [0., 0.], [0., 0.]])
null_basis = information(null_j, np.eye(2), np.eye(
    3)).reduction(max_lost_dofs=0).model_state_basis_mat
assert null_basis.shape == (2, 1)
for method in ("li", "gn", "lm", "lm_cg"):
    pair = [workspace(evaluator=lambda state: (null_j @ state, null_j))
            for _ in range(2)]
    for item in pair:
        item.model_state_covmat = covariance(np.eye(2))
        item.measurement_vec_error_covmat = covariance(np.eye(3))
    pair[0].OEM(method=method, max_iter=80, stop_dx=1e-10)
    pair[1].ReducedOEM(measurement_basis_mat=np.eye(
        3), model_state_basis_mat=null_basis, method=method, max_iter=80, stop_dx=1e-10)
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
curved_basis = information(jac_at_prior, sa, se).reduction(rank=1).model_state_basis_mat
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
    ws.ReducedOEM(measurement_basis_mat=np.eye(
        3), model_state_basis_mat=curved_basis, method=method, max_iter=80, stop_dx=1e-12)
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
b = reduction.model_state_basis_mat
a = j @ b
reduced_precision = np.eye(1) + a.T @ np.linalg.solve(se, a)
z = np.linalg.solve(reduced_precision, a.T @
                    np.linalg.solve(se, y - j @ xa - [.25, -.5, 1]))
gain = b @ np.linalg.solve(reduced_precision, a.T @ np.linalg.inv(se))
for method in methods:
    ws = workspace()
    ws.ReducedOEM(measurement_basis_mat=reduction.measurement_basis_mat,
                  model_state_basis_mat=b, method=method, max_iter=80, stop_dx=1e-10)
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

# Generate both inputs entirely through workspace methods and consume them as
# workspace variables. No Python basis construction or explicit inputs to OEM.
for rank in (None, 1, 2):
    for method in methods:
        ws = workspace()
        ws.measurement_jac = j
        ws.ReducedOEMBasisCalc()
        if rank is not None:
            ws.ReducedOEMBasisReduce(rank=rank)
        B = np.array(ws.model_state_basis_mat)
        C = np.array(ws.measurement_basis_mat)
        jr = C @ j @ B
        expected_gain = B @ np.linalg.solve(np.eye(B.shape[1]) + jr.T @ jr, jr.T) @ C
        expected_state = xa + expected_gain @ (y - j @ xa - [.25, -.5, 1])
        ws.ReducedOEM(method=method, max_iter=80, stop_dx=1e-10)
        np.testing.assert_allclose(ws.model_state_vec,
                                   expected_state, atol=2e-5)
        np.testing.assert_allclose(ws.measurement_gain_mat, expected_gain, atol=1e-10)
        basis_results[rank][method] = (
            np.array(ws.model_state_vec), np.array(ws.measurement_vec_fit),
            expected_state, j @ expected_state + [.25, -.5, 1],
        )

# Dropping an exactly unobservable mode preserves the solution and uncertainty.
for method in methods:
    pair = [workspace(evaluator=lambda state: (null_j @ state, null_j))
            for _ in range(2)]
    for item in pair:
        item.model_state_covmat = covariance(np.eye(2))
        item.measurement_vec_error_covmat = covariance(np.eye(3))
    pair[0].OEM(method=method, max_iter=80, stop_dx=1e-10)
    pair[1].measurement_jac = null_j
    pair[1].ReducedOEMBasisCalc()
    pair[1].ReducedOEMBasisReduce()  # Remove the zero-information mode automatically.
    assert np.asarray(pair[1].model_state_basis_mat).shape == (2, 1)
    pair[1].ReducedOEM(method=method, max_iter=80, stop_dx=1e-10)
    np.testing.assert_allclose(pair[1].model_state_vec,
                               pair[0].model_state_vec, atol=2e-5)
    np.testing.assert_allclose(
        pair[1].measurement_gain_mat, pair[0].measurement_gain_mat, atol=1e-10)

# Discard a weak nonzero mode using its total information budget, then retrieve.
loss_limit = report.reduction(rank=1).discarded_information_bits + 1e-8
for method in methods:
    ws = workspace()
    ws.measurement_jac = j
    ws.ReducedOEMBasisCalc()
    ws.ReducedOEMBasisReduce(max_lost_information_bits=loss_limit)
    assert np.asarray(ws.model_state_basis_mat).shape == (2, 1)
    ws.ReducedOEM(method=method, max_iter=80, stop_dx=1e-10)
    np.testing.assert_allclose(ws.model_state_vec, xa + b @ z, atol=2e-5)
    np.testing.assert_allclose(ws.measurement_gain_mat, gain, atol=1e-10)

for method in ("li", "lm"):
    ws = workspace(evaluator=lambda state: (np.zeros(3), np.zeros((3, 2))))
    ws.measurement_jac = np.zeros((3, 2))
    ws.ReducedOEMBasisCalc()
    ws.ReducedOEMBasisReduce()
    ws.ReducedOEM(method=method, max_iter=80, stop_dx=1e-10)
    np.testing.assert_allclose(ws.model_state_vec, xa, atol=1e-12)
    np.testing.assert_allclose(ws.measurement_gain_mat, np.zeros((2, 3)), atol=1e-12)

# Dropping an exactly unobservable mode preserves the solution and uncertainty.
null_report = information([[2., 0.]], np.eye(2), np.eye(1))
null_reduction = null_report.reduction(max_lost_dofs=0)
assert null_reduction.rank == 1
np.testing.assert_allclose(null_reduction.posterior_covariance(), np.diag([.2, 1.]))

# A start-cost exit at the prior should use its already evaluated full Jacobian.
calls = []
ws = workspace(calls=calls)
ws.ReducedOEM(measurement_basis_mat=np.eye(3), model_state_basis_mat=b,
              method="lm", max_start_cost=1e-100)
assert ws.oem_diagnostics.status == arts.OptimalEstimationStatus.StartCostLimit
assert len(calls) == 1
np.testing.assert_allclose(ws.model_state_vec, xa)
np.testing.assert_allclose(ws.measurement_jac, j)

# LI's full-state restoration can fail after the solver itself succeeds.
ws = workspace(fail_at=lambda state: not np.array_equal(state, xa))
ws.ReducedOEM(measurement_basis_mat=np.eye(3), model_state_basis_mat=b, method="li")
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
        ws.ReducedOEM(model_state_basis_mat=B, measurement_basis_mat=C,
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
    ws.ReducedOEM(model_state_basis_mat=B, measurement_basis_mat=C,
                  method=method, max_iter=80, stop_dx=1e-10, **normalization)
    np.testing.assert_allclose(ws.model_state_vec, affine_reference(B, C)[0], atol=2e-5)

# Full-state information survives measurement compression for this affine model.
complete = report.reduction(rank=2)
assert complete.measurement_basis_mat.shape == (2, 3)
for method in methods:
    full, reduced = workspace(), workspace()
    full.OEM(method=method, max_iter=80, stop_dx=1e-10)
    reduced.ReducedOEM(model_state_basis_mat=complete.model_state_basis_mat,
                       measurement_basis_mat=complete.measurement_basis_mat,
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
    ws.ReducedOEM(model_state_basis_mat=B, measurement_basis_mat=C,
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
            model_state_basis_mat=np.eye(2), measurement_basis_mat=C, method="li")
    except RuntimeError:
        assert not calls
    else:
        raise AssertionError("Invalid measurement reduction accepted")

# Both matrices are required, even when one side uses identity.
for arguments in (dict(model_state_basis_mat=np.eye(2)), dict(measurement_basis_mat=np.eye(3))):
    try:
        workspace().ReducedOEM(method="li", **arguments)
    except (RuntimeError, TypeError):
        pass
    else:
        raise AssertionError("Missing reduction matrix accepted")


def test_grouped_measurements():
    from scipy import sparse

    jac = np.array([[1., 2.], [2., -1.], [-2., -4.]])
    offset = np.array([.25, -.5, 1.])
    for noise in (np.diag([1., 2., .5]), se):
        for method in methods:
            pair = [workspace(evaluator=lambda x: (jac @ x + offset, jac))
                    for _ in range(2)]
            for item in pair:
                item.measurement_vec_error_covmat = covariance(noise)
            full, reduced = pair
            full.OEM(method=method, max_iter=80, stop_dx=1e-10)
            reduced.measurement_jac = jac
            reduced.model_state_basis_mat = np.eye(2)
            reduced.measurement_basis_matCalc()
            reduced.ReducedOEM(method=method, max_iter=80, stop_dx=1e-10)
            assert reduced.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, reduced.oem_diagnostics
            np.testing.assert_allclose(
                reduced.model_state_vec, full.model_state_vec, atol=2e-5)
            np.testing.assert_allclose(
                reduced.measurement_gain_mat, full.measurement_gain_mat, atol=1e-10)
            np.testing.assert_allclose(
                reduced.measurement_vec_fit, full.measurement_vec_fit, atol=4e-5)

    # Explicit sparse C with off-diagonal noise blocks exercises sparse covariance
    # projection, including the implicit transpose of cross-correlation blocks.
    for storage in (arts.Matrix, lambda a: arts.Sparse(sparse.csr_matrix(a))):
        ws = workspace()
        ws.measurement_vec_error_covmat = arts.CovarianceMatrix()
        ws.measurement_vec_error_covmat.blocks = [
            arts.Block(arts.Range(0, 1), arts.Range(0, 1), (0, 0), storage(se[:1, :1])),
            arts.Block(arts.Range(1, 2), arts.Range(1, 2), (1, 1), storage(se[1:, 1:])),
            arts.Block(arts.Range(0, 1), arts.Range(1, 2), (0, 1), storage(se[:1, 1:])),
        ]
        c = np.array([[.5, 0., .5], [0., 1., 0.]])
        ws.ReducedOEM(model_state_basis_mat=np.eye(2),
                      measurement_basis_mat=arts.Sparse(sparse.csr_matrix(c)), method="li")
        state, gain = affine_reference(np.eye(2), c)
        np.testing.assert_allclose(ws.model_state_vec, state, atol=1e-10)
        np.testing.assert_allclose(ws.measurement_gain_mat, gain, atol=1e-10)

    # Many independent groups: storing either C or its projected diagonal noise
    # densely would require large quadratic arrays. The useful Jacobian stays thin.
    count = 10000
    jac = np.repeat(np.column_stack(
        (np.ones(count // 2), np.arange(count // 2))), 2, axis=0)
    jac[:, 1] /= count
    ws = workspace(evaluator=lambda x: (jac @ x, jac))
    ws.measurement_vec = jac @ [1., -.5]
    ws.measurement_jac = jac
    ws.model_state_basis_mat = np.eye(2)
    ws.measurement_vec_error_covmat = arts.CovarianceMatrix()
    ws.measurement_vec_error_covmat.blocks = [arts.Block(
        arts.Range(0, count), arts.Range(0, count), (0, 0),
        arts.Sparse(sparse.eye(count, format="csr")))]
    ws.measurement_basis_matCalc()
    assert ws.measurement_basis_mat.shape == (count // 2, count)
    assert ws.measurement_basis_mat.matrix.tocsr().nnz == count
    ws.ReducedOEM(method="li")
    expected = np.linalg.solve(np.linalg.inv(sa) + jac.T @ jac,
                               np.linalg.solve(sa, xa) + jac.T @ np.asarray(ws.measurement_vec))
    np.testing.assert_allclose(ws.model_state_vec, expected, atol=1e-10)


test_grouped_measurements()


def print_timings(method="lm", repeats=20):
    """Time complete workspace calls on the same small affine problem."""
    prepared = workspace()
    prepared.measurement_jac = j
    started = perf_counter()
    prepared.ReducedOEMBasisCalc()
    basis_seconds = perf_counter() - started

    cases = [("OEM (2 states, 3 measurements)", workspace(), "OEM", full_results[method])]
    for rank in (None, 2, 1):
        ws = workspace()
        ws.model_state_basis_mat = np.array(prepared.model_state_basis_mat)
        ws.measurement_basis_mat = np.array(prepared.measurement_basis_mat)
        ws.oem_basis_singular_values = prepared.oem_basis_singular_values
        if rank is not None:
            ws.ReducedOEMBasisReduce(rank=rank)
        states = np.asarray(ws.model_state_basis_mat).shape[1]
        measurements = np.asarray(ws.measurement_basis_mat).shape[0]
        label = (f"ReducedOEM ({states} state{'s' if states != 1 else ''}, "
                 f"{measurements} measurement{'s' if measurements != 1 else ''})")
        cases.append((label, ws, "ReducedOEM", basis_results[rank][method][:2]))

    options = dict(method=method, max_iter=80, stop_dx=1e-10)

    def retrieve(case):
        _, ws, name, expected = case
        # Reset outside the timer; each call must recompute its fit and Jacobian.
        ws.model_state_vec = []
        ws.measurement_vec_fit = []
        ws.measurement_jac = arts.Matrix()
        solve = getattr(ws, name)
        started = perf_counter()
        solve(**options)
        elapsed = perf_counter() - started
        assert ws.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, ws.oem_diagnostics
        np.testing.assert_allclose(ws.model_state_vec, expected[0], atol=2e-5)
        np.testing.assert_allclose(ws.measurement_vec_fit, expected[1], atol=2e-5)
        return elapsed

    for case in cases:
        retrieve(case)  # Exclude one warm-up per case.
    timings = [[] for _ in cases]
    for repeat in range(repeats):
        # Rotate which case runs first to distribute ordering effects.
        for offset in range(len(cases)):
            index = (repeat + offset) % len(cases)
            timings[index].append(retrieve(cases[index]))

    print(
        f"\nAffine {method.upper()} comparison: {repeats} runs per case after one warm-up; rotating order.")
    print("Same prior and empty fit/Jacobian each run; complete OEM calls including forward model and gain.")
    print("Workspace/basis setup, input resets, result checks and plots are outside the timer.")
    print("This 2-state, 3-measurement problem mainly measures call overhead; see the wind case for physical timings.")
    print(
        f"Basis calculation from the existing Jacobian (one call): {1e3*basis_seconds:.3f} ms")
    baseline = float(np.median(timings[0]))
    for (label, _, _, (state, _)), times in zip(cases, timings):
        median = float(np.median(times))
        delta = median - baseline
        difference = np.max(abs(state - full_results[method][0]))
        print(f"  {label}: median {1e3*median:.3f} ms, "
              f"range {1e3*min(times):.3f}–{1e3*max(times):.3f} ms; "
              f"vs OEM {1e3*delta:+.3f} ms ({100*delta/baseline:+.1f}%)")
        print(
            f"    state={np.array2string(state, precision=7)}, max |state − OEM|={difference:.3g}")
    print("Positive timing differences mean slower. Rank 1 discards information and changes the fitted state.")
    # Timings are reported, not used as pass/fail thresholds.


def plot(full_results, basis_results, method="lm"):
    """Compare existing affine runs; the truncated case has its own reference."""
    import matplotlib.pyplot as plt

    cases = (
        (None, "Full bases: 2 states, 3 measurements", "tab:blue", "o"),
        (2, "Reduced to 2 states, 2 measurements", "tab:orange", "s"),
        (1, "Reduced to 1 state, 1 measurement", "tab:green", "^"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(13, 11))
    state_index, measurement_index = np.arange(len(xa)), np.arange(len(y))
    full_state, full_fit = full_results[method]
    axes[0, 0].plot(state_index, xa, "k:", label="A priori")
    axes[0, 0].plot(state_index, full_state, "kx-", label="Full OEM")
    axes[0, 1].plot(measurement_index, y, "k+", ms=10, label="Measurements")
    axes[0, 1].plot(measurement_index, full_fit, "kx-", label="Full OEM")
    for rank, label, color, marker in cases:
        state, fit, _, _ = basis_results[rank][method]
        for ax, index, values in (
            (axes[0, 0], state_index, state),
            (axes[0, 1], measurement_index, fit),
            (axes[1, 0], state_index, state - full_state),
            (axes[1, 1], measurement_index, fit - full_fit),
        ):
            ax.plot(index, values, color=color, marker=marker, fillstyle="none",
                    linestyle="--", label=label)
        # Each case is checked against its own independent analytic solution.
        # Nonzero differences from full OEM for rank 1 are truncation, not errors.
        state_error, fit_error = [], []
        for solver in methods:
            state, fit, reference_state, reference_fit = basis_results[rank][solver]
            state_error.append(np.max(abs(state - reference_state)))
            fit_error.append(np.max(abs(fit - reference_fit)))
        for ax, values in zip(axes[2], (state_error, fit_error)):
            ax.plot(np.arange(len(methods)), values, color=color, marker=marker,
                    fillstyle="none", linestyle="--", label=label)

    axes[0, 0].set_title(f"Retrieved state ({method.upper()})")
    axes[0, 1].set_title(f"Fitted measurements ({method.upper()})")
    axes[0, 0].set_ylabel("Synthetic state value")
    axes[0, 1].set_ylabel("Synthetic measurement value")
    for col, (index, coordinate) in enumerate(((state_index, "State"),
                                              (measurement_index, "Measurement"))):
        for row in (0, 1):
            axes[row, col].set_xticks(index)
            axes[row, col].set_xlabel(f"{coordinate} index")
        axes[1, col].axhline(0, color="k", linewidth=.5)
        axes[1, col].set_title("Difference from full OEM")
        axes[1, col].set_ylabel(f"ReducedOEM − OEM {coordinate.lower()}")
        axes[2, col].set_title("Error against each case's analytic solution")
        axes[2, col].set_ylabel(f"Maximum absolute {coordinate.lower()} error")
        axes[2, col].set_ylim(bottom=0)
        axes[2, col].ticklabel_format(
            axis="y", style="sci", scilimits=(0, 0), useOffset=False)
        axes[2, col].set_xticks(np.arange(len(methods)),
                                methods, rotation=45, ha="right")
        axes[2, col].set_xlabel("OEM method")
    for ax in axes.flat:
        ax.grid(True, alpha=.3)
    axes[0, 0].legend(fontsize="small")
    axes[0, 1].legend(fontsize="small")
    fig.suptitle("ReducedOEM affine regression: complete bases and rank truncation")
    fig.tight_layout()
    return fig


# Compare full retrievals and individual damped steps, with nonzero starts.
for state_basis in (np.eye(2), np.diag([2., .5]),
                    np.array([[1., .3], [.2, 1.]]), np.array([[1.], [.4]])):
    for diagonal_prior in (False, True):
        for method in methods:
            for iterations in ((1, 80) if method in ("lm", "lm_cg") else (80,)):
                outputs = []
                for sparse in (False, True):
                    ws = workspace()
                    if diagonal_prior:
                        ws.model_state_covmat = arts.CovarianceMatrix()
                        ws.model_state_covmat.blocks = [arts.Block(
                            arts.Range(0, 2), arts.Range(0, 2), (0, 0),
                            arts.Sparse(scipy.sparse.diags([4., 2.]).tocsr()))]
                    ws.model_state_vec = np.asarray(
                        ws.model_state_vec_apriori) + state_basis @ np.full(state_basis.shape[1], .1)
                    basis = arts.Sparse(scipy.sparse.csr_matrix(
                        state_basis)) if sparse else state_basis
                    ws.ReducedOEM(model_state_basis_mat=basis, measurement_basis_mat=np.eye(3),
                                  method=method, max_iter=iterations, stop_dx=1e-10, clear_matrices=0)
                    assert ws.oem_diagnostics.status != arts.OptimalEstimationStatus.Error, ws.oem_diagnostics
                    outputs.append((np.array(ws.model_state_vec),
                                   np.array(ws.measurement_gain_mat)))
                for dense, sparse in zip(*outputs):
                    np.testing.assert_allclose(
                        dense, sparse, rtol=1e-8, atol=1e-9, err_msg=method)

# Truncating SVD-ordered sparse bases preserves their storage and values.
ws = workspace()
ws.model_state_basis_mat = arts.Sparse(scipy.sparse.eye(2).tocsr())
ws.measurement_basis_mat = arts.Sparse(scipy.sparse.eye(3).tocsr())
ws.oem_basis_singular_values = [2., 1.]
ws.ReducedOEMBasisReduce(rank=1)
assert ws.model_state_basis_mat.is_sparse
np.testing.assert_array_equal(np.array(ws.model_state_basis_mat), [[1.], [0.]])


if __name__ == "__main__":
    print_timings()
    if "ARTS_HEADLESS" not in os.environ:
        plot(full_results, basis_results)
        import matplotlib.pyplot as plt
        plt.show()
