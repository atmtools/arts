"""Independent linear Gaussian references for the read-only OEM setup report."""

import copy
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyarts3 as pyarts

from pyarts3.retrieval import information, information_from_workspace

arts = pyarts.arts
K = np.array([[1, 2], [2, -1], [1, 1]], dtype=float)
SA = np.array([[4, 1], [1, 2]], dtype=float)
SE = np.array([[1, 0.2, 0], [0.2, 2, 0.3], [0, 0.3, 0.5]])
# Exact rational Gaussian elimination of (Sa^-1 + K^T Se^-1 K)^-1.
# This oracle is independent of ARTS and the report's factorization/SVD.
POSTERIOR = np.array([[5494, -1563], [-1563, 2551]], dtype=float) / 18575


def rejects(operation):
    try:
        operation()
    except (RuntimeError, ValueError, TypeError, MemoryError) as error:
        assert str(error), "Validation error has no explanation"
    else:
        raise AssertionError("Invalid information report input was accepted")


def covariance(values):
    values = np.asarray(values, dtype=float)
    size = values.shape[0]
    result = arts.CovarianceMatrix()
    result.blocks = [
        arts.Block(
            arts.Range(0, size), arts.Range(0, size), (0, 0), arts.Matrix(values)
        )
    ]
    return result


def covariance_snapshot(value):
    # XML includes stored inverse blocks; the public blocks property alone
    # cannot detect a report silently creating or discarding the inverse cache.
    with TemporaryDirectory() as directory:
        path = Path(directory) / "covariance.xml"
        value.savexml(str(path))
        return path.read_text()


def assert_modes(report, jacobian, prior, error, posterior):
    modes = report.state_modes
    singular = report.singular_values
    n = prior.shape[0]
    m = error.shape[0]
    assert modes.shape == (n, n)
    assert singular.shape == (n,)
    assert report.measurement_modes.shape == (m, min(m, n))
    np.testing.assert_allclose(modes @ modes.T, prior, atol=1e-13, rtol=1e-12)
    whitened = np.linalg.solve(np.linalg.cholesky(prior), modes)
    np.testing.assert_allclose(whitened.T @ whitened, np.eye(n), atol=1e-12)
    np.testing.assert_allclose(
        (modes / (1 + singular**2)) @ modes.T,
        posterior,
        atol=1e-13,
        rtol=1e-12,
    )
    projected = jacobian @ modes
    np.testing.assert_allclose(
        projected.T @ np.linalg.solve(error, projected),
        np.diag(singular**2),
        atol=1e-11,
        rtol=1e-11,
    )
    np.testing.assert_allclose(
        report.measurement_modes.T @ report.measurement_modes,
        np.eye(min(m, n)),
        atol=1e-12,
    )


def test_correlated_reference():
    report = information(K, SA, SE, state_labels=["temperature", "water"])
    np.testing.assert_allclose(
        report.prior_standard_deviation, np.sqrt(np.diag(SA)), rtol=1e-13
    )
    np.testing.assert_allclose(
        report.posterior_standard_deviation,
        np.sqrt(np.diag(POSTERIOR)),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        report.variance_reduction,
        1 - np.diag(POSTERIOR) / np.diag(SA),
        atol=1e-13,
    )
    # Degrees of freedom is also trace(G K), independent of the SVD route.
    gain = POSTERIOR @ K.T @ np.linalg.inv(SE)
    np.testing.assert_allclose(
        report.degrees_of_freedom, np.trace(gain @ K), rtol=1e-12
    )
    np.testing.assert_allclose(
        report.information_bits,
        0.5 * np.log2(np.linalg.det(SA) / np.linalg.det(POSTERIOR)),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        report.mode_variance_reduction,
        report.singular_values**2 / (1 + report.singular_values**2),
        rtol=1e-13,
    )
    assert_modes(report, K, SA, SE, POSTERIOR)
    assert report.innovation_chi_square is None
    description = report.describe()
    assert "temperature" in description and "water" in description
    assert "local" in description.lower()
    assert str(report) == description


def test_scalar_and_null_modes():
    for sensitivity, reduction in [(1.0, 0.5), (3.0, 0.9)]:
        report = information([[sensitivity]], [[4.0]], [[4.0]])
        np.testing.assert_allclose(report.singular_values, [sensitivity])
        np.testing.assert_allclose(report.mode_variance_reduction, [reduction])
        np.testing.assert_allclose(report.variance_reduction, [reduction])
        np.testing.assert_allclose(report.degrees_of_freedom, reduction)
        np.testing.assert_allclose(
            report.posterior_standard_deviation, [2 / np.sqrt(1 + sensitivity**2)]
        )

    zero = information(np.zeros_like(K), SA, SE)
    np.testing.assert_array_equal(zero.singular_values, [0, 0])
    np.testing.assert_allclose(zero.posterior_standard_deviation, np.sqrt(np.diag(SA)))
    np.testing.assert_allclose(zero.variance_reduction, [0, 0], atol=1e-14)
    assert zero.degrees_of_freedom == 0
    assert zero.information_bits == 0
    assert_modes(zero, np.zeros_like(K), SA, SE, SA)

    # Two state combinations are completely unobserved when n > m.
    jacobian = np.array([[1.0, 0, 0]])
    prior = np.diag([4.0, 9.0, 16.0])
    error = np.array([[4.0]])
    underdetermined = information(jacobian, prior, error)
    np.testing.assert_allclose(underdetermined.singular_values, [1, 0, 0])
    np.testing.assert_allclose(underdetermined.mode_variance_reduction, [0.5, 0, 0])
    np.testing.assert_allclose(
        underdetermined.posterior_standard_deviation, [np.sqrt(2), 3, 4]
    )
    assert_modes(underdetermined, jacobian, prior, error, np.diag([2.0, 9.0, 16.0]))


def test_units_and_permutations():
    base = information(K, SA, SE)
    state_scale = np.diag([1e-20, 1e3])
    measurement_scale = np.diag([1e-5, 1e2, 1e-1])
    scaled = information(
        measurement_scale @ K @ np.linalg.inv(state_scale),
        state_scale @ SA @ state_scale,
        measurement_scale @ SE @ measurement_scale,
    )
    np.testing.assert_allclose(scaled.singular_values, base.singular_values, rtol=1e-11)
    np.testing.assert_allclose(
        scaled.posterior_standard_deviation,
        np.diag(state_scale) * base.posterior_standard_deviation,
        rtol=1e-11,
    )
    np.testing.assert_allclose(scaled.variance_reduction, base.variance_reduction)

    state_order = [1, 0]
    measurement_order = [2, 0, 1]
    permuted = information(
        K[np.ix_(measurement_order, state_order)],
        SA[np.ix_(state_order, state_order)],
        SE[np.ix_(measurement_order, measurement_order)],
    )
    np.testing.assert_allclose(permuted.singular_values, base.singular_values)
    np.testing.assert_allclose(
        permuted.posterior_standard_deviation,
        base.posterior_standard_deviation[state_order],
    )

    # Relative information alone does not express absolute achievable accuracy.
    common = information(K, 100 * SA, 100 * SE)
    np.testing.assert_allclose(common.singular_values, base.singular_values)
    np.testing.assert_allclose(common.degrees_of_freedom, base.degrees_of_freedom)
    np.testing.assert_allclose(common.information_bits, base.information_bits)
    np.testing.assert_allclose(
        common.posterior_standard_deviation, 10 * base.posterior_standard_deviation
    )


def test_innovation():
    observation = np.array([2.0, -1.0, 1.5])
    prediction = K @ np.array([0.5, -0.25]) + np.array([0.25, -0.5, 1.0])
    residual = observation - prediction
    expected = residual @ np.linalg.solve(K @ SA @ K.T + SE, residual)
    report = information(
        K, SA, SE, measurement=observation, prior_prediction=prediction
    )
    np.testing.assert_allclose(report.innovation_chi_square, expected, rtol=1e-12)
    assert report.innovation_expected_mean == 3
    exact = information(K, SA, SE, measurement=prediction, prior_prediction=prediction)
    assert exact.innovation_chi_square == 0


def test_native_block_layout():
    # Block IDs do not determine coordinate order. These correlations are
    # stored in the upper block triangle but the lower physical triangle.
    prior = arts.CovarianceMatrix()
    prior.blocks = [
        arts.Block(arts.Range(1, 1), arts.Range(1, 1), (0, 0), arts.Matrix([[2.0]])),
        arts.Block(arts.Range(0, 1), arts.Range(0, 1), (2, 2), arts.Matrix([[4.0]])),
        arts.Block(arts.Range(1, 1), arts.Range(0, 1), (0, 2), arts.Matrix([[1.0]])),
    ]
    before = covariance_snapshot(prior)
    report = information(K, prior, covariance(SE))
    np.testing.assert_allclose(
        report.posterior_standard_deviation, np.sqrt(np.diag(POSTERIOR))
    )
    assert_modes(report, K, SA, SE, POSTERIOR)
    assert covariance_snapshot(prior) == before


def test_diagonal_and_memory_limit():
    diagonal = information(K, np.diag(SA), np.diag(SE))
    dense = information(K, np.diag(np.diag(SA)), np.diag(np.diag(SE)))
    np.testing.assert_allclose(diagonal.singular_values, dense.singular_values)
    np.testing.assert_allclose(
        diagonal.posterior_standard_deviation, dense.posterior_standard_deviation
    )

    # A large diagonal measurement covariance must not create m by m arrays.
    # The row sensitivities have a simple independent posterior reference.
    count = 10000
    jacobian = np.tile(np.eye(2), (count // 2, 1))
    report = information(
        jacobian, [4.0, 9.0], np.ones(count), max_dense_elements=500000
    )
    np.testing.assert_allclose(
        report.posterior_standard_deviation,
        1 / np.sqrt(np.array([0.25, 1 / 9]) + count / 2),
        rtol=1e-12,
    )
    assert report.measurement_modes.shape == (count, 2)
    rejects(lambda: information(K, SA, SE, max_dense_elements=1))

    # Mode information and posterior uncertainty remain representable even
    # when squaring the signal/noise ratio would overflow double precision.
    strong = information([[1e200]], [1.0], [1.0])
    assert np.isfinite(strong.information_bits)
    np.testing.assert_allclose(strong.information_bits, np.log2(1e200))
    np.testing.assert_allclose(
        strong.posterior_standard_deviation, [1e-200], atol=0, rtol=1e-12
    )
    # The final whitened Jacobian is finite even if applying the large prior
    # square root first overflows. Common covariance scaling leaves its
    # spectrum unchanged and only rescales physical uncertainties.
    strong_scaled = information([[1e200]], [1e300], [1e300])
    np.testing.assert_allclose(strong_scaled.singular_values, strong.singular_values)
    np.testing.assert_allclose(strong_scaled.information_bits, strong.information_bits)
    np.testing.assert_allclose(
        strong_scaled.posterior_standard_deviation, [1e-50], atol=0, rtol=1e-12
    )
    weak = information([[1e-10]], [1.0], [1.0])
    np.testing.assert_allclose(
        weak.information_bits, 0.5 * np.log1p(1e-20) / np.log(2), atol=0, rtol=1e-12
    )


def test_read_only_workspace():
    ws = pyarts.Workspace()
    ws.measurement_jac = K
    ws.model_state_covmat = covariance(SA)
    ws.measurement_vec_error_covmat = covariance(SE)
    ws.measurement_vec = [2.0, -1.0, 1.5]
    # Deliberately does not represent F(xa). The wrapper must not assume it does.
    ws.measurement_vec_fit = [99.0, 99.0, 99.0]
    before_jacobian = np.array(ws.measurement_jac, copy=True)
    before_prior = copy.deepcopy(ws.model_state_covmat)
    before_error = copy.deepcopy(ws.measurement_vec_error_covmat)
    before_prior_storage = covariance_snapshot(ws.model_state_covmat)
    before_error_storage = covariance_snapshot(ws.measurement_vec_error_covmat)
    report = information_from_workspace(ws)
    direct = information(K, SA, SE)
    np.testing.assert_allclose(report.singular_values, direct.singular_values)
    assert report.innovation_chi_square is None
    np.testing.assert_array_equal(ws.measurement_jac, before_jacobian)
    for before, current in (
        (before_prior, ws.model_state_covmat),
        (before_error, ws.measurement_vec_error_covmat),
    ):
        assert len(before.blocks) == len(current.blocks)
        for old, new in zip(before.blocks, current.blocks):
            np.testing.assert_array_equal(old.matrix, new.matrix)
    np.testing.assert_array_equal(ws.measurement_vec_fit, [99, 99, 99])
    assert covariance_snapshot(ws.model_state_covmat) == before_prior_storage
    assert covariance_snapshot(ws.measurement_vec_error_covmat) == before_error_storage
    report_with_residual = information_from_workspace(
        ws,
        prior_prediction=[0.25, 0.75, 1.25],
    )
    assert report_with_residual.innovation_chi_square is not None

    # The standard measurement helper supplies both sparse variance and inverse.
    # Read-only analysis must preserve that existing inverse as well.
    ws.measurement_sensor = [arts.SensorObsel() for _ in range(3)]
    ws.measurement_vec_error_covmatConstant(value=0.25)
    before_error_storage = covariance_snapshot(ws.measurement_vec_error_covmat)
    cached = information_from_workspace(ws)
    np.testing.assert_allclose(
        cached.singular_values, information(K, SA, 0.25 * np.eye(3)).singular_values
    )
    assert covariance_snapshot(ws.measurement_vec_error_covmat) == before_error_storage


def test_validation():
    # NumPy permits complex-to-real casting with a warning. This interface
    # must reject it instead of analyzing silently changed physical inputs.
    rejects(lambda: information(K.astype(complex) + 1j, SA, SE))
    rejects(lambda: information(K, SA.astype(complex) + 1j, SE))
    rejects(lambda: information(K, SA, SE.astype(complex) + 1j))
    rejects(
        lambda: information(
            K,
            SA,
            SE,
            measurement=np.array([1 + 1j, 2, 3]),
            prior_prediction=[1, 2, 3],
        )
    )
    rejects(
        lambda: information(
            K,
            SA,
            SE,
            measurement=[1, 2, 3],
            prior_prediction=np.array([1, 2 + 1j, 3]),
        )
    )
    for invalid in (
        np.array([[-1.0, 0], [0, -2]]),
        np.array([[1.0, 2], [2, 1]]),
        np.array([[1.0, 1], [1, 1]]),
        np.array([[1.0, 0.1], [0.2, 1]]),
        np.array([[1.0, np.nan], [np.nan, 1]]),
        np.array([[1.0, np.inf], [np.inf, 1]]),
    ):
        rejects(lambda: information(K, invalid, SE))
        rejects(lambda: information(K, covariance(invalid), SE))
    for invalid in ([0.0, 1], [-1.0, 1], [np.nan, 1], [np.inf, 1]):
        rejects(lambda: information(K, invalid, SE))
    for jacobian, prior, error in (
        (K, np.eye(3), SE),
        (K, SA, np.eye(2)),
        (np.ones(2), SA, SE),
        (np.empty((0, 2)), SA, np.empty((0, 0))),
        (np.ones((3, 0)), np.empty((0, 0)), SE),
        (np.full_like(K, np.nan), SA, SE),
        (np.full_like(K, np.inf), SA, SE),
    ):
        rejects(lambda: information(jacobian, prior, error))
    rejects(lambda: information(K, SA, SE, state_labels=["only one"]))
    rejects(lambda: information(K, SA, SE, measurement=[1, 2, 3]))
    rejects(lambda: information(K, SA, SE, prior_prediction=[1, 2, 3]))
    rejects(
        lambda: information(K, SA, SE, measurement=[1, 2], prior_prediction=[1, 2, 3])
    )
    rejects(
        lambda: information(
            K, SA, SE, measurement=[1, np.nan, 3], prior_prediction=[1, 2, 3]
        )
    )
    for tolerance in (0, -1, np.nan, np.inf):
        rejects(lambda: information(K, SA, SE, relative_tolerance=tolerance))
    for limit in (0, -1, 1.5, np.nan):
        rejects(lambda: information(K, SA, SE, max_dense_elements=limit))


for test in (
    test_correlated_reference,
    test_scalar_and_null_modes,
    test_units_and_permutations,
    test_innovation,
    test_native_block_layout,
    test_diagonal_and_memory_limit,
    test_read_only_workspace,
    test_validation,
):
    test()
