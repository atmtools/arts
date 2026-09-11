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


def check_reductions():
    report = information(K, SA, SE)
    complete = report.reduction(rank=2)
    np.testing.assert_allclose(complete.posterior_covariance(), POSTERIOR, rtol=1e-12)
    assert complete.discarded_degrees_of_freedom == 0
    assert complete.discarded_information_bits == 0
    assert not complete.model_state_basis_mat.flags.writeable
    assert not complete.measurement_basis_mat.flags.writeable
    for rank in (1, 2):
        reduction = report.reduction(rank=rank)
        b, c = reduction.model_state_basis_mat, reduction.measurement_basis_mat
        np.testing.assert_allclose(
            b.T @ np.linalg.solve(SA, b), np.eye(rank), atol=1e-12)
        np.testing.assert_allclose(c @ SE @ c.T, np.eye(rank), atol=1e-12)
        np.testing.assert_allclose(
            c @ K @ b, np.diag(report.singular_values[:rank]), atol=1e-12)
    reduced = report.reduction(rank=1)
    # Independent posterior for the truncated forward matrix, with the full prior.
    b = reduced.model_state_basis_mat
    projected_j = K @ b @ b.T @ np.linalg.inv(SA)
    reference = np.linalg.inv(np.linalg.inv(SA) + projected_j.T @
                              np.linalg.solve(SE, projected_j))
    np.testing.assert_allclose(reduced.posterior_covariance(), reference, rtol=1e-12)
    discarded = report.state_modes[:, 1:]
    prior_coordinates = np.linalg.solve(SA, discarded)
    np.testing.assert_allclose(
        prior_coordinates.T @ reduced.posterior_covariance() @ prior_coordinates, [[1.]], atol=1e-12
    )
    assert "Retaining 1 of 2" in str(reduced)

    # Equal weak modes: DOFS=1, but preserving 90% of it requires 90 modes.
    broad = information(np.eye(100) / np.sqrt(99), np.eye(100), np.eye(100))
    np.testing.assert_allclose(broad.degrees_of_freedom, 1)
    assert broad.reduction(max_lost_dofs=0.100001).rank == 90
    spectrum = information(np.diag([3., 1., .1, 0.]), np.eye(4), np.eye(4))
    assert spectrum.reduction(max_lost_dofs=.01).rank == 2
    assert spectrum.reduction(max_lost_information_bits=.01).rank == 2
    assert spectrum.reduction(
        max_lost_dofs=.6, max_lost_information_bits=.001).rank == 3
    assert spectrum.reduction(max_lost_dofs=0, max_lost_information_bits=0).rank == 3
    null = information(np.zeros((2, 3)), np.eye(
        3), np.eye(2)).reduction(max_lost_dofs=0)
    assert null.rank == 1
    complete_null = information(np.zeros((2, 3)), np.eye(3),
                                np.eye(2)).reduction(rank=3)
    assert complete_null.model_state_basis_mat.shape == (3, 3)
    assert complete_null.measurement_basis_mat.shape == (2, 2)
    np.testing.assert_array_equal(null.posterior_covariance(), np.eye(3))
    for kwargs in ({}, {"rank": 0}, {"rank": 3}, {"rank": 1.5}, {"rank": True},
                   {"rank": 1, "max_lost_dofs": 0}, {"max_lost_dofs": -1},
                   {"max_lost_information_bits": np.nan}, {"max_lost_dofs": np.inf}):
        rejects(lambda: report.reduction(**kwargs))


check_reductions()


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
    reduction = report.reduction(rank=n)
    b, c = reduction.model_state_basis_mat, reduction.measurement_basis_mat
    q = min(m, n)
    expected = np.zeros((q, n))
    np.fill_diagonal(expected, singular[:q])
    np.testing.assert_allclose(c @ jacobian @ b, expected, atol=1e-10)
    np.testing.assert_allclose(c @ error @ c.T, np.eye(q), atol=1e-12)
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


def test_workspace_bases():
    def check(jacobian, prior, noise, rank, *, prior_input=None, noise_input=None):
        ws = pyarts.Workspace()
        ws.measurement_jac = jacobian
        ws.model_state_covmat = covariance(
            prior) if prior_input is None else prior_input
        ws.measurement_vec_error_covmat = covariance(
            noise) if noise_input is None else noise_input
        before_prior = covariance_snapshot(ws.model_state_covmat)
        before_noise = covariance_snapshot(ws.measurement_vec_error_covmat)
        ws.ReducedOEMBasisCalc()
        m, n = jacobian.shape
        full_b = np.array(ws.model_state_basis_mat)
        full_c = np.array(ws.measurement_basis_mat)
        singular = np.array(ws.oem_basis_singular_values)
        assert full_b.shape == (n, n)
        assert full_c.shape == (m, m)
        assert singular.shape == (min(m, n),)
        np.testing.assert_allclose(full_b @ full_b.T, prior, atol=1e-12)
        inverse_c = np.linalg.solve(full_c, np.eye(m))
        np.testing.assert_allclose(inverse_c @ inverse_c.T, noise, atol=1e-12)
        np.testing.assert_allclose(full_c @ noise @ full_c.T, np.eye(m), atol=1e-12)
        sigma = np.zeros((m, n))
        np.fill_diagonal(sigma, singular)
        np.testing.assert_allclose(full_c @ jacobian @ full_b, sigma, atol=1e-12)
        ws.ReducedOEMBasisReduce(rank=rank)
        b, c = np.array(ws.model_state_basis_mat), np.array(ws.measurement_basis_mat)
        q = min(rank, m)
        assert b.shape == (n, rank)
        assert c.shape == (q, m)
        np.testing.assert_allclose(
            b.T @ np.linalg.solve(prior, b), np.eye(rank), atol=1e-12)
        np.testing.assert_allclose(c @ noise @ c.T, np.eye(q), atol=1e-12)
        report = information(jacobian, prior, noise)
        np.testing.assert_allclose(
            singular, report.singular_values[:min(m, n)], atol=1e-12)
        expected = report.reduction(rank=rank)
        np.testing.assert_allclose(ws.oem_basis_lost_dofs,
                                   expected.discarded_degrees_of_freedom, atol=1e-12)
        np.testing.assert_allclose(
            ws.oem_basis_lost_information_bits, expected.discarded_information_bits, atol=1e-12)
        np.testing.assert_allclose(np.linalg.svd(c @ jacobian @ b, compute_uv=False),
                                   report.singular_values[:q], atol=1e-12)
        if rank == n:
            np.testing.assert_allclose(b @ b.T, prior, atol=1e-12)
            jr = c @ jacobian @ b
            gain = b @ np.linalg.solve(np.eye(rank) + jr.T @ jr, jr.T) @ c
            reference = np.linalg.solve(np.linalg.inv(prior) + jacobian.T @ np.linalg.solve(noise, jacobian),
                                        np.linalg.solve(noise, jacobian).T)
            np.testing.assert_allclose(gain, reference, atol=1e-12)
        elif np.count_nonzero(report.singular_values) == n:
            # Singular-vector signs are arbitrary, so compare subspaces.
            np.testing.assert_allclose(
                b @ b.T, expected.model_state_basis_mat @ expected.model_state_basis_mat.T, atol=1e-12)
        assert covariance_snapshot(ws.model_state_covmat) == before_prior
        assert covariance_snapshot(ws.measurement_vec_error_covmat) == before_noise
        np.testing.assert_array_equal(ws.measurement_jac, jacobian)
        np.testing.assert_array_equal(ws.model_state_basis_mat, full_b[:, :rank])
        np.testing.assert_array_equal(ws.measurement_basis_mat, full_c[:q, :])
        np.testing.assert_array_equal(ws.oem_basis_singular_values, singular)
        return ws

    for rank in (1, 2):
        check(K, SA, SE, rank)
    for storage in (arts.Matrix, arts.Sparse):
        # A connected prior with block IDs in reversed coordinate order.
        prior = arts.CovarianceMatrix()
        prior.blocks = [
            arts.Block(arts.Range(1, 1), arts.Range(1, 1), (0, 0), storage([[2.]])),
            arts.Block(arts.Range(0, 1), arts.Range(0, 1), (2, 2), storage([[4.]])),
            arts.Block(arts.Range(1, 1), arts.Range(0, 1), (0, 2), storage([[1.]])),
        ]
        noise = arts.CovarianceMatrix()
        noise.blocks = [arts.Block(arts.Range(
            0, 3), arts.Range(0, 3), (0, 0), storage(SE))]
        check(K, SA, SE, 1, prior_input=prior, noise_input=noise)
    check(np.array([[1., 0, 0]]), np.diag([4., 9., 16.]), np.array([[4.]]), 3)
    check(np.zeros((3, 2)), SA, SE, 1)
    check(np.eye(2), np.eye(2), np.eye(2), 2)  # repeated singular values
    check(np.array([[2.]]), np.array([[4.]]), np.array([[9.]]), 1)

    mixed = arts.CovarianceMatrix()
    mixed.blocks = [
        arts.Block(arts.Range(0, 1), arts.Range(0, 1), (0, 0), arts.Sparse([[1.]])),
        arts.Block(arts.Range(1, 2), arts.Range(1, 2),
                   (1, 1), arts.Matrix([[2., .3], [.3, .5]])),
    ]
    check(K, SA, np.array([[1., 0, 0], [0, 2., .3], [0, .3, .5]]), 2, noise_input=mixed)
    cached = pyarts.Workspace()
    cached.measurement_sensor = [arts.SensorObsel() for _ in range(3)]
    cached.measurement_vec_error_covmatConstant(value=.25)
    check(K, SA, .25 * np.eye(3), 2, noise_input=cached.measurement_vec_error_covmat)

    # Diagonal noise uses a diagonal factor; the saved full basis is square.
    count = 100
    ws = pyarts.Workspace()
    ws.measurement_jac = np.tile(np.eye(2), (count // 2, 1))
    ws.model_state_covmat = covariance(SA)
    ws.measurement_vec_error_covmat = arts.CovarianceMatrix()
    from scipy import sparse
    ws.measurement_vec_error_covmat.blocks = [arts.Block(
        arts.Range(0, count), arts.Range(0, count), (0, 0),
        arts.Sparse(sparse.eye(count, format="csr")))]
    ws.ReducedOEMBasisCalc()
    assert np.asarray(ws.measurement_basis_mat).shape == (count, count)
    ws.ReducedOEMBasisReduce(rank=1)
    c = np.asarray(ws.measurement_basis_mat)
    np.testing.assert_allclose(c @ c.T, [[1.]], atol=1e-12)

    # Failure leaves the previously generated bases intact.
    ws = check(K, SA, SE, 1)
    b, c = np.array(ws.model_state_basis_mat), np.array(ws.measurement_basis_mat)
    lost = (float(ws.oem_basis_lost_dofs), float(ws.oem_basis_lost_information_bits))
    for kwargs in ({"rank": 0}, {"rank": 3}, {"rank": -2},
                   {"rank": 1, "max_lost_dofs": 0},
                   {"rank": 1, "max_lost_information_bits": 0},
                   {"max_lost_dofs": -.1}, {"max_lost_dofs": np.nan},
                   {"max_lost_information_bits": -2}, {"max_lost_information_bits": np.inf}):
        rejects(lambda: ws.ReducedOEMBasisReduce(**kwargs))
        np.testing.assert_array_equal(ws.model_state_basis_mat, b)
        np.testing.assert_array_equal(ws.measurement_basis_mat, c)
        assert (float(ws.oem_basis_lost_dofs), float(
            ws.oem_basis_lost_information_bits)) == lost
    singular = np.array(ws.oem_basis_singular_values)
    for jacobian, prior, noise in ((K, np.eye(3), SE), (K, SA, np.eye(2)),
                                   (K, [[1., 2.], [2., 1.]], SE),
                                   (np.full((3, 2), np.nan), SA, SE)):
        ws.measurement_jac = jacobian
        ws.model_state_covmat = covariance(prior)
        ws.measurement_vec_error_covmat = covariance(noise)
        rejects(lambda: ws.ReducedOEMBasisCalc())
        np.testing.assert_array_equal(ws.model_state_basis_mat, b)
        np.testing.assert_array_equal(ws.measurement_basis_mat, c)
        np.testing.assert_array_equal(ws.oem_basis_singular_values, singular)


test_workspace_bases()


def test_workspace_mode_selection():
    def select(jacobian, **options):
        m, n = jacobian.shape
        ws = pyarts.Workspace()
        ws.measurement_jac = jacobian
        ws.model_state_covmat = covariance(np.eye(n))
        ws.measurement_vec_error_covmat = covariance(np.eye(m))
        ws.ReducedOEMBasisCalc()
        ws.ReducedOEMBasisReduce(**options)
        b, c = np.asarray(ws.model_state_basis_mat), np.asarray(
            ws.measurement_basis_mat)
        rank = b.shape[1]
        np.testing.assert_allclose(b.T @ b, np.eye(rank), atol=1e-12)
        np.testing.assert_allclose(c @ c.T, np.eye(min(rank, m)), atol=1e-12)
        # Match the established report's loss-budget selection.
        limits = {key: value for key, value in options.items()
                  if key != "rank" and value != -1}
        if not limits:
            limits = {"max_lost_information_bits": 0}
        expected = information(jacobian, np.eye(n), np.eye(m)).reduction(**limits)
        assert rank == expected.rank, (options, rank, expected.rank)
        np.testing.assert_allclose(ws.oem_basis_lost_dofs,
                                   expected.discarded_degrees_of_freedom, atol=1e-12)
        np.testing.assert_allclose(
            ws.oem_basis_lost_information_bits, expected.discarded_information_bits, atol=1e-12)
        return rank

    spectrum = np.diag([3., 1., .1, 0.])
    assert select(spectrum) == 3
    assert select(spectrum, max_lost_dofs=.01) == 2
    assert select(spectrum, max_lost_information_bits=.01) == 2
    assert select(spectrum, max_lost_dofs=.6, max_lost_information_bits=.001) == 3
    assert select(spectrum, max_lost_dofs=0, max_lost_information_bits=0) == 3
    assert select(spectrum, rank=-1, max_lost_dofs=-
                  1, max_lost_information_bits=-1) == 3
    assert select(np.zeros((2, 3))) == 1  # ReducedOEM still requires one coefficient.
    assert select(np.array([[1., 0., 0.]])) == 1

    # Budgets apply to the SUM of discarded contributions, not each mode.
    broad = np.eye(100) / np.sqrt(99)
    assert select(broad, max_lost_dofs=.100001) == 90

    # Never judge importance only relative to the strongest singular value.
    strong = np.diag([1e100, 1e-6, 0.])
    assert select(strong) == 2
    assert select(strong, max_lost_information_bits=1e-11) == 1

    # These redundant directions contain no individual zero entries. A small
    # budget also removes roundoff residuals in mathematically null modes.
    assert select(np.ones((3, 4)), max_lost_information_bits=1e-12) == 1

    # Recalculation restores the state null space removed by selection.
    from scipy import sparse
    ws = pyarts.Workspace()
    n = 200
    ws.measurement_jac = np.r_[1., np.zeros(n - 1)].reshape(1, n)
    ws.model_state_covmat = arts.CovarianceMatrix()
    ws.model_state_covmat.blocks = [arts.Block(
        arts.Range(0, n), arts.Range(0, n), (0, 0),
        arts.Sparse(sparse.eye(n, format="csr")))]
    ws.measurement_vec_error_covmat = covariance([[1.]])
    ws.ReducedOEMBasisCalc()
    assert np.asarray(ws.model_state_basis_mat).shape == (n, n)
    ws.ReducedOEMBasisReduce()
    b = np.asarray(ws.model_state_basis_mat)
    c = np.asarray(ws.measurement_basis_mat)
    assert b.shape == (n, 1) and c.shape == (1, 1)
    np.testing.assert_allclose(
        b @ c, np.r_[1., np.zeros(n - 1)].reshape(n, 1), atol=1e-12)
    rejects(lambda: ws.ReducedOEMBasisReduce(rank=n))
    ws.ReducedOEMBasisCalc()
    ws.ReducedOEMBasisReduce(rank=n)
    full_b = np.asarray(ws.model_state_basis_mat)
    np.testing.assert_allclose(full_b @ full_b.T, np.eye(n), atol=1e-12)


test_workspace_mode_selection()


def test_workspace_basis_reselection():
    for jacobian, prior, noise in ((K, SA, SE), (K.T, SE, SA)):
        ws = pyarts.Workspace()
        ws.measurement_jac = jacobian
        ws.model_state_covmat = covariance(prior)
        ws.measurement_vec_error_covmat = covariance(noise)
        ws.ReducedOEMBasisCalc()
        b = np.array(ws.model_state_basis_mat)
        c = np.array(ws.measurement_basis_mat)
        singular = np.array(ws.oem_basis_singular_values)
        m, n = jacobian.shape
        report = information(jacobian, prior, noise)

        # Selection uses only the saved decomposition, even if the current
        # Jacobian/covariances no longer form a valid setup for another SVD.
        ws.measurement_jac = [[np.nan]]
        ws.model_state_covmat = covariance([[1.]])
        ws.measurement_vec_error_covmat = covariance([[1.]])
        for rank in (n, 1, 1):
            ws.ReducedOEMBasisReduce(rank=rank)
            np.testing.assert_array_equal(ws.model_state_basis_mat, b[:, :rank])
            np.testing.assert_array_equal(ws.measurement_basis_mat, c[:min(rank, m), :])
            expected = report.reduction(rank=rank)
            np.testing.assert_allclose(
                ws.oem_basis_lost_dofs, expected.discarded_degrees_of_freedom, atol=1e-12)
            np.testing.assert_allclose(
                ws.oem_basis_lost_information_bits, expected.discarded_information_bits, atol=1e-12)
            np.testing.assert_array_equal(ws.oem_basis_singular_values, singular)
        rejects(lambda: ws.ReducedOEMBasisReduce(rank=n))
        rejects(lambda: ws.ReducedOEMBasisReduce(max_lost_information_bits=0))
        np.testing.assert_array_equal(ws.model_state_basis_mat, b[:, :1])
        np.testing.assert_array_equal(ws.measurement_basis_mat, c[:1, :])
        # Explicit restoration permits selecting more modes again.
        ws.model_state_basis_mat = b
        ws.measurement_basis_mat = c
        ws.ReducedOEMBasisReduce(rank=n)

        # Malformed saved inputs are rejected before replacing active outputs.
        ws.ReducedOEMBasisReduce(rank=1)
        old_b = np.array(ws.model_state_basis_mat)
        old_c = np.array(ws.measurement_basis_mat)
        for bad_spectrum in ([], [1.], [1., 2.], [-1., -2.], [np.inf, 0.], [np.nan, 0.]):
            ws.oem_basis_singular_values = bad_spectrum
            rejects(lambda: ws.ReducedOEMBasisReduce())
            np.testing.assert_array_equal(ws.model_state_basis_mat, old_b)
            np.testing.assert_array_equal(ws.measurement_basis_mat, old_c)
        ws.oem_basis_singular_values = singular
        for name, saved in (("model_state_basis_mat", old_b), ("measurement_basis_mat", old_c)):
            setattr(ws, name, np.zeros((0, 0)))
            rejects(lambda: ws.ReducedOEMBasisReduce())
            assert np.asarray(getattr(ws, name)).shape == (0, 0)
            setattr(ws, name, saved)


test_workspace_basis_reselection()


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
    reduction = report.reduction(rank=1)
    assert reduction.measurement_basis_mat.shape == (1, count)
    np.testing.assert_allclose(reduction.measurement_basis_mat @
                               reduction.measurement_basis_mat.T, [[1.]], atol=1e-12)
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
    ws.measurement_noise_scales = arts.Vector()
    ws.measurement_vec_error_covmatNormalization(
        normalization=ws.measurement_noise_scales
    )
    np.testing.assert_allclose(ws.measurement_noise_scales, np.sqrt(np.diag(SE)))
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


def test_workspace_state_labels():
    ws = pyarts.Workspace()
    ws.atm_field = arts.AtmField()
    ws.surf_field = arts.SurfaceField()
    ws.surf_field.ellipsoid = [1, 1]
    ws.abs_bands = arts.AbsorptionBands()
    ws.measurement_sensor = arts.ArrayOfSensorObsel()
    ws.subsurf_field = arts.SubsurfaceField()
    for key, count, value in (("t", 2, 280.), ("H2O", 3, .01)):
        ws.atm_field[key] = arts.GriddedField3(
            data=np.full((count, 1, 1), value),
            grid_names=["Altitude", "Latitude", "Longitude"],
            grids=[np.arange(count, dtype=float), [0], [0]],
        )
    ws.RetrievalInit()
    ws.RetrievalAddSpeciesVMR(species="H2O", matrix=np.eye(3))
    ws.RetrievalAddTemperature(matrix=np.eye(2))
    ws.RetrievalFinalizeDiagonal()
    ws.measurement_jac = np.arange(15, dtype=float).reshape(3, 5) / 10
    ws.measurement_vec_error_covmat = covariance(np.eye(3))
    ws.measurement_vec = [0, 0, 0]
    report = information_from_workspace(ws)
    for target in ws.jac_targets.atm:
        name = f"atm.{target.type}"
        start, count = target.x_start, target.x_size
        assert (name, start, count) in report.state_blocks
        assert report.state_labels[start:start+count] == tuple(
            f"{name}[{i}]" for i in range(count))
        assert f"x[{start}:{start+count}]: {name}" in report.describe()
    custom = information_from_workspace(ws, state_labels=list("abcde"))
    assert custom.state_labels == tuple("abcde")
    assert custom.state_blocks == report.state_blocks
    np.testing.assert_array_equal(custom.state_modes, report.state_modes)
    # Stale ranges must not silently attach a field name to the wrong column.
    ws.measurement_jac = np.ones((3, 2))
    rejects(lambda: information_from_workspace(ws))


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
    test_workspace_state_labels,
    test_validation,
):
    test()
