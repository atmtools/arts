"""Covariance checks and local information analysis before or after OEM.

The information spectrum follows Nesser et al. (2021), Sect. 2.2,
https://doi.org/10.5194/amt-14-5521-2021. These functions neither execute an
agenda nor change a workspace, its covariance matrices, or its settings.
"""

from dataclasses import dataclass, field, replace
import operator

import numpy as np
from scipy import linalg, sparse

from . import arts

__all__ = ["InformationReport", "ReductionReport", "information", "information_from_workspace"]


def _real_array(value, name):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must contain real values")
    return np.asarray(value, dtype=float)


def _readonly(value):
    value = np.array(value, dtype=float, copy=True)
    value.flags.writeable = False
    return value


def _mode_information_bits(singular_values):
    # Avoid overflow for strong modes and cancellation for weak modes.
    weak = singular_values <= 1
    result = np.empty_like(singular_values)
    result[weak] = 0.5 * np.log1p(singular_values[weak] ** 2)
    strong = singular_values[~weak]
    result[~weak] = np.log(strong) + 0.5 * np.log1p((1 / strong) ** 2)
    return result / np.log(2.0)


@dataclass(frozen=True)
class ReductionReport:
    """Fixed state/measurement reductions and their local linear information loss.

    Pass both reduction matrices to ``Workspace.ReducedOEM``. Losses describe the Jacobian
    used for the information report, not a bound on nonlinear retrieval error.
    The underlying state modes are shared with that read-only report.
    """

    model_state_red_mat: np.ndarray = field(repr=False)
    measurement_red_mat: np.ndarray = field(repr=False)
    retained_degrees_of_freedom: float
    discarded_degrees_of_freedom: float
    retained_information_bits: float
    discarded_information_bits: float
    _state_modes: np.ndarray = field(repr=False)
    _posterior_factors: np.ndarray = field(repr=False)

    @property
    def rank(self):
        """Number of retained columns, and hence reduced state variables."""
        return self.model_state_red_mat.shape[1]

    def posterior_covariance(self):
        """Return the local full-state covariance, retaining discarded priors.

        This allocates an n by n matrix. It is exact for a linear Gaussian
        model when all nonzero modes are retained. Otherwise it is the
        posterior of the truncated linear model. It is not recomputed at the
        state returned by a subsequent ReducedOEM call.
        """
        scaled = self._state_modes / self._posterior_factors
        return _readonly(scaled @ scaled.T)

    def __str__(self):
        return (
            f"Retaining {self.rank} of {self.model_state_red_mat.shape[0]} state modes; "
            f"discarding {self.discarded_degrees_of_freedom:.6g} DOFS and "
            f"{self.discarded_information_bits:.6g} bits (local linear analysis)"
        )


@dataclass(frozen=True)
class InformationReport:
    """Snapshot of information implied by a Jacobian and assumed covariances.

    All arrays are independent, read-only snapshots. ``singular_values`` has
    one entry per state mode, including unobserved null modes. Columns of
    ``state_modes`` are perturbations in the supplied state coordinates:
    ``state_modes @ state_modes.T`` equals the prior covariance. Mode signs
    and bases within repeated singular values are not unique.

    ``measurement_modes`` contains the left singular vectors in whitened
    measurement coordinates, not physical measurement units. Its number of
    columns is the smaller of the measurement and state dimensions.

    Standard deviations and marginal ``variance_reduction`` refer to the
    original state coordinates. ``mode_variance_reduction`` refers to the
    independent prior-normalized modes. Correlated priors can improve a
    parameter indirectly through information about other parameters.

    ``state_labels`` identify individual state coordinates. Workspace reports
    also include ``state_blocks`` tuples of (field name, x_start, x_size).
    These indices describe flattened retrieval coordinates, not physical grids.

    Innovation statistics are present only when both measurements and an
    explicit prior prediction were supplied. Their interpretation assumes
    independent Gaussian prior and observation errors and a valid linear
    approximation around the prior.
    """

    singular_values: np.ndarray
    mode_variance_reduction: np.ndarray
    degrees_of_freedom: float
    information_bits: float
    prior_standard_deviation: np.ndarray
    posterior_standard_deviation: np.ndarray
    variance_reduction: np.ndarray
    state_modes: np.ndarray
    measurement_modes: np.ndarray
    state_labels: tuple[str, ...]
    prior_correlation_condition: float
    measurement_correlation_condition: float
    innovation_chi_square: float | None
    innovation_expected_mean: int | None
    _noise: "_CovarianceFactor" = field(repr=False)
    # (field label, x_start, x_size), in state-vector order.
    state_blocks: tuple[tuple[str, int, int], ...] = ()

    def reduction(self, rank=None, *, max_lost_dofs=None, max_lost_information_bits=None):
        """Select leading modes for ReducedOEM by rank or absolute loss limits.

        Supply either an explicit integer rank in 1..n, or one or both finite,
        nonnegative loss limits. Limits select the smallest rank satisfying
        all supplied limits, retaining at least one column even for zero
        information. Limits are absolute DOFS/bits, not fractions. Nothing is
        rounded from the total DOFS; a broad weak spectrum can require many
        more modes than that total suggests. Zero limits discard only modes
        with zero information in the computed spectrum.

        The returned matrices are B = L_a V_r and C = U_q.T L_e^-1, where
        q = min(rank, m). The reduced prior and noise covariances are identity.
        C contains weighted combinations, not a selection of physical channels.
        """
        limits = (max_lost_dofs, max_lost_information_bits)
        supplied = [limit is not None for limit in limits]
        if (rank is None and not any(supplied)) or (rank is not None and any(supplied)):
            raise ValueError("Supply either rank or information-loss limits")
        n = len(self.singular_values)
        bits = _mode_information_bits(self.singular_values)
        dofs = self.mode_variance_reduction
        # Reverse sums preserve small discarded tails beside large leading modes.
        tails = [np.r_[np.cumsum(values[::-1])[::-1], 0.0] for values in (dofs, bits)]
        if rank is not None:
            if isinstance(rank, (bool, np.bool_)):
                raise TypeError("rank must be an integer, not a boolean")
            rank = operator.index(rank)
            if not 1 <= rank <= n:
                raise ValueError(f"rank must be between 1 and {n}")
        else:
            eligible = np.ones(n + 1, dtype=bool)
            eligible[0] = False
            for limit, tail in zip(limits, tails):
                if limit is None:
                    continue
                limit = float(limit)
                if not np.isfinite(limit) or limit < 0:
                    raise ValueError("Information-loss limits must be finite and nonnegative")
                eligible &= tail <= limit
            rank = int(np.flatnonzero(eligible)[0])
        factors = np.ones(n)
        factors[:rank] = np.hypot(1.0, self.singular_values[:rank])
        return ReductionReport(
            model_state_red_mat=_readonly(self.state_modes[:, :rank]),
            measurement_red_mat=_readonly(self._noise.solve_left(
                self.measurement_modes[:, :rank], transpose=True
            ).T),
            retained_degrees_of_freedom=float(np.sum(dofs[:rank])),
            discarded_degrees_of_freedom=float(tails[0][rank]),
            retained_information_bits=float(np.sum(bits[:rank])),
            discarded_information_bits=float(tails[1][rank]),
            _state_modes=self.state_modes,
            _posterior_factors=_readonly(factors),
        )

    def describe(self, max_states=10):
        """Explain the report, showing at most ``max_states`` state entries."""
        max_states = operator.index(max_states)
        if max_states < 0:
            raise ValueError("max_states must be nonnegative")
        n = len(self.singular_values)
        m = self.measurement_modes.shape[0]
        lines = [
            f"Local OEM information: {m} measurements, {n} state variables",
            f"Degrees of freedom for signal: {self.degrees_of_freedom:.6g} of {n}",
            f"Expected information gain: {self.information_bits:.6g} bits",
            "Modes with signal/noise >= 1: "
            f"{np.count_nonzero(self.singular_values >= 1)} of {n}",
            "Correlation condition numbers (independent of coordinate units): "
            f"prior {self.prior_correlation_condition:.6g}, "
            f"measurement {self.measurement_correlation_condition:.6g}",
        ]
        if self.innovation_chi_square is not None:
            lines.append(
                f"Prior innovation chi-square: {self.innovation_chi_square:.6g}; "
                f"expected mean {self.innovation_expected_mean} "
                "under the linear Gaussian model"
            )
        if self.state_blocks:
            lines.extend(["", "State-vector fields (Python slices):"])
            for label, start, size in self.state_blocks:
                lines.append(f"  x[{start}:{start + size}]: {label} ({size} entries)")
        if max_states:
            lines.extend(["", "Linear analysis", "State: prior SD -> posterior SD; variance removed"])
            for label, prior, posterior, reduction in zip(
                self.state_labels[:max_states],
                self.prior_standard_deviation[:max_states],
                self.posterior_standard_deviation[:max_states],
                self.variance_reduction[:max_states],
            ):
                lines.append(
                    f"  {label}: {prior:.6g} -> {posterior:.6g}; {100 * reduction:.3g}%"
                )
            if n > max_states:
                lines.append(
                    f"  ... {n - max_states} further states are available in the arrays"
                )

        return "\n".join(lines)

    def __str__(self):
        return self.describe()

    def __repr__(self):
        return self.describe()

    def plot(self):
        """Return ``(figure, axes)`` with mode and state uncertainty diagnostics.

        Both panels are dimensionless, so quantities with different physical
        units can be compared. This does not call ``matplotlib.pyplot.show``.
        """
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
        n = len(self.singular_values)
        axes[0].plot(np.arange(1, n + 1), self.mode_variance_reduction, "o-")
        axes[0].axhline(0.5, color="0.6", linestyle="--", label="Signal/noise = 1")
        axes[0].set(
            xlabel="State mode (strongest first)",
            ylabel="Fraction of mode variance removed",
            ylim=(0, 1.05),
        )
        axes[0].legend()
        axes[1].plot(
            np.arange(n),
            self.posterior_standard_deviation / self.prior_standard_deviation,
            "o",
        )
        axes[1].axhline(1, color="0.6", linestyle="--")
        axes[1].set(
            xlabel="State index", ylabel="Posterior SD / prior SD", ylim=(0, 1.05)
        )
        if n <= 15:
            axes[1].set_xticks(np.arange(n), self.state_labels, rotation=45, ha="right")
        fig.suptitle(
            f"Local information: {self.degrees_of_freedom:.3g} degrees of freedom"
        )
        return fig, axes


class _CovarianceFactor:
    """Square-root factors on disjoint coordinate subsets, with diagonal fast paths."""

    def __init__(self, components, size, standard_deviation, condition):
        self.components = components
        self.size = size
        self.standard_deviation = standard_deviation
        self.condition = condition

    def multiply_left(self, values):
        result = np.empty_like(values)
        for indices, factor in self.components:
            if factor.ndim == 1:
                result[indices] = factor[:, None] * values[indices]
            else:
                result[indices] = factor @ values[indices]
        return result

    def multiply_right(self, values):
        result = np.empty_like(values)
        for indices, factor in self.components:
            if factor.ndim == 1:
                result[:, indices] = values[:, indices] * factor
            else:
                result[:, indices] = values[:, indices] @ factor
        return result

    def solve_left(self, values, *, transpose=False):
        result = np.empty_like(values)
        for indices, factor in self.components:
            if factor.ndim == 1:
                divisor = factor if values.ndim == 1 else factor[:, None]
                result[indices] = values[indices] / divisor
            else:
                result[indices] = linalg.solve_triangular(
                    factor, values[indices], lower=True, check_finite=False,
                    trans="T" if transpose else "N",
                )
        return result


def _factor_covariance(value, size, name, tolerance, max_dense_elements):
    if not isinstance(value, arts.CovarianceMatrix):
        array = _real_array(value, name)
        if array.shape == (size,):
            payload = arts.Sparse(sparse.diags(array, format="csr"))
        elif array.shape == (size, size):
            if size * size > max_dense_elements:
                raise ValueError(f"{name}: dense covariance exceeds max_dense_elements")
            payload = arts.Matrix(array)
        else:
            raise ValueError(
                f"{name} must have shape ({size}, {size}) or ({size},) variances"
            )
        value = arts.CovarianceMatrix()
        value.blocks = [
            arts.Block(arts.Range(0, size), arts.Range(0, size), (0, 0), payload)
        ]
    try:
        value.validate(
            expected_size=size,
            relative_tolerance=tolerance,
            max_dense_elements=max_dense_elements,
        )
    except (RuntimeError, ValueError) as error:
        raise ValueError(f"{name}: {error}") from error

    blocks = value.blocks
    diagonals = {}
    neighbors = {}
    data = {}
    for block in blocks:
        i, j = block.indices
        payload = block.matrix
        data[i, j] = (
            payload.tocsr()
            if isinstance(payload, arts.Sparse)
            else np.array(payload, copy=True)
        )
        neighbors.setdefault(i, set()).add(j)
        neighbors.setdefault(j, set()).add(i)
        if i == j:
            r = block.row_range
            diagonals[i] = np.arange(r.offset, r.offset + r.extent)

    components = []
    standard_deviation = np.empty(size)
    remaining = set(diagonals)
    smallest, largest = np.inf, 0.0
    dense_elements = 0
    while remaining:
        pending = [min(remaining)]
        group = set()
        while pending:
            i = pending.pop()
            if i not in group:
                group.add(i)
                pending.extend(neighbors[i] - group)
        remaining -= group
        group = sorted(group, key=lambda i: diagonals[i][0])
        indices = np.concatenate([diagonals[i] for i in group])
        if len(group) == 1:
            block = data[group[0], group[0]]
            diagonal = block.diagonal()
            nonzero = (
                block.count_nonzero()
                if sparse.issparse(block)
                else np.count_nonzero(block)
            )
            if nonzero == len(diagonal):
                sigma = np.sqrt(diagonal)
                standard_deviation[indices] = sigma
                components.append((indices, sigma))
                smallest, largest = min(smallest, 1.0), max(largest, 1.0)
                continue
        count = len(indices)
        dense_elements += count * count
        if dense_elements > max_dense_elements:
            raise ValueError(
                f"{name}: connected covariance factors exceed max_dense_elements"
            )
        matrix = np.zeros((count, count))
        offsets = {}
        offset = 0
        for i in group:
            offsets[i] = slice(offset, offset + len(diagonals[i]))
            offset += len(diagonals[i])
        for (i, j), block in data.items():
            if i not in offsets or j not in offsets:
                continue
            part = block.toarray() if sparse.issparse(block) else block
            matrix[offsets[i], offsets[j]] = part
            if i != j:
                matrix[offsets[j], offsets[i]] = part.T
        sigma = np.sqrt(np.diag(matrix))
        correlation = matrix / sigma[:, None] / sigma[None, :]
        # Validation has already rejected material asymmetry. Remove only the
        # accepted roundoff discrepancy before factorization/eigenanalysis.
        correlation = 0.5 * (correlation + correlation.T)
        factor = sigma[:, None] * linalg.cholesky(
            correlation, lower=True, check_finite=False
        )
        eigenvalues = linalg.eigvalsh(correlation, check_finite=False)
        smallest, largest = min(smallest, eigenvalues[0]), max(largest, eigenvalues[-1])
        standard_deviation[indices] = sigma
        components.append((indices, factor))
    condition = largest / smallest if smallest > 0 else np.inf
    return _CovarianceFactor(components, size, standard_deviation, float(condition))


def information(
    measurement_jac,
    model_state_covmat,
    measurement_vec_error_covmat,
    *,
    state_labels=None,
    measurement=None,
    prior_prediction=None,
    relative_tolerance=1e-10,
    max_dense_elements=10_000_000,
):
    r"""Validate covariance inputs and calculate local retrieval information.

    Parameters
    ----------
    measurement_jac : ``array_like``
        Jacobian with shape ``(m, n)`` in the actual retrieved coordinates
        and measurement units.
        It must describe the state at which this local analysis is intended.
    model_state_covmat, measurement_vec_error_covmat
        Prior and observation-error covariances, each a
        :class:`~pyarts3.arts.CovarianceMatrix` or ``array_like``.
        Square arrays are full
        covariances; one-dimensional arrays contain independent *variances*,
        not standard deviations. Native blocks and stored inverses are checked.
    state_labels : ``sequence of str``
        Optional labels, one per state coordinate in Jacobian-column order.
    measurement, prior_prediction : ``array_like``
        Optional arrays with shape ``(m,)``. Supply both to compare
        measurements with an explicit :math:`F(\vec{x}_a)`. The
        Jacobian must then be appropriate around :math:`\vec{x}_a`. A spectrum fitted at
        another state is not a prior prediction.
    relative_tolerance : float
        Relative, scaled tolerance for covariance symmetry/inverse checks.
    max_dense_elements : int
        Limit on estimated dense analysis storage and on dense covariance
        factors. The initial implementation uses an exact SVD and returns all
        n state modes. LAPACK can require additional temporary workspace.
        Independent measurement variances do not allocate an m-by-m matrix.

    Returns
    -------
    InformationReport
        Dimensionless information spectrum, physical-coordinate uncertainty
        estimates, and optional prior innovation statistics. This function
        does not evaluate the forward model or modify any of its inputs.

    Notes
    -----
    With covariance factors
    :math:`\mathbf{S}_a=\mathbf{L}_a\mathbf{L}_a^{\top}` and
    :math:`\mathbf{S}_\epsilon=\mathbf{L}_\epsilon\mathbf{L}_\epsilon^{\top}`,
    the singular values of the whitened Jacobian

    .. math::

        \widetilde{\mathbf{J}} = \mathbf{L}_\epsilon^{-1}\mathbf{J}\mathbf{L}_a

    are the mode signal-to-noise ratios. The posterior mode variance fractions
    are :math:`1/(1+s_i^2)`. Values apply to the assumed linear Gaussian model;
    high information does not establish correct uncertainties.
    No explicit covariance or posterior-precision inverse is formed.
    """
    max_dense_elements = operator.index(max_dense_elements)
    if max_dense_elements <= 0:
        raise ValueError("max_dense_elements must be positive")
    if not np.isfinite(relative_tolerance) or not 0 < relative_tolerance < 1:
        raise ValueError("relative_tolerance must be finite and between zero and one")
    jacobian = _real_array(measurement_jac, "measurement_jac")
    if jacobian.ndim != 2 or min(jacobian.shape) == 0:
        raise ValueError("measurement_jac must be a nonempty two-dimensional matrix")
    if not np.all(np.isfinite(jacobian)):
        raise ValueError("measurement_jac must contain only finite values")
    m, n = jacobian.shape
    if 3 * m * n + 4 * n * n > max_dense_elements:
        raise ValueError(
            "The exact information analysis exceeds max_dense_elements; "
            "reduce the state/measurement selection or explicitly raise the limit"
        )
    labels = (
        tuple(f"x[{i}]" for i in range(n))
        if state_labels is None
        else tuple(state_labels)
    )
    if len(labels) != n or not all(isinstance(label, str) for label in labels):
        raise ValueError("state_labels must contain one string per Jacobian column")
    if (measurement is None) != (prior_prediction is None):
        raise ValueError("Supply both measurement and prior_prediction, or neither")
    residual = None
    if measurement is not None:
        measured = _real_array(measurement, "measurement")
        predicted = _real_array(prior_prediction, "prior_prediction")
        if measured.shape != (m,) or predicted.shape != (m,):
            raise ValueError(
                "measurement and prior_prediction must match the Jacobian rows"
            )
        if not np.all(np.isfinite(measured)) or not np.all(np.isfinite(predicted)):
            raise ValueError("measurement and prior_prediction must be finite")
        with np.errstate(over="ignore", invalid="ignore"):
            residual = measured - predicted

    prior = _factor_covariance(
        model_state_covmat,
        n,
        "model_state_covmat",
        relative_tolerance,
        max_dense_elements,
    )
    noise = _factor_covariance(
        measurement_vec_error_covmat,
        m,
        "measurement_vec_error_covmat",
        relative_tolerance,
        max_dense_elements,
    )
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        whitened = noise.solve_left(prior.multiply_right(jacobian))
        # Left and right operations commute. A different order can avoid an
        # overflowing intermediate when their physical scales cancel.
        if not np.all(np.isfinite(whitened)):
            whitened = prior.multiply_right(noise.solve_left(jacobian))
    if not np.all(np.isfinite(whitened)):
        raise ValueError(
            "Whitened Jacobian is not finite; "
            "check covariance scales and coordinate units"
        )
    u, observed_s, vt = linalg.svd(whitened, full_matrices=m < n, check_finite=False)
    if not np.all(np.isfinite(observed_s)):
        raise ValueError("Information spectrum exceeds numerical range")
    singular_values = np.zeros(n)
    singular_values[: len(observed_s)] = observed_s
    factors = np.hypot(1.0, singular_values)
    mode_reduction = (singular_values / factors) ** 2
    mode_bits = _mode_information_bits(singular_values)
    state_modes = prior.multiply_left(vt.T)
    posterior_sd = np.hypot.reduce(state_modes / factors, axis=1)
    marginal_reduction = np.clip(
        1.0 - (posterior_sd / prior.standard_deviation) ** 2, 0, 1
    )

    innovation = None
    if residual is not None:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            whitened_residual = noise.solve_left(residual)
        if not np.all(np.isfinite(whitened_residual)):
            raise ValueError(
                "Whitened prior innovation is not finite; check values and scales"
            )
        projected = u.T @ whitened_residual
        complement = whitened_residual - u @ projected if m > n else np.zeros(m)
        innovation = float(
            np.dot(complement, complement)
            + np.sum((projected / factors[: len(observed_s)]) ** 2)
        )
        if not np.isfinite(innovation):
            raise ValueError(
                "Prior innovation statistic overflowed; check values and scales"
            )

    return InformationReport(
        singular_values=_readonly(singular_values),
        mode_variance_reduction=_readonly(mode_reduction),
        degrees_of_freedom=float(np.sum(mode_reduction)),
        information_bits=float(np.sum(mode_bits)),
        prior_standard_deviation=_readonly(prior.standard_deviation),
        posterior_standard_deviation=_readonly(posterior_sd),
        variance_reduction=_readonly(marginal_reduction),
        state_modes=_readonly(state_modes),
        measurement_modes=_readonly(u),
        state_labels=labels,
        prior_correlation_condition=prior.condition,
        measurement_correlation_condition=noise.condition,
        innovation_chi_square=innovation,
        innovation_expected_mean=m if innovation is not None else None,
        _noise=noise,
    )


def _workspace_state_metadata(ws, size):
    """Use target offsets, not category order, to identify Jacobian columns."""
    labels = [[] for _ in range(size)]
    blocks = []
    if ws.has("jac_targets"):
        for category in ("atm", "surf", "subsurf", "line", "sensor", "error"):
            for target in getattr(ws.jac_targets, category):
                start, count = int(target.x_start), int(target.x_size)
                if start < 0 or count <= 0 or start + count > size:
                    raise ValueError(
                        f"Jacobian target {category}.{target.type} has invalid state "
                        f"range [{start}:{start + count}] for {size} columns; "
                        "finalize targets and recompute the Jacobian"
                    )
                name = f"{category}.{target.type}"
                blocks.append((name, start, count))
                for offset in range(count):
                    labels[start + offset].append(f"{name}[{offset}]")
    # Overlapping targets can intentionally share a state coordinate.
    return (tuple(" / ".join(names) if names else f"x[{i}]"
                  for i, names in enumerate(labels)),
            tuple(sorted(blocks, key=lambda block: (block[1], block[0]))))


def information_from_workspace(ws, *, prior_prediction=None, **options):
    """Analyze the Jacobian and covariances already present in a workspace.

    No agenda is executed. Supply ``prior_prediction`` explicitly to enable
    innovation checking using ``ws.measurement_vec``. The caller must ensure
    that the Jacobian describes the intended state and coordinates.
    Field labels and ``state_blocks`` (name, start, size) are inferred from
    ``jac_targets``. Offsets refer to the supplied Jacobian's columns; field
    indices are flattened target indices, not altitude or physical units.
    Keys alone do not identify logarithmic or other coordinate transforms.
    Missing target metadata leaves generic ``x[i]`` labels. Explicit
    ``state_labels`` override the automatic per-coordinate labels.
    Remaining options are forwarded to :func:`information`.
    """
    if "measurement" in options:
        raise TypeError(
            "information_from_workspace takes measurements from ws.measurement_vec"
        )

    assert ws.has("measurement_jac"), "Jacobian not present in workspace"
    assert ws.has("model_state_covmat"), "Prior covariance not present in workspace"
    assert ws.has("measurement_vec_error_covmat"), "Measurement error covariance not present in workspace"
    assert ws.has("measurement_vec") or prior_prediction is not None, "Measurement vector not present in workspace"

    labels, blocks = _workspace_state_metadata(ws, np.asarray(ws.measurement_jac).shape[1])
    if options.get("state_labels") is None:
        options["state_labels"] = labels
    report = information(
        ws.measurement_jac,
        ws.model_state_covmat,
        ws.measurement_vec_error_covmat,
        measurement=ws.measurement_vec if prior_prediction is not None else None,
        prior_prediction=prior_prediction,
        **options,
    )
    return replace(report, state_blocks=blocks)
