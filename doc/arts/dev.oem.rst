.. _sec-development-oem:

Developing the OEM interface
============================

The user guide is in :ref:`sec-user-oem` and the mathematical formulation is
in :ref:`Sec OEM`.  This page separates the initial
interface cleanup from further work on numerical algorithms and public APIs.
The statistical objective and existing method names should remain stable
unless a change is explicitly documented and tested.

Responsibilities and regression baseline
----------------------------------------

``src/m_oem.cc`` validates inputs and connects workspace inputs and outputs
to the inverse problem.  ``src/oem.h`` contains the agenda adapter, solver
wrappers, and logging.  ``3rdparty/invlib/src/invlib/map.cpp`` drives
the optimization and computes costs and convergence measures.  The LM
step and linear solvers live below that layer.  Covariance storage and
inverse operations are implemented in
``src/core/jacobian/covariance_matrix.cc``; retrieval covariance setup
also involves ``src/m_covmat.cc``.

The first cleanup centralizes method selection and LM configuration while
retaining the public defaults.  Direct and CG LM methods now use all six
settings, the supplied ``stop_dx``, and the same diagonal of prior precision
for damping.  Unsupported ``li_m`` and ``gn_m`` are rejected before the
forward model runs.  State normalization is rejected for measurement-space
methods instead of applying state scales to a measurement-space system.

Measurement-space convergence uses the state-space Hessian metric
(Rodgers 5.30), so it remains well-defined when the numbers of measurements
and states differ.  Diagnostics come from the formulation that actually
ran.  LM history is collected independently of terminal output, and aliases
have the same initial-cost behavior.  Final nonlinear Jacobians are refreshed
at the returned state before calculating the gain.  Outputs that would
otherwise retain an earlier run's gain or errors are cleared, and the
``clear_matrices`` policy also applies to a skipped retrieval.

``src/tests/test_oem_methods.cc`` exercises the workspace interface using
small deterministic forward models.  CTest registers the individual cases
as ``cpp.fast.oem.<method>``, plus settings and validation cases.  Each
supported spelling, including ``ml`` and ``ml_cg``, is part of the contract.
These tests run without the optional MATLAB integration required by the
legacy invlib-level tests in ``src/tests/test_oem.cc``.  They complement
the atmospheric retrieval examples under ``tests/core/vmr``.  Exact initial
solutions and underdetermined problems are covered; a zero right-hand side
now returns a zero CG solution directly.

The numerical reference for an affine model should come from an independently
formed linear posterior solution.  Check the state, fitted measurement,
cost decomposition, gain, and output dimensions; agreement between two
methods alone can preserve a shared bug.  Include correlated covariances,
unequal state and measurement dimensions, a non-prior starting state, and
equivalent scaled and unscaled systems.  Nonlinear cases need to test
accepted-state Jacobians and meaningful damping changes, rather than only
checking that a residual decreased.

Bounded inner iterations
------------------------

An outer ``max_iter`` does not bound the work of an inner linear solve or
LM trial search.  Both now enforce independent limits and raise an error
when they cannot complete within those limits.  OEM maps errors caught
during inversion to status 9 and records the explanation in ``errors``.
Exhausting the existing LM damping range retains its status 2 behavior.
Gauss--Newton must propagate linear-solver exceptions with their nested
cause.  Returning an empty step after catching an error hides the failure
and can cause invalid vector operations in the measurement-space formulation.

The native ``ConjugateGradient`` constructor and both preconditioned
variants accept ``max_iterations`` after the verbosity argument, with a
default of 1000.  The solve loop counts iterations locally, so a custom
convergence predicate cannot disable the bound.  Checks for non-finite
arithmetic and non-positive curvature reject invalid or broken solves
instead of treating them as converged.  These checks do not replace
covariance validation.  OEM retains its relative residual tolerance of
``1e-10`` and uses the native iteration limit.

The native LM optimizer provides ``get_maximum_trials()`` and
``set_maximum_trials()``; the positive trial limit defaults to 100 per
outer iteration.  A local counter bounds rejected-step retries, and a
damping update must make progress before another trial is attempted.
This catches multiplication that rounds back to the current damping.
The trial budget counts linear solves, including any additional undamped
solve used to check stationarity.
Neither limit is currently an OEM workspace argument or part of the six
named LM damping settings.  Preserve the termination regressions when
changing convergence predicates, trial acceptance, or damping updates.

LM acceptance and stop outcomes
-------------------------------

``LMStopReason`` records termination independently of the damping value.
``None`` means the optimizer can continue; ``Stationary`` is successful
termination.  ``DampingLimit`` returns a zero step and maps to workspace
status 2.  ``TrialLimit``, ``DampingStalled``, ``LinearSolverFailure``, and
``NumericalFailure`` accompany exceptions, which OEM maps to status 9
when caught during inversion.  Ordinary convergence after an accepted
step still uses the configured convergence criterion and damping gate.
Both that criterion and ``Stationary`` map to status 0; outer iteration
exhaustion retains status 1.

All MAP formulations consult the optimizer's explicit outcome before
testing the returned step for convergence.  A zero step returned after
damping exhaustion must not count as convergence.  Rejected trials must
never update the accepted state, and damping history must contain physical
damping values, not a numerical sentinel such as ``maximum + 1``.  At a
large finite maximum, that addition can round back to the maximum itself.
The ordinary convergence gate returns a zero tolerance while damping
exceeds its limit, so the strict comparison also rejects a rounded-zero
state change.

The generic helper in ``optimization/minimize.h`` also honors optional
``stop_iteration()`` and ``converged()`` hooks.  A reported optimizer
failure returns before testing the criterion or updating the state.  A
stationary optimizer's verified step is applied, but success still requires
``J.criterion`` to meet the helper's own tolerance.  Iteration exhaustion
returns one.  Retain native regressions for these outcomes and for custom
minimizers without either hook.

``MAPBase::model_cost_scale()`` returns two: MAP reports the full quadratic
cost, while its normal equations use the half-gradient and half-Hessian.
The LM prediction includes this scale; an exact affine model therefore
has an actual-to-predicted reduction ratio of one.  Generic objectives
without this optional method retain scale one and must provide mutually
consistent cost, gradient, and Hessian definitions.  Keep regressions for
both conventions when changing the reduction formula or acceptance
thresholds.

LM compares reductions against ``32 * epsilon * max(abs(old_cost),
abs(new_cost))`` before dividing them.  If actual or predicted reduction
is unresolved at that scale, it checks the undamped normal equations once
within the current step.  Numerical stationarity requires a finite,
nonnegative undamped decrement below ``n * stop_dx``, a finite,
nonnegative predicted reduction within the current cost's roundoff
scale, and a candidate cost indistinguishable from the current cost at
that scale.  This step can converge independently of the damping gate:
its size is checked without damping.  An exactly zero gradient is also
stationary.  These checks must not substitute the damped step for the
undamped decrement; very strong damping can hide a large remaining error.
Resolved acceptance requires finite, positive actual and predicted
reductions and a ratio of at least 0.5.

Keep workspace regressions for all four LM spellings at damping
``1e20``, tight affine tolerances, and an exactly stationary state with
nonzero cost.  Native tests also cover explicit stop reasons, objective
scaling, rejected steps, and bounded or stalled damping updates.  Check
the state, fitted measurement, costs, gain, and damping history as well
as the final status.

Remaining numerical work
------------------------

The following items require separate implementation and regression work.

1. **Expose inner-solver controls and diagnostics.**  Provide public
   settings for relative CG tolerance, maximum linear iterations, and
   maximum LM trials, and report linear iterations and residuals.  Keep
   these separate from outer ``max_iter`` and the six damping controls.
   Extend numerical coverage to nearly zero right-hand sides and
   ill-conditioned positive-definite systems; bounded termination alone
   does not establish the accuracy of a difficult solve.

2. **Extend public stop diagnostics.**  The workspace retains its numeric
   status mapping; native LM stop reasons are not a separate workspace
   output.  A future result could expose the precise stop reason and
   distinguish successful completion of a one-step linear solve from
   exhaustion of an iterative convergence budget.

3. **Replace dense normalization and explicit gain inversion.**  Store
   state scales as a vector and apply diagonal products without allocating
   an :math:`n\times n` dense scaling matrix.  Compute the gain by solving
   the posterior precision system for its right-hand sides.  Measure peak
   memory as well as runtime, and retain the option to skip gain production.
   Consider factorization reuse and a preconditioner appropriate to each
   formulation only after correctness is established.

4. **Make agenda state and failure behavior explicit.**  Validation,
   covariance inversion, and the initial agenda evaluation can still throw
   directly; status 9 covers errors caught during inversion.  Trial forward
   model evaluations mutate atmosphere, sensor, and surface workspace data.
   Define which state and diagnostics remain valid after an exception or a
   rejected trial.  A future result should distinguish the last accepted
   state from a failed trial without treating NaNs as the only status record.

Covariance validation and construction
--------------------------------------

``CovarianceMatrix.validate(expected_size=-1, relative_tolerance=1e-10,
max_dense_elements=10_000_000)`` checks an existing physical covariance
without modifying it.  Validation covers block dimensions and ranges,
coverage, unique upper-triangular block identifiers, non-null matrices,
finite values, positive variances, symmetry, and positive definiteness.
Symmetry and supplied-inverse comparisons use variance-scaled coordinates
so mixed units do not let a large-variance quantity hide an invalid
small-variance block.  Positive definiteness is checked over connected
components of the complete covariance, including cross-block correlations.

Each represented inverse component must be complete and consistent with
its connected covariance component.  Whole independent components may
remain uncached during validation.  In normalized coordinates, the matrix-product
residual is compared with ``relative_tolerance`` times the component size;
the diagonal fast path uses ``relative_tolerance`` directly.
``compute_inverse`` uses the same checks and inverts uncached independent
components.  This catches invalid
physical covariances on the existing OEM inversion path without imposing
the standalone analysis's memory limit on existing retrievals.  It preserves
the inverse-only representation used internally for LM damping; that
representation is not accepted as a physical covariance by ``validate``
or the standalone information report.

Replacing blocks or accessing mutable blocks invalidates cached inverses.
Adding a covariance block invalidates the affected connected components
while retaining independent cached components.  Shared matrix aliases
can still outlive those access points; validation before inversion must
continue to detect a stale supplied inverse.  Do not silently overwrite
asymmetric covariance entries during inversion.  The covariance-addition
helper in ``src/m_covmat.cc`` now inserts the optional supplied inverse
as an inverse and validates its inputs before mutation.

Sparse diagonal validation and inversion retain diagonal storage.
Explicit validation of non-diagonal connected components requires dense
work; check
``max_dense_elements`` before allocating each component matrix.  This
guard bounds individual dense arrays, not total peak memory.  Keep
regressions for diagonal and correlated blocks, mixed scales, incomplete
coverage, inverse-only inputs, stale caches, and supplied inverse
consistency.  Validation must not silently repair a statistical model.

Useful construction helpers would accept standard deviations and an optional
correlation matrix or named correlation kernel, then construct covariance
blocks.  Parameter names should state whether they accept a standard
deviation, variance, covariance, or precision.  A report could identify each
block by retrieval target and include its units, variance range, and
factorization or conditioning diagnostics.  Covariance repair, such as
adding a diagonal term or changing correlations, must be an explicit user
choice because it changes the inference.

Standalone information analysis
--------------------------------

``python/src/pyarts3/retrieval.py`` provides ``information`` and
``information_from_workspace`` independently of the OEM optimizer.
The workspace adapter reads the current Jacobian and the two physical
covariances.  Neither entry point executes a forward agenda, runs OEM,
updates a workspace, or selects new covariance values.  This separation
lets a user compare proposed measurement and uncertainty models with a
fixed Jacobian before a retrieval, and examine local sensitivities after
one.  The mathematical definitions belong in :ref:`sec-oem-information`;
usage and interpretation belong in :ref:`sec-user-oem-information`.

The report uses covariance factors to obtain a dimensionless Jacobian,
then an exact singular value decomposition.  ``singular_values`` and
``mode_variance_reduction`` include all state directions, padding the
unobserved modes with zero singular values when measurements are fewer
than states.  ``state_modes`` contains physical state directions scaled
to unit prior uncertainty; ``measurement_modes`` refers to whitened
measurement coordinates.  Retain the distinction between mode and
marginal variance reduction.  Mode signs and bases inside a degenerate
subspace are not stable identifiers for comparisons between runs.

Use factor solves rather than forming covariance inverses for this
analysis.  Diagonal measurement covariances must not allocate a full
measurement-square matrix.  Guard dense work with ``max_dense_elements``;
even a diagonal prior still needs a state-square basis when all modes
are returned.  The analysis guards an estimate of ``3*m*n + 4*n*n``
elements as well as the cumulative dense factor sizes; LAPACK can need
additional work arrays.  Future truncated or operator-based analyses must state
which parts of the spectrum and uncertainty they approximate, and must
not label retained-mode information as the complete information content.

``information`` requires both an explicit measurement and
``prior_prediction`` for ``innovation_chi_square``.
``information_from_workspace`` reads the workspace measurement only when
given ``prior_prediction``; it does not accept a measurement override.
Do not infer the prior prediction from a workspace's fitted measurement:
its generating state is not known to this API.  The
chi-squared reference distribution assumes a linear model, a prediction
at the prior mean, and the stated independent Gaussian error sources.
The spectrum itself is independent of the measured residual.

Regression oracles should include diagonal analytic systems, correlated
covariances, unequal dimensions, exact null directions, consistent
coordinate changes, and posterior covariance reconstructed from the
returned modes.  Cross-check the innovation statistic against an
independently formed innovation covariance.  Keep workspace immutability
and dense-memory guards covered.  Plot uncertainty ratios when combining
state elements with different units; retain absolute standard deviations
and labels in the numerical report.

Helping users choose settings
-----------------------------

Named LM settings are implemented by
:class:`~pyarts3.arts.OEMLMSettings`.  Its keyword-only constructor and
mutable fields expose ``initial_damping``, ``decrease_factor``,
``increase_factor``, ``maximum_damping``, ``damping_threshold``, and
``convergence_damping_limit`` in the existing vector order.  The native
settings representation provides shared validation for the Python object
and the legacy vector used by ``OEM``.  Both direct and CG optimizers must
continue to use that same mapping.
The type lives in ``src/core/jacobian/oem_settings.h`` and
``src/core/jacobian/oem_settings.cc`` without an invlib dependency;
``src/python_interface/py_retrieval.cpp`` binds the Python interface.

The Python constructor validates all six settings, and each field setter
validates a temporary copy before replacing the stored object.  An invalid
edit must preserve the previous values and report the affected setting.
This also avoids losing useful validation errors through nanobind's
implicit-conversion error handling.  Related field changes must either
keep each intermediate configuration valid or use a replacement object
constructed with the desired keyword arguments together.

``validate()`` checks a configured object, ``as_vector()`` validates and
exports its six values, and ``from_vector()`` validates and imports an
existing vector.  Conversion and consumption at the workspace call boundary
retain validation as well.  ``describe()`` explains the active values;
the representation displays every field.  Keep these descriptions, field
validation, and the positional mapping consistent when extending the API.
The Python object is not a workspace group; workspace and XML storage
continue to use the vector representation.

The object's defaults are an explicit starting configuration.  The empty
default of the ``OEM`` argument is unchanged, so existing calls do not
silently select new damping values.  Do not describe the threshold as a
hard minimum: it controls both restart after rejection and switching a
proposed decrease to zero.  The convergence damping limit gates the
existing state-step criterion using the updated damping.  Preserve the
accepted/rejected-step regressions when changing either behavior.
A preset claiming particular convergence or performance properties needs
documented assumptions and representative retrieval benchmarks.

``src/tests/test_oem_methods.cc`` tests the native settings and solver
behavior.  ``tests/core/jac/oem_lm_settings.py`` covers the Python interface
and its conversion through the actual workspace call.  Keep both named
and legacy paths covered for all four LM spellings, including correlated
priors and nonlinear rejected trials; comparing complete results and
damping histories catches changes that an endpoint-only test would miss.

Further API work should keep existing scripts usable:

* Represent iteration family, linear solver, and formulation separately
  internally.  A single method descriptor table can drive parsing,
  supported-name documentation, and test enumeration.  Existing string
  spellings and aliases can remain as the public compatibility layer.
* Provide a setup report with state and measurement sizes, estimated
  matrix storage, chosen coordinates and scales, covariance diagnostics,
  and solver settings.  Recommend a method with a reason while leaving
  the explicit choice available.  Dimensions can suggest a formulation;
  they cannot establish linearity or predict nonlinear convergence.
* Add an optional local Jacobian check using representative perturbations
  in prior-scaled coordinates, and report discrepancies by measurement
  and retrieval target.  Start with a deterministic small-model example
  before applying it to expensive radiative-transfer agendas.
* Return structured diagnostics with accepted iteration costs, state-step
  measures, damping, trial rejections, linear residuals, and forward-model
  call counts.  Preserve the five-element diagnostic vector as an adapter.
  A user can then distinguish poor model linearity, a difficult linear
  solve, an exhausted iteration budget, and a mismatch in assumed errors.

Method recommendations should be assessed on linear, weakly nonlinear,
strongly nonlinear, correlated, poorly scaled, and large sparse cases.
Compare final objective and state as well as runtime, peak memory, and
forward-model calls.  Automatic changes to prior or measurement covariance
require scientific assumptions that cannot be inferred from solver progress
alone.  Keep those assumptions visible in both the setup and uncertainty
reports.
