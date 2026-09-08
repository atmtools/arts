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

Remaining numerical work
------------------------

The following items require separate implementation and regression work.

1. **Make CG termination explicit.**  Expose relative tolerance and maximum
   linear iterations, and report linear iterations and residuals.  Handle
   non-finite arithmetic, zero curvature, and loss of positive definiteness
   explicitly.  Extend the exact-initial-state regression to nearly zero
   right-hand sides and ill-conditioned positive-definite systems.
   An outer ``max_iter`` does not bound the work of an inner linear solve.

2. **Audit LM's predicted cost reduction.**  ``MAP::cost_function`` returns
   the full quadratic objective, while the gradient and Hessian assembled
   for a step correspond to half that objective.  The reduction ratio in
   ``levenberg_marquardt.cpp`` uses these step quantities without matching
   the factor of two.  An exact linear model therefore gives a ratio of two
   instead of one.  Characterize accepted and rejected trials before
   correcting this: changing the ratio affects damping trajectories and
   the established acceptance thresholds, even though the minimizer's
   objective is unchanged.

   Also handle near-zero predicted reduction and roundoff stagnation.  In
   an affine regression with ``stop_dx=1e-12``, LM can reach the reference
   solution while the preceding step's convergence measure is still above
   tolerance.  The next nearly zero step then gives an unreliable reduction
   ratio and can exhaust the damping limit despite an optimal state.  The
   baseline uses ``stop_dx=1e-9`` with independent tight state and cost
   checks.  A stricter-tolerance regression needs a defined numerical
   termination policy before changing the acceptance logic.

3. **Use explicit stop reasons.**  LM currently signals its damping limit
   by setting gamma to ``gamma_max + 1``.  At sufficiently large finite
   values, this rounds back to ``gamma_max``.  Replace numerical sentinels
   with a stop reason carried by the optimizer.  Keep the legacy numeric
   status as a compatibility mapping.  Likewise, distinguish successful
   completion of a one-step linear solve from exhaustion of an iterative
   convergence budget.

4. **Replace dense normalization and explicit gain inversion.**  Store
   state scales as a vector and apply diagonal products without allocating
   an :math:`n\times n` dense scaling matrix.  Compute the gain by solving
   the posterior precision system for its right-hand sides.  Measure peak
   memory as well as runtime, and retain the option to skip gain production.
   Consider factorization reuse and a preconditioner appropriate to each
   formulation only after correctness is established.

5. **Make agenda state and failure behavior explicit.**  Validation,
   covariance inversion, and the initial agenda evaluation can still throw
   directly; status 9 covers errors caught during inversion.  Trial forward
   model evaluations mutate atmosphere, sensor, and surface workspace data.
   Define which state and diagnostics remain valid after an exception or a
   rejected trial.  A future result should distinguish the last accepted
   state from a failed trial without treating NaNs as the only status record.

Covariance validation and construction
--------------------------------------

The current covariance inverse implementation uses generic matrix inversion;
it does not establish that a covariance is a valid statistical model.
A validation API should report dimensions, block coverage, finite entries,
positive variances, symmetry, and positive definiteness.  Check connected
blocks and use matrix structure so that validation does not require
unconditionally materializing a full dense covariance.

Two storage issues need dedicated tests.  In ``src/m_covmat.cc``, the optional
inverse branch of ``add_diagonal_covmat`` inserts the covariance again
instead of adding the supplied inverse.  Separately, an inverse-only
``CovarianceMatrix`` is not interchangeable with one containing covariance
blocks: measurement-space methods also multiply by the covariance, and
missing forward blocks can silently yield zero contributions.  Specify
supported storage forms, reject incomplete representations, and make cache
invalidation reliable when covariance blocks change.  Test supplied inverse
consistency against complete correlated covariances, not independent inverses
of their component blocks.

Useful construction helpers would accept standard deviations and an optional
correlation matrix or named correlation kernel, then construct covariance
blocks.  Parameter names should state whether they accept a standard
deviation, variance, covariance, or precision.  A report could identify each
block by retrieval target and include its units, variance range, and
factorization or conditioning diagnostics.  Covariance repair, such as
adding a diagonal term or changing correlations, must be an explicit user
choice because it changes the inference.

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
