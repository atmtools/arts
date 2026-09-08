.. _sec-user-oem:

Configuring an optimal-estimation retrieval
###########################################

This guide covers method selection, settings, covariance setup, and retrieval
diagnostics.  The mathematical formulation is in :ref:`Sec OEM`.

Choosing a method
=================

:meth:`~pyarts3.workspace.Workspace.OEM` always minimizes the same objective,
including the a priori term.  Choosing a method changes how this objective
is minimized.  It does not change the statistical meaning of the two
covariance matrices.

There are two decisions: whether the forward model needs nonlinear iteration,
and how to solve the linear system inside each iteration.  In the table,
``n`` is the number of retrieved state elements and ``m`` is the
number of measurements.

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - ``method``
     - Iteration and linear solver
     - When to consider it
   * - ``li``
     - One Gauss--Newton step; direct solve in state space
     - The forward model is linear over the relevant state range.
   * - ``li_cg``
     - One step; conjugate gradient (CG) in state space
     - A large linear problem where a direct solve is expensive.
   * - ``li_cg_m``
     - One step; CG in measurement space
     - A linear problem with substantially fewer measurements than states.
   * - ``gn``
     - Gauss--Newton; direct solve in state space
     - A mildly nonlinear problem with a plausible starting state.
   * - ``gn_cg``
     - Gauss--Newton; CG in state space
     - The same nonlinear problem when the linear solve is the bottleneck.
   * - ``gn_cg_m``
     - Gauss--Newton; CG in measurement space
     - Fewer measurements than states, after checking against ``gn`` on a
       smaller representative problem.
   * - ``lm``
     - Levenberg--Marquardt (LM); direct solve in state space
     - Gauss--Newton steps are too large or the starting state is uncertain.
   * - ``lm_cg``
     - LM; CG in state space
     - LM is needed and the direct linear solve is expensive.

``ml`` is an alias for ``lm``, and ``ml_cg`` is an alias for ``lm_cg``.
In particular, ``ml`` does not select maximum-likelihood estimation or
remove the prior.  ``li_m`` and ``gn_m`` are unsupported.

Start development of a retrieval with a direct method on a small, representative
case.  Use ``li`` only after checking linearity; a single step applied to a
nonlinear model is only a local approximation.  Use ``gn`` when this
linearization is reliable and try ``lm`` when its unconstrained steps fail.
Check the Jacobian and units before tuning damping to compensate for poor
convergence.

The state-space system has size ``n`` by ``n``; the measurement-space
system has size ``m`` by ``m``.  Size is only a first guide to performance:
covariance structure, conditioning, and forward-model cost also matter.
CG has a fixed internal relative residual tolerance of ``1e-10``;
``stop_dx`` controls the outer retrieval iteration, not this linear solve.
OEM still stores the measurement Jacobian, and computing the gain matrix
requires additional dense matrices.  ``clear_matrices=1`` skips the gain
calculation and returns empty Jacobian and gain matrices when those outputs
are not needed.

Choosing covariance matrices
============================

Set the state coordinates, measurement units, and ordering first.  Then
construct :attr:`~pyarts3.workspace.Workspace.model_state_covmat`
(shape ``n`` by ``n``) and
:attr:`~pyarts3.workspace.Workspace.measurement_vec_error_covmat`
(shape ``m`` by ``m``) in those coordinates.
Both must be finite, symmetric, and positive definite for these retrievals.
Positive diagonal entries alone do not guarantee positive definiteness.

The diagonal of a covariance matrix contains **variances**, not standard
deviations or inverse variances.  An uncertainty of 5 K corresponds to
25 K\ :sup:`2`.  An absolute VMR uncertainty of ``1e-6`` corresponds
to a variance of ``1e-12``.  An off-diagonal entry has units equal to
the product of the units of its two coordinates.

For independent measurement errors with standard deviation ``sigma`` in the
units of ``measurement_vec``, the constant-noise helper takes the variance:

.. code-block:: python

   ws.measurement_vec_error_covmatConstant(value=sigma**2)

This helper assigns the same variance to every measurement and zero
cross-covariance.  Use a full or block covariance when channels have different
uncertainties or correlated errors.  Include uncertain forward-model
contributions that are represented as observation error; add covariance
contributions only when their errors are independent.  If a nuisance parameter
is retrieved explicitly, account for it in the state and prior without also
counting the same uncertainty independently in the measurement covariance.

The prior covariance describes uncertainty about the state before using this
measurement.  Increasing it weakens the pull towards the prior; increasing
the measurement covariance weakens the pull towards the measurement.
Choose these uncertainties from the measurement and prior error models.
Changing them to obtain a preferred fit changes the inference and its reported
uncertainty.  A parameter that must be fixed should be excluded from the state;
zero variance would make the covariance singular.

Construct correlated covariances from standard deviations and a valid
correlation matrix, as described in :ref:`sec-oem-covariance`.  A correlation
length can express how nearby profile elements vary together.  Long correlations
constrain differences between those elements and affect the retrievable
resolution.  Perfect correlation creates a singular matrix.  Correlations
between retrieved quantities also belong in the prior covariance when
supported by the prior model; independently populated diagonal blocks assume
those quantities are uncorrelated.  A dense diagonal *block* can itself contain
correlated elements.

The inverse covariance is the **precision** matrix.  Do not pass a precision
matrix as a covariance.  For correlated variables, its diagonal generally
differs from the reciprocals of the covariance diagonal.  ARTS can store
covariance blocks and inverse blocks separately; a provided inverse must represent
the same complete covariance, including its correlations.  Supplying the
covariance and allowing ARTS to compute its inverse avoids having two
potentially inconsistent descriptions.

Coordinate changes and numerical scaling
----------------------------------------

The prior must use the coordinates actually retrieved.  Relative retrievals
use dimensionless relative variances, and logarithmic retrievals use variances
in log coordinates.  A 10 percent uncertainty in a relative coordinate has
variance ``0.1**2``; this is also a small-error approximation for a log-relative
coordinate.  For large uncertainties, define the intended prior distribution
in the transformed coordinates explicitly.  The covariance transformation is
given in :ref:`sec-oem-covariance`.  Changing measurement units also requires
transforming the measurement covariance and Jacobian consistently.

``model_state_covmat_normalization`` is a separate numerical scaling of
the linear system.  It does not change the retrieved coordinate system,
the supplied covariances, or the mathematical objective.  Leave it empty
to disable scaling, or supply ``n`` finite, strictly positive scales.
Prior standard deviations (the square roots of the prior variances)
are a useful initial choice for mixed units
such as temperature and absolute VMR.  This scales state *increments*,
so a zero prior mean is not a reason to use a zero scale.  Check that the
retrieved state agrees with the unscaled solution within numerical accuracy.
State normalization is unsupported for ``li_cg_m`` and ``gn_cg_m``;
these methods solve in measurement space.

Setting LM damping
==================

LM adds damping to limit the size of a Gauss--Newton step.  Larger gamma
penalizes larger steps, with the penalty scaled by the diagonal of the prior
precision matrix.  A zero value gives the Gauss--Newton step.  LM adjusts
damping by comparing actual and predicted cost changes, and may try several
forward-model evaluations within one outer iteration.  The damped system is
defined in :ref:`sec-oem-damping`.

``lm_ga_settings`` is required for all LM names; its default is an empty
vector, not a preset.  Direct and CG variants use the same six entries:

.. list-table::
   :header-rows: 1
   :widths: 8 24 68

   * - Index
     - Suggested local name
     - Meaning
   * - 0
     - ``gamma_start``
     - Initial damping, at least zero and no greater than the maximum.
   * - 1
     - ``decrease_factor``
     - Divisor when damping is reduced; must be greater than one.
   * - 2
     - ``increase_factor``
     - Multiplier when damping is increased; must be greater than one.
   * - 3
     - ``gamma_max``
     - Positive maximum damping.  Failure to find an acceptable step at
       this value stops the retrieval.
   * - 4
     - ``gamma_threshold``
     - Positive restart value when a zero-damping step fails.  When a
       proposed decrease would fall below this value, damping becomes zero.
       Must not exceed the maximum.
   * - 5
     - ``gamma_convergence``
     - Nonnegative upper damping limit for enabling the ordinary
       ``stop_dx`` criterion.  This refers to the current damping after its
       update, not the lowest value used in an earlier accepted iteration.

All entries must be finite.  Reduction is a **division** by entry 1:
setting it to 0.5 would increase damping and is invalid.  Not every accepted
step reduces damping; the adjustment also depends on the quality of the
local model and whether trial steps were rejected.

For an already configured retrieval, an explicit starting point is:

.. code-block:: python

   gamma_start = 10.0
   decrease_factor = 2.0
   increase_factor = 2.0
   gamma_max = 100.0
   gamma_threshold = 1.0
   gamma_convergence = 0.0
   ws.OEM(
       method="lm",
       max_iter=20,
       lm_ga_settings=[
           gamma_start, decrease_factor, increase_factor,
           gamma_max, gamma_threshold, gamma_convergence,
       ],
   )

These values are an example to assess on representative cases, not universal
defaults.  A zero convergence limit enables the ordinary stopping test once
damping reaches zero.  A larger limit permits termination while damping
still constrains the steps; a small damped step can then hide a remaining
distance to the minimum.  Raise the initial damping if early trial states
are consistently too aggressive.  If damping repeatedly reaches its maximum,
check the forward model, Jacobian, state parameterization, and covariance
scales before increasing that maximum.

Checking the result
===================

``oem_diagnostics`` has five entries with zero-based indices:

.. list-table::
   :header-rows: 1

   * - Index
     - Meaning
   * - 0
     - Status: 0 means the convergence criterion was met, 1 means the
       iteration budget was exhausted, 2 means the LM damping limit was
       reached, 9 means an error was caught during inversion, and 99 means
       the starting cost exceeded ``max_start_cost``.
   * - 1
     - Starting total cost, divided by the number of measurements.
   * - 2
     - Final total cost, divided by the number of measurements.
   * - 3
     - Final measurement contribution to the cost, divided by the number of measurements.
   * - 4
     - Number of outer iterations.

Unavailable values are NaN.  ``max_start_cost`` is a limit on the total
cost at the starting state; its default is infinity.  A value at or below
zero disables this limit and can skip computation of the starting cost
for non-LM methods when progress output is off.  It is not an upper bound
on the final residual.  The ``li`` variants perform exactly one step and
can return status 1 even for an exact linear solution, because they do not
take a second step to establish iterative convergence.

``stop_dx`` defaults to 0.01 and controls the weighted state-step convergence
measure described in :ref:`sec-oem-convergence`.  It does not set a relative
change in cost or an unweighted state difference.  Lower it for a stricter
convergence requirement and increase ``max_iter`` if useful progress continues
at the iteration limit.

An excessively small ``stop_dx`` can demand changes below the numerical
accuracy of the forward model or linear solve.  LM may then exhaust its
damping limit even when the state is already close to the solution.  Inspect
costs, residuals, and damping history before increasing the maximum gamma.

``lm_ga_history`` records the starting damping and the updated damping
after each outer iteration, including when ``display_progress=0``.
Unused trailing entries are NaN.  It does not record every rejected trial
step.  Non-LM methods return an empty history.  Inspect ``errors`` when
status 9 is returned; invalid inputs can also raise an exception before
iteration starts.

Convergence establishes a numerical stopping condition.  Also inspect the
state, structured measurement residuals, and prior departures.  The
reported costs are normalized by the number of measurements, not by
residual degrees of freedom, so the final measurement cost need not be one.
For correlated measurement errors, whiten residuals using a covariance
factorization, as described in :ref:`sec-oem-covariance`.  Dividing each
channel only by its standard deviation leaves the cross-channel
correlations in place.

With ``clear_matrices=0``, the gain and averaging kernel help distinguish
measurement information from prior constraints.  First call
:meth:`~pyarts3.workspace.Workspace.measurement_averaging_kernelCalc`.
The contributions to retrieval uncertainty can then be calculated with
:meth:`~pyarts3.workspace.Workspace.measurement_vec_error_covmat_observation_systemCalc`
and :meth:`~pyarts3.workspace.Workspace.model_state_covmat_smoothing_errorCalc`.
Their mathematical definitions and relation to posterior covariance are in
:ref:`sec-oem-uncertainty`.  The observation contribution alone is not the
full posterior covariance.  The `overview of uncertainty reporting in
atmospheric retrievals <https://amt.copernicus.org/articles/13/4393/2020/>`_
discusses how to report the role of prior information and smoothing in an
uncertainty budget.

Retrieval transformations
=========================

ARTS provides built-in retrieval transformations described in
:doc:`concept.oem`.  Custom
transformations can be assigned directly to any Jacobian target from Python by
providing its three operators:

.. code-block:: python

  target = ws.jac_targets.atm[-1]
  target.transform_state = lambda t, field: A @ (t - b)
  target.inverse_state = lambda x, field: A_inv @ x + b
  target.inverse_jacobian = lambda J, x, field: J @ A_inv

Here ``field`` is the complete owning field or data object, allowing mappings
that need information beyond the target itself.  Each callable must return a
vector or matrix with the same shape as the target block.  The example supports
a general invertible affine transformation; ``A`` need not be diagonal or
orthogonal.  The same interface can express bounded and other reversible
functional transformations.

See :doc:`concept.oem` for the forward/inverse transformation definitions and
the Jacobian chain rule used by these operators.
