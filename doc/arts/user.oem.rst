.. _sec-user-oem:

Configuring an optimal-estimation retrieval
###################################################

This guide covers method selection, settings, covariance setup, and retrieval
diagnostics.  The mathematical formulation is in :ref:`Sec OEM`.

Retrieval data and ownership
===========================

``ws.oem`` holds the observations, prior state, covariances, current state,
Jacobian, fitted measurements, basis information and retrieval results.
``oemCalc`` and ``oemCalcReduced`` update this object. Their method, iteration
and damping controls remain named call arguments. Normalization vectors are
members of ``ws.oem``; empty vectors disable normalization.

For target-based setup, ``oemInit`` starts an empty OEM object and resets
``jac_targets``. The ``oemAdd`` methods add targets and store their pending
covariance blocks in ``oem.covmat_diagonal_blocks``. Call
``jac_targetsFinalize`` when target offsets are needed for forward calculations
or ``model_state_vecFromData``. After setting the prior, observations and
measurement-error covariance, ``oemFinalizeDiagonal`` assembles
``oem.model_state_covmat`` and checks the complete setup. Success sets
``oem.checked``; missing or invalid inputs cause a detailed error and leave it
unchecked.
Use ``oemStateCovmatCorrelateConstant`` and ``oemMeasurementCovmatConstant`` to
construct or edit the owned covariances.

Prepare the prior with ``model_state_vecFromData``, then use ``oemSetApriori``
to consume ``model_state_vec``. ``oemSetMeasurement`` consumes the observations
from ``measurement_vec``. These methods only replace their respective inputs;
they retain covariances, bases and previous modeling results. A size change
makes the object unchecked; replacing values at the same size preserves
``checked``. If dimensions change, update the associated data and check again.

Forward calculations and state-mapping methods still use the primitive workspace
vectors and Jacobian. Retrieval results are read from ``ws.oem.model_state_vec``,
``ws.oem.measurement_vec_fit`` and ``ws.oem.measurement_jac``. A calculation does
not overwrite the independent ``ws.model_state_vec``, ``ws.measurement_vec_fit``
or ``ws.measurement_jac`` variables.

For a complete numerical input problem, ``oemInitFromData`` consumes
``model_state_vec`` as the prior, ``measurement_vec`` as the observations,
``model_state_covmat`` and ``measurement_vec_error_covmat``. It checks dimensions
before consumption and replaces the whole OEM object, leaving it unchecked.
It requires no physical fields and does not change the target mapping.
Consumed workspace variables are left empty.

Neither setup route imports ``measurement_vec_fit`` or ``measurement_jac``:
these are modeling results. OEM calculates them as needed. For information
analysis before retrieval, supply the Jacobian separately in
``oem.measurement_jac``.

``oemRestoreApriori`` restores only the retrieved physical quantities.
The OEM object, including the previous fitted state and diagnostics, is const
and remains unchanged. This supports preparing the next measurement from the
same physical prior while retaining the previous result.

Both covariance members can also be assigned explicitly.

``ws.oem.clear_auxiliary()`` (or ``oemClearAuxiliary``) releases fitted
measurements, Jacobian, gain, averaging kernel, error-covariance results,
basis spectrum, loss summaries, pending covariance blocks and prepared covariance caches. It preserves
the observations, prior, current state, covariances, selected bases,
normalization, diagnostics and ``checked``. It is also allowed on checked data.
Either calculation can run again afterwards;
the discarded products will be recomputed as needed. Basis selection needs
a spectrum, so rerun ``oemBasisCalc`` before selecting modes again.
``ws.oem.clear()`` instead resets the entire object, including its inputs.

Checking manually edited data
=============================

``oem.checked`` is read-only in Python. Call ``oem.check()`` to validate
numerical inputs and covariances. It reports all detected problems together,
leaving the object unchecked on failure. ``ws.oemCheck()`` also checks agreement
with the finalized target mapping and can be used inside an agenda.

After a successful check, replacing a Python attribute raises an error asking
for ``oem.uncheck()``. Unchecking preserves the data and permits replacement;
call ``check()`` again after editing. Reading members and changing values in
existing arrays or covariance objects remain allowed, including through NumPy
views. Such edits are not tracked automatically: recheck when needed.
``oemInit`` and XML loading produce unchecked data. Both ``oemCalc`` and
``oemCalcReduced`` validate unchecked data before preparing covariances or
changing outputs. Successful validation marks the object checked; failures
produce a detailed report. Repeated calculations reuse that status until an
operation invalidates it. Call ``ws.oemCheck()`` to force revalidation, including
after in-place edits.

Choosing a method
=========================

:meth:`~pyarts3.workspace.Workspace.oemCalc` always minimizes the same objective,
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
   * - ``li_m``
     - One step; direct solve in measurement space
     - Substantially fewer measurements than states; avoids CG tolerance tuning.
   * - ``li_cg``
     - One step; conjugate gradient (CG) in state space
     - A large linear problem where a direct solve is expensive.
   * - ``li_cg_m``
     - One step; CG in measurement space
     - A linear problem with substantially fewer measurements than states.
   * - ``gn``
     - Gauss--Newton; direct solve in state space
     - A mildly nonlinear problem with a plausible starting state.
   * - ``gn_m``
     - Gauss--Newton; direct solve in measurement space
     - Substantially fewer measurements than states; avoids CG tolerance tuning.
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
remove the prior.

Start development of a retrieval with a direct method on a small, representative
case.  Use ``li`` only after checking linearity; a single step applied to a
nonlinear model is only a local approximation.  Use ``gn`` when this
linearization is reliable and try ``lm`` when its unconstrained steps fail.
Check the Jacobian and units before tuning damping to compensate for poor
convergence.

The state-space system has size ``n`` by ``n``; the measurement-space
system has size ``m`` by ``m``.  Size is only a first guide to performance:
covariance structure, conditioning, and forward-model cost also matter.
CG uses ``cg_tolerance=1e-10`` by default. ``cg_max_iter=0`` selects a
budget of :math:`\max(1000,2d)` iterations per linear solve, where :math:`d`
is the dimension of the selected state- or measurement-space system. Set a
positive ``cg_max_iter`` for an explicit budget. ``stop_dx`` and ``max_iter``
control the outer retrieval iteration independently.
An unconverged CG step is rejected. GN and LI stop with
``LinearSolverLimit``; LM retries with more damping and uses that status if
the final solve also exhausts its budget. Numerical breakdown such as
non-positive curvature gives ``Error`` with an explanation in
``oem.diagnostics.errors``. Check covariance validity and scaling in that case.
OEM still stores the measurement Jacobian, and computing the gain matrix
requires additional dense matrices.  ``clear_matrices=1`` skips the gain
calculation and returns empty Jacobian and gain matrices when those outputs
are not needed.

Reducing forward-model work
-----------------------------------

OEM reuses simulations and Jacobians for exactly matching states within a
retrieval.  Continuing Gauss--Newton iterations obtain the simulation and
Jacobian together; LM trials need only the simulation until derivatives
are needed at an accepted state.  Derivative calculations use ``jac_targets``
and do nothing when it is empty.  State updates use the separate, complete
``model_state_targets``.  ``clear_matrices=1`` also avoids a final derivative
evaluation needed solely for retained matrix outputs.

OEM passes both target sets by reference to ``inversion_iterate_agenda``.
Custom agendas should use ``UpdateModelStates`` for state mapping and pass
``jac_targets`` to radiative transfer.  Measurement-error values use
``model_state_targets`` even when no derivatives are requested.
For standalone state updates outside OEM, supply the mapping explicitly:
``ws.UpdateModelStates(model_state_targets=ws.jac_targets)``.
Direct calls to the inversion agenda must supply ``model_state_targets`` and
use either the full ``jac_targets`` or an empty ``JacobianTargets``.

The forward model should be repeatable for the same state and configuration.
If supplying an initial ``measurement_vec_fit`` and ``measurement_jac``,
both must correspond to
the supplied starting state and current forward-model configuration.
Clear these outputs after changing that configuration to force reevaluation.

To inspect which workspace variables an agenda copies, use:

.. code-block:: python

   print(ws.inversion_iterate_agenda.document())
   # Inspect this helper too when it is used by the inversion agenda.
   print(ws.measurement_inversion_agenda.document())

The listing distinguishes shared variables from copied inputs that the
agenda modifies internally.  Python tuple-style operators can additionally
copy large fields when converting their arguments and returned values.
For a custom Python forward model, a ``CallbackOperator`` updating its
declared workspace outputs in place can avoid returning unchanged model
fields through that tuple interface.  Measure performance with a
representative retrieval before redesigning an agenda around copy costs.

Choosing covariance matrices
====================================

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

   ws.oemMeasurementCovmatConstant(value=sigma**2)

Supply observations with ``oemSetMeasurement`` first; the helper takes its
size from ``oem.measurement_vec``.

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
the same covariance, including its correlations.  Inverse blocks must cover
an entire group of coordinates connected by correlations when supplied
for that group.  Independent groups may leave their inverses uncomputed.
Supplying the covariance and allowing ARTS to compute its inverse avoids
having two potentially inconsistent descriptions.

Check a constructed covariance explicitly before using it:

.. code-block:: python

   ws.oem.model_state_covmat.validate(
       expected_size=len(ws.oem.model_state_vec_apriori)
   )
   ws.oem.measurement_vec_error_covmat.validate(
       expected_size=len(ws.oem.measurement_vec)
   )

Validation checks the represented covariance, including its block layout,
finite entries, symmetry, positive definiteness, and consistency with a
supplied inverse.  It raises an error when a check fails.  The default
``relative_tolerance=1e-10`` controls numerical comparisons; it does not
repair the matrix.  Successful validation establishes that the covariance
is numerically admissible.  The uncertainty values, correlations, units,
and ordering still need to match the intended physical problem.
The optional ``max_dense_elements=10_000_000`` bounds each dense connected
component needed for validation.  Independent diagonal errors are checked
without allocating a full dense covariance.

Coordinate changes and numerical scaling
------------------------------------------------

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
The state-sized normalization setting is unsupported for ``li_m``, ``gn_m``,
``li_cg_m`` and ``gn_cg_m``. These methods accept ``measurement_vec_normalization``: empty disables
scaling (the default); otherwise supply one finite positive :math:`D_{ii}` per
measurement. Noise standard deviations, :math:`D_{ii}=\sqrt{S_{\epsilon,ii}}`, are a useful
choice. The measurement-space system is scaled by their inverses. This
uses only vector scaling and preserves the statistical objective. For CG, the
relative residual tolerance applies to the scaled system. For correlated
errors this is diagonal scaling, not full whitening.

Inspect the noise standard deviations with::

    ws.oemMeasurementCovmatNormalization()
    ws.oemCheck()
    ws.oemCalc(method="gn_cg_m")

The output Vector contains :math:`D_{ii}` in measurement units. Its workspace name is
chosen by the caller. Computing this vector alone does not enable scaling;
pass it explicitly to OEM as shown above. Nonempty measurement scaling is
rejected for state-space methods.

.. _sec-user-oem-information:

Correlating temperature and log-water on the same grid
==============================================================

For two atmospheric retrieval targets with diagonal marginal covariances
and identical altitude, latitude and longitude grids, use:

.. code-block:: python

   ws.oemInit()
   ws.oemAddTemperature(
       matrix=np.diag(np.full(nlevels, 3.0**2)), d=1e-3)
   ws.oemAddSpeciesVMR(
       species="H2O", matrix=np.diag(np.full(nlevels, 0.2**2)), d=1e-7)
   ws.jac_targetsFinalize()
   ws.jac_targetsToggleLogarithmicAtmTarget(key="H2O")
   ws.model_state_vecFromData()
   ws.oemSetApriori()
   ws.oemSetMeasurement()
   ws.oemMeasurementCovmatConstant(value=measurement_variance)
   ws.oemFinalizeDiagonal()
   ws.oemStateCovmatCorrelateConstant(
       target1="temperature", target2="H2O", correlation=0.6)

Here ``measurement_vec`` must already contain the observations, and
``measurement_variance`` is their common error variance. The prior standard deviations are 3 K and 0.2 in natural-log water
VMR.  Positive correlation expresses a preference for warmer-than-prior
states to have more water than the prior at the same grid point.  It does
not impose a temperature-to-water conversion or alter either marginal
variance.  The supplied water variance must already describe log-water;
toggling the logarithmic target does not convert a VMR covariance.

For absorption models using numerical derivatives, choose ``d`` in physical
kelvin or VMR units, even for a logarithmic retrieval target.  The example
uses 0.001 K and 1e-7 VMR for its PWR98 model.  Check derivative stability
when adapting these perturbations to a different atmosphere.

The coefficient must be strictly between -1 and 1; perfect correlation
would make the covariance singular.  Calling the method again replaces
the pair's correlation; zero removes it.  Other pairs remain unchanged.
The complete covariance must remain positive definite, so adding another
pair can fail even when each coefficient individually lies in this range.
Failure leaves the covariance unchanged.  Dense and sparse diagonal
marginals are supported; targets must have one state coordinate per grid
point.  General spatial cross-covariances require explicit blocks.
The helper uses the covariance validator's default dense-component size
limit; a sparse cross block does not make validation or inversion fully sparse.

Checking what the measurements can constrain
====================================================

Use :func:`~pyarts3.retrieval.information` to examine a Jacobian together with
the assumed prior and measurement covariances.  It can run before a
retrieval, using a Jacobian already evaluated at a representative state.
It validates its inputs and produces a report without executing ``oemCalc``
or a forward-model agenda:

.. code-block:: python

   import numpy as np
   from pyarts3.retrieval import information

   # A dimensionless example: two quantities, two measurements.
   jacobian = np.diag([2.0, 0.1])
   prior_covariance = np.eye(2)
   measurement_covariance = np.eye(2)
   report = information(
       jacobian,
       prior_covariance,
       measurement_covariance,
       state_labels=["first quantity", "second quantity"],
   )
   print(report.describe())
   figure, axes = report.plot()

Covariance inputs can be square arrays, ARTS ``CovarianceMatrix`` objects,
or one-dimensional arrays of diagonal **variances**.  For example,
``measurement_covariance = np.ones(2)`` represents the same independent
unit errors as ``np.eye(2)`` above.  The diagonal representation avoids
allocating a full measurement covariance when there are many independent
channels.

Here the first quantity loses 80 percent of its prior variance, while
the second loses approximately 1 percent.  Converging the optimizer more
tightly cannot create sensitivity to the second quantity.  The example
has about 0.81 degrees of freedom for signal despite having two state
elements and two measurements.

For an existing workspace, use its current Jacobian and covariance values:

.. code-block:: python

   from pyarts3.retrieval import information_from_workspace

   report = information_from_workspace(ws)
   print(report)

This reads ``measurement_jac``, ``model_state_covmat``, and
``measurement_vec_error_covmat``.  It does not refresh the Jacobian or
change the workspace.  The Jacobian's columns must match the prior's
state coordinates and ordering, and its rows must match the measurement
covariance.  After a nonlinear retrieval, the report describes information
near the retrieved state; at the prior it describes the initial local
problem.  The two reports can differ because the sensitivities change.
Each report is a snapshot; create a new one after changing its inputs.

Read the report's quantities as follows:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Quantity
     - Interpretation
   * - ``singular_values``
     - Sensitivity of independent state patterns relative to prior and
       measurement uncertainty.  Values above one indicate stronger
       measurement constraints than prior constraints in that pattern.
   * - ``mode_variance_reduction``
     - Fraction of prior variance removed for each pattern, from zero
       to one.  Zero includes state directions the measurements cannot see.
   * - ``degrees_of_freedom``
     - Sum of those reductions: the effective number of state quantities
       constrained by the measurement.  It need not be an integer.
   * - ``information_bits``
     - Reduction of Gaussian uncertainty volume expressed as entropy
       in bits.  It is a different summary from degrees of freedom.
   * - ``prior_standard_deviation``, ``posterior_standard_deviation``
     - Uncertainties of individual state elements in their retrieval units.
       These assume the stated covariances and local Jacobian.
   * - ``variance_reduction``
     - Fraction of prior variance removed for each individual state element.
       Correlations make these different from the independent-mode values.
   * - ``prior_correlation_condition``, ``measurement_correlation_condition``
     - Condition numbers after removing coordinate units and variance scales.
       Large values identify nearly dependent error patterns; they do not
       measure the information content or replace covariance validation.

``state_modes[:, i]`` gives a state pattern in physical retrieval
coordinates, scaled to unit prior uncertainty.  A weak pattern may mix
several parameters, so inspecting Jacobian columns individually can miss
it.  ``measurement_modes`` describes patterns in whitened measurement
coordinates; with correlated errors these mix the original channels.
Signs of modes are arbitrary, and equally informative modes do not have
a unique orientation.  The spectrum and its connection to reduced
retrieval spaces are discussed by :cite:t:`nesser:21`; the expressions
used here are in :ref:`sec-oem-information`.

The plot shows variance reduction by mode and the ratio of posterior to
prior standard deviation by state element.  Ratios allow state elements
with different units to share an axis.  Use the numerical arrays and
``state_labels`` to make plots suitable for a particular profile or target.

For an optional check of the prior prediction, give the workspace helper
a prediction explicitly evaluated at the prior mean.  It then reads the
measurement from ``ws.oem.measurement_vec``:

.. code-block:: python

   report = information_from_workspace(
       ws,
       prior_prediction=prediction_at_prior,
   )
   print(report.innovation_chi_square)

``prediction_at_prior`` must contain ``F(xa)`` in measurement units, and
the Jacobian must be appropriate around that prior state for this check.
The helper cannot establish which state produced a stored simulation;
it never substitutes ``measurement_vec_fit`` automatically.  For the
linear Gaussian model with the stated uncertainties, this statistic has
mean equal to the number of measurements.  It is not a target value for
each individual realization or the fitted OEM measurement cost.  A large
value is a reason to inspect units, model biases, outliers, and uncertainty
assumptions together.
The array-based ``information`` function requires both ``measurement``
and ``prior_prediction`` explicitly when requesting this check.

These checks leave covariance choices explicit.  Compare scientifically
plausible uncertainty models or candidate measurement sets while keeping
track of what changed.  Increasing assumed prior uncertainty can increase
reported information without adding measurements.  A weak mode suggests
examining measurement coverage, state parameterization, or independent
prior knowledge; changing a covariance only to make the report look
better changes the question being answered.

The analysis uses dense linear algebra and includes all state modes,
including unobserved directions.  ``max_dense_elements`` defaults to
10,000,000 and bounds an estimate of dense analysis storage and the dense
covariance factors.  The estimate includes ``3*m*n + 4*n*n`` elements
for ``m`` measurements and ``n`` states.  Additional library work arrays
mean this is not an absolute bound on peak memory.  In particular,
retaining all state modes requires storage proportional to the square
of the state size.  Reduce the analysis size or raise the limit
deliberately when the guard rejects a large problem.

Setting LM damping
==========================

LM adds damping to limit the size of a Gauss--Newton step.  Larger gamma
penalizes larger steps, with the penalty scaled by the diagonal of the prior
precision matrix.  A zero value gives the Gauss--Newton step.  LM adjusts
damping by comparing actual and predicted cost changes, and may try several
forward-model evaluations within one outer iteration.  The damped system is
defined in :ref:`sec-oem-damping`.

Each outer iteration allows at most 100 LM linear solves, including an
additional undamped solve when needed to check stationarity.  This internal
limit is separate from ``max_iter`` and ``maximum_damping``.  Reaching the
trial limit, or failing to increase damping after a rejection because of
floating-point rounding, stops the retrieval with status ``Error`` and an
explanation in ``oem_diagnostics.errors``.  Increasing ``max_iter`` does not change the
trial limit.

Use :class:`~pyarts3.arts.LevenbergMarquardtSettings` to give the damping controls names.
For an already configured retrieval:

.. code-block:: python

   from pyarts3.arts import LevenbergMarquardtSettings

   damping = LevenbergMarquardtSettings(
       initial_damping=10.0,
       decrease_factor=2.0,
       increase_factor=2.0,
       maximum_damping=100.0,
       damping_threshold=1.0,
       convergence_damping_limit=0.0,
   )
   print(damping.describe())
   ws.oemCheck()
   ws.oemCalc(method="lm", max_iter=20, lm_ga_settings=damping)

These are the defaults of ``LevenbergMarquardtSettings()``.  They provide a visible
starting configuration to assess on representative retrievals.  Check
the forward model, Jacobian, and covariance assumptions before using
damping changes to address convergence problems.

The constructor accepts keyword arguments only, so each override states
what it changes.  ``print(damping)`` displays all six names and values;
``damping.describe()`` explains the configured behavior in words.
The same object works with ``lm``, ``ml``, ``lm_cg``, and ``ml_cg``.

The named controls are:

.. list-table::
   :header-rows: 1
   :widths: 8 30 62

   * - Setting name
     - Meaning
   * - ``initial_damping``
     - Initial damping, at least zero and no greater than the maximum.
   * - ``decrease_factor``
     - Divisor when damping is reduced; must be greater than one.
   * - ``increase_factor``
     - Multiplier when damping is increased; must be greater than one.
   * - ``maximum_damping``
     - Positive maximum damping.  Failure to find an acceptable step at
       this value stops the retrieval.
   * - ``damping_threshold``
     - Positive restart value when a step with damping below this value
       fails.  When a proposed decrease would fall below this value,
       damping becomes zero.
       Must not exceed the maximum.
   * - ``convergence_damping_limit``
     - Nonnegative upper damping limit for enabling the ordinary
       ``stop_dx`` criterion.  This refers to the current damping after its
       update, not the lowest value used in an earlier accepted iteration.

All entries must be finite.  ``initial_damping`` and
``convergence_damping_limit`` may be zero; the maximum and threshold must
be positive.  Both factors must be greater than one.  The initial damping
and threshold must each be no greater than the maximum.

Construction and each field edit validate all six settings immediately.
An invalid edit raises an error naming the affected setting and leaves
the object unchanged.  You can also call ``validate()`` explicitly;
``oemCalc`` checks the values again:

.. code-block:: python

   damping.initial_damping = 20.0
   damping.validate()

For related changes, construct a replacement with the desired keyword
arguments together.  When editing fields individually, keep each
intermediate configuration valid.  For example, raise ``maximum_damping``
before setting ``initial_damping`` above the old maximum:

.. code-block:: python

   damping.maximum_damping = 1000.0
   damping.initial_damping = 200.0

How the controls interact
---------------------------------

``decrease_factor`` is a **divisor**: a value of 2 halves the damping when
a reduction is made.  Increasing this factor removes damping faster.
``increase_factor`` is a multiplier: a value of 2 doubles the damping
after rejection, up to the maximum.  Increasing this factor tries more
strongly damped steps sooner.  A factor of 0.5 is invalid for either setting.

``damping_threshold`` has two roles.  A proposed reduction below it sets
damping to zero, allowing a Gauss--Newton step.  If a step with damping
below the threshold fails, the next trial restarts at the threshold.
It is not a minimum damping: an initial value below it is allowed.
For example, with a threshold of 1 and a decrease factor of 2, successive
reductions from 10 give 5, 2.5, 1.25, and then 0.  Not every accepted step
causes a reduction.  The local model must predict the cost change well
enough, and a step accepted after a rejection keeps the damping used for
its final trial.

``convergence_damping_limit`` gates the ordinary ``stop_dx`` check.  It
does not set the accuracy of that check.  The default of zero enables the
check once damping reaches zero.  A positive limit permits termination
while damping still constrains the steps; a small damped step can then
hide a remaining distance to the minimum.  The gate uses the updated
damping, so a step calculated with damping 10 and then reduced to 5 can
pass a limit of 5.  The iteration history records this updated value.
LM can also recognize numerical stationarity using an undamped step,
independently of this gate.  A step made tiny only by strong damping does
not establish stationarity.

Use the behavior of representative retrievals to guide changes:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Observation
     - What to examine or change
   * - Early trial states are too aggressive.
     - Check the Jacobian and state parameterization, then try increasing
       ``initial_damping``.  Damping cannot enforce physical bounds.
   * - Convergence is reported while damping remains large.
     - Check ``convergence_damping_limit``.  A value of zero requires
       damping to be removed before the ordinary stopping check is enabled.
   * - Useful progress continues at the iteration limit.
     - Increase ``max_iter`` and compare the final state and costs.
       This increases the outer iteration budget; each LM iteration can
       contain several trial evaluations.
   * - Damping repeatedly reaches its maximum.
     - Inspect the forward model, Jacobian, state parameterization, and
       covariance scales.  If the state is already near the solution,
       ``stop_dx`` may demand changes below numerical accuracy.
       Raising the maximum alone does not resolve these causes.

Damping controls how the solver approaches a minimum.  The covariance
matrices define the statistical problem being solved.  Choose them from
the prior and measurement error models; changing them to cure a convergence
problem also changes the inferred state and uncertainty.

Storing settings
---------------------------------

``LevenbergMarquardtSettings`` is a workspace group with named printing and
XML storage. Pass it directly as ``lm_ga_settings``. OEM uses the default
``LevenbergMarquardtSettings()`` when the argument is omitted. The
six-element input shorthand remains supported in Python:

.. code-block:: python

   damping = LevenbergMarquardtSettings([10, 2, 2, 100, 1, 0])
   ws.oemCheck()
   ws.oemCalc(method="lm", lm_ga_settings=[10, 2, 2, 100, 1, 0])

The order is ``initial_damping``, ``decrease_factor``, ``increase_factor``,
``maximum_damping``, ``damping_threshold``, ``convergence_damping_limit``.
Exactly six values are required and validated. This is a Python-only input
conversion; neither settings nor diagnostics provide ``as_vector()``.

Checking the result
===========================

``oem_diagnostics`` is an ``OptimalEstimationDiagnostics`` object. Access
fields by name, for example ``ws.oem.diagnostics.status``.
It contains the following fields:

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``status``
     - ``OptimalEstimationStatus``: ``NotRun``, ``Converged``,
       ``IterationLimit``, ``DampingLimit``, ``Error``, or ``StartCostLimit``.
       ``Converged`` also includes LM numerical stationarity.
   * - ``initial_cost``
     - Starting total cost, divided by the number of measurements.
   * - ``final_cost``
     - Final total cost, divided by the number of measurements.
   * - ``measurement_cost``
     - Final measurement contribution to the cost, divided by the number of measurements.
   * - ``iterations``
     - Integer number of completed outer iterations, initially zero.
   * - ``lm_ga_history``
     - Starting and updated LM damping values.
   * - ``errors``
     - Errors and warnings recorded by OEM.

Unavailable costs are NaN.  ``max_start_cost`` is a limit on the total
cost at the starting state; its default is infinity.  A value at or below
zero disables this limit and can skip computation of the starting cost
for non-LM methods when progress output is off.  It is not an upper bound
on the final residual.  The ``li`` variants perform exactly one step and
can return status ``IterationLimit`` even for an exact linear solution, because they do not
take a second step to establish iterative convergence.

``stop_dx`` defaults to 0.01 and controls the weighted state-step convergence
measure described in :ref:`sec-oem-convergence`.  It does not set a relative
change in cost or an unweighted state difference.  Lower it for a stricter
convergence requirement and increase ``max_iter`` if useful progress continues
at the iteration limit.

Near a solution, differences between computed costs can be too small to
resolve reliably.  LM then checks an undamped step against ``stop_dx``
and the cost's floating-point resolution.  If both checks establish
stationarity, it returns status ``Converged`` without exhausting the damping range.
An exactly zero gradient also establishes stationarity, including when
the minimum cost is nonzero.  Very large damping alone cannot pass these
checks: for example, ``initial_damping=maximum_damping=1e20`` can return
status ``DampingLimit`` with an unchanged state when no acceptable step is found.

An excessively small ``stop_dx`` can still demand accuracy beyond the
forward model or linear solve.  These stationarity checks do not measure
forward-model noise or establish a global minimum.  Inspect costs,
residuals, and damping history before increasing the maximum gamma.

``oem_diagnostics.lm_ga_history`` records the starting damping and the updated damping
after each outer iteration, including when ``display_progress=0``.
Unused trailing entries are NaN.  It does not record every rejected trial
step.  Non-LM methods return an empty history.  Inspect ``oem_diagnostics.errors`` when
status ``Error`` is returned; invalid inputs can also raise an exception before
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
:meth:`~pyarts3.workspace.Workspace.oemAveragingKernelCalc`.
The contributions to retrieval uncertainty can then be calculated with
:meth:`~pyarts3.workspace.Workspace.oemObservationErrorCalc`
and :meth:`~pyarts3.workspace.Workspace.oemSmoothingErrorCalc`.
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

Direct measurement-space solvers
----------------------------------------

``li_m`` and ``gn_m`` assemble and factor an :math:`m\times m` matrix, where :math:`m` is the
number of measurements. They use the same measurement-space update as the
CG variants, with one step for ``li_m`` and iteration for ``gn_m``. Optional
``measurement_vec_normalization`` scales the direct system as well.
With measurement Jacobian :math:`\mathbf{J}`, assembly computes
:math:`\mathbf{B}=\mathbf{S}_a\mathbf{J}^{\top}` once, forms
:math:`\mathbf{M}=\mathbf{J}\mathbf{B}+\mathbf{S}_\epsilon` with matrix-matrix
operations, and reuses :math:`\mathbf{B}` for the state update. This avoids a state-sized
normal matrix but uses additional :math:`n\times m` temporary storage. The CG variants
retain lazy assembly and are the lower-memory alternative.
Small :math:`m` relative to the state dimension is therefore a useful starting point
for method choice, not a guarantee of improved runtime. Requested gain-matrix
output still uses the existing state-space postprocessing, so the complete
retrieval can retain state-sized costs.

Covariance storage
-------------------------------------------------

Covariance preparation is automatic when calling ``oemCalc``; it does not require
an additional workspace call or a factorization setting. Continue choosing
``method="lm"`` or ``method="lm_cg"`` explicitly. Matrix and Sparse are input
storage choices, not solver choices: a diagonal Matrix can use diagonal
solves, while a correlated Sparse component without inverse blocks can use dense Cholesky.
If any inverse blocks are supplied, the current preparation selects the
inverse-application path for the whole covariance and fills missing component
inverses. It does not yet mix supplied inverses and Cholesky factors within
one covariance. In particular,
``oemMeasurementCovmatConstant`` supplies both a sparse diagonal
covariance and its inverse. LM currently requests explicit prior precision;
the measurement covariance is where diagonal/factorized solves are most
useful.

Without ``ARTS_HEADLESS``, the script plots the prior, manipulated starting
state, fitted states and fitted measurements. Four overlapping marker shapes
show the equivalent covariance representations for each correlation pattern.

Identifying the state coordinates in an information report
------------------------------------------------------------------

``pyarts3.retrieval.information_from_workspace(ws)`` uses ``jac_targets`` to
label the state coordinates, for example ``atm.H2O[0]``. It uses each target's
``x_start`` and ``x_size``, not the order of target categories. The printed
report also lists field slices such as ``x[0:3]: atm.H2O (3 entries)``.
``report.state_blocks`` contains immutable ``(name, start, size)`` tuples;
``report.state_labels`` contains one label per Jacobian column. Small-state
uncertainty plots use these labels automatically.

The index in a field label is its flattened target index, not an altitude.
Keys alone do not specify units or logarithmic/relative state transformations.
Supply ``state_labels=[...]`` when more specific per-coordinate labels are
needed. Missing metadata retains generic ``x[i]`` labels; target ranges outside
the Jacobian columns are rejected rather than attaching misleading names.
The function does not finalize targets, run an agenda or modify the workspace.

ReducedOEM
-----------

``oemCalcReduced`` runs the forward model at full size and reduces the state and
measurement coordinates used by the solver. Both matrices must be supplied:
``model_state_basis_mat`` has full-state rows and reduced-state columns;
``measurement_basis_mat`` has reduced-measurement rows and full-measurement
columns. Their columns and rows, respectively, must be linearly independent.
Their normalization is arbitrary: ReducedOEM transforms both covariances.
They can be stored as workspace variables or passed explicitly.

Given a current full Jacobian and the covariances, prepare the full bases,
then choose which modes to retain::

    # Uses measurement_jac, model_state_covmat and measurement_vec_error_covmat
    ws.oemBasisCalc()
    print(ws.oem.basis_singular_values)
    ws.oemBasisReduce(max_lost_information_bits=0.01)
    print(ws.oem.basis_lost_dofs, ws.oem.basis_lost_information_bits)
    ws.oem.model_state_vec = []  # Start at the prior
    ws.oemCheck()
    ws.oemCalcReduced(method="lm")

``oemBasisCalc`` sets ``model_state_basis_mat``,
``measurement_basis_mat`` and ``oem_basis_singular_values``. It requires
only the Jacobian and covariances; it does not run a forward agenda or need
a measurement vector. With ``full_matrices=1`` (the default), both full bases are square and include all null-space
directions. This step only changes coordinates and loses no information.
The state basis reconstructs the full prior covariance when multiplied by
its transpose. For the measurement covariance, it is the inverse basis
times its inverse transpose that reconstructs the original covariance.

Use ``oemBasisCalc(full_matrices=0)`` to keep only :math:`\min(m,n)` modes
in each basis. This avoids a dense measurement-square basis when there are
many more measurements than states. It omits extra null-space directions of
the supplied Jacobian. If there are more states than measurements, the state
basis no longer reconstructs the entire prior covariance: omitted directions
retain prior uncertainty. Nonlinear sensitivity can change, so assess this
reduction at representative states.

``oemBasisReduce`` retains leading columns and rows in place in
``model_state_basis_mat`` and ``measurement_basis_mat``. The example discards
at most 0.01 bits of total local Gaussian information, choosing the smallest
retained rank that meets the budget. Both transformed covariances are identity
matrices. The two loss outputs report actual discarded DOFS and bits, also
when choosing an explicit rank.

Selection overwrites both basis matrices, while preserving the full spectrum
for reporting total information loss. Further calls can remove more modes.
To restore removed modes, restore saved copies of both matrices or rerun
``oemBasisCalc``. Keep the spectrum and both bases from the same
decomposition together.

Calling ``ws.oemBasisReduce()`` removes only modes with zero computed
information. Tiny nonzero values caused by roundoff require a positive
budget. Individual zero entries in a covariance or Jacobian do not identify
redundant directions: correlations can mix coordinates, and even a matrix
without zero entries can have a null space.

Use ``max_lost_dofs`` to limit the total discarded DOFS instead, or supply
both limits to require both. These are absolute totals over all discarded
modes, not per-mode thresholds or percentages. A value of -1 leaves a
limit unset. With both unset, the information-bit budget is zero.

An explicit ``rank=r`` retains exactly ``r`` leading state modes, from 1
through the full state size. It cannot be combined with loss limits. The
default ``rank=-1`` enables automatic selection. At least one coefficient
is retained even when all modes are uninformative, because ``oemCalcReduced``
requires a nonempty state. The measurement basis retains the smaller of
the selected rank and the measurement count. Inspect
``ws.oem.model_state_basis_mat.shape[1]`` for the selected rank.

For no reduction of either dimension, call ``oemCalcReduced`` directly after
``oemBasisCalc``, without calling ``oemBasisReduce``.

Selecting the full state rank through ``oemBasisReduce`` can still
remove measurement null-space rows when there are more measurements than
states. The matrices before selection also retain those rows.

The Jacobian must use the same state coordinates and measurement units as
the covariances. Bases are computed at that Jacobian's linearization point
and stay fixed during the retrieval. Recompute them after changing the
covariances or choosing a different linearization point. Basis signs are
arbitrary; compare the represented subspaces or retrieval results rather
than individual signs.

Diagonal measurement covariances remain inexpensive to factor. Correlated
components are factored separately. Saving complete bases requires dense
state-square and measurement-square matrices, as well as the whitened
Jacobian and SVD working storage. There is no configured allocation cutoff.
Selection replaces the basis matrices with their retained columns and rows,
reducing subsequent OEM storage and work. Retain separate copies yourself if
you need to restore removed modes without repeating the decomposition.

The Python information report also generates both bases and describes the
local information retained and discarded::

    report = pyarts3.retrieval.information_from_workspace(ws)
    reduction = report.reduction(rank=r)
    ws.oem.model_state_vec = []  # Start at the prior
    ws.oem.model_state_basis_mat = reduction.model_state_basis_mat
    ws.oem.measurement_basis_mat = reduction.measurement_basis_mat
    ws.oemCheck()
    ws.oemCalcReduced(method="lm")

Alternatively select the smallest rank meeting an information-loss budget at
the report's linearization point::

    reduction = report.reduction(max_lost_dofs=0.05, max_lost_information_bits=0.01)
    print(reduction)  # Rank and discarded DOFS/bits

Assign ``reduction.model_state_basis_mat`` and
``reduction.measurement_basis_mat`` to the corresponding members of ``ws.oem``. With arrays already
available outside a workspace, use
``pyarts3.retrieval.information(J, Sa, Se).reduction(rank=r)``.

These budgets are absolute DOFS and bits, not percentages. With two budgets,
both must be satisfied. At least one state column and one measurement row
are retained even when all modes are uninformative. Zero loss allows
discarding only zero-information modes in the computed spectrum.

DOFS is a sum of fractional mode contributions, not a count of modes. For
example, 100 modes each contributing 0.01 DOFS have total DOFS 1, but retaining
90 percent requires 90 modes. Avoid choosing rank by rounding total DOFS.

The helper's measurement reduction combines and noise-whitens channels.
To reduce only the state, supply an identity measurement matrix instead::

    ws.oem.model_state_basis_mat = reduction.model_state_basis_mat
    ws.oem.measurement_basis_mat = np.eye(len(ws.oem.measurement_vec))
    ws.oemCheck()
    ws.oemCalcReduced(method="lm")

Conversely, an identity state matrix leaves the state dimension unchanged.
Rows of an identity measurement matrix can select physical channels, but
channel selection can lose information that weighted combinations retain.
The measurement covariance is transformed with the same matrix, including
correlations between channels.

Measurement reduction forms signed weighted combinations of channels;
it does not crop the frequency grid. ``oemCalcReduced`` still
computes the full forward spectrum and Jacobian on every required agenda
call, then projects them. Its measurement reduction saves solver work;
skipping far-wing radiative-transfer calculations would require a separate
change to the sensor or frequency grid. Full measurement outputs remain
available after the retrieval.

All OEM method choices remain available. An explicit starting state must lie
in the affine subspace through the prior. Optional normalization vectors
refer to the reduced dimensions; the usual OEM method restrictions apply.
The reduction matrices remain fixed during the call. A reduction chosen at
one state can become unsuitable for a strongly nonlinear problem.

State, fit, Jacobian and gain outputs retain their full dimensions.
The physical agenda uses the real, full model-state targets. Full Jacobians
are computed before projection, so this saves solver work but does not reduce
Jacobian construction or storage. Input fit/Jacobian caches, if supplied,
must describe the starting state, as for OEM.

Diagnostic initial/final costs and ``max_start_cost`` use the original
measurements, covariances and measurement count, including discarded
residuals. Convergence and iteration progress use the reduced objective.
Consequently, a converged reduced retrieval can still have a large full
measurement cost.

``reduction.posterior_covariance()`` reconstructs the full-state covariance
of the truncated local linear model, preserving discarded prior uncertainty.
It allocates a full square matrix and uses the report's original Jacobian,
not a subsequently retrieved nonlinear state. For the helper's matched
singular-mode reductions in a linear problem this agrees with the sum of the
workspace smoothing-error and observation-error covariance contributions.
For arbitrary reductions that error sum can also include coupling from
discarded modes; see :ref:`sec-reduced-oem`.

Measurement-only grouping
^^^^^^^^^^^^^^^^^^^^^^^^^

``oemMeasurementBasisCalc`` creates only ``measurement_basis_mat``, using
the current Jacobian and measurement error covariance. It automatically groups
channels whose complete Jacobian rows are proportional, including opposite
signs. Matching one parameter's derivative is insufficient when other
parameters are retrieved. Matching is exact after row normalization; similar
but distinct sensitivities remain separate.

For diagonal measurement noise the result uses sparse storage, combining
channels according to their sensitivities and noise variances. Repeated
observations still contribute their extra precision. With correlated noise,
the covariance solve may produce a dense projection. Use an identity
``model_state_basis_mat`` with ``oemCalcReduced`` to keep the full state.

Both ``model_state_basis_mat`` and ``measurement_basis_mat`` hold a
``BlockMatrix``, accepting dense ``Matrix`` or ``Sparse`` input. The
``is_sparse`` property reports the storage type and
``shape`` inspects dimensions without converting the data. For a sparse result,
``matrix.tocsr()`` exposes SciPy sparse data; conversion to a NumPy array
explicitly allocates dense storage.

These groups preserve the supplied linear model's retrieval information.
They can retain more measurements than the SVD method because general linear
dependencies between different sensitivity directions are not combined.
An all-zero group is retained as one measurement. No state modes, singular
spectrum or loss outputs are generated; ``oemBasisReduce`` applies to
the SVD outputs, not to these groups. Recalculate after changing the noise
covariance or the Jacobian used for grouping. A fixed grouping remains a
local approximation for nonlinear retrievals.

Sparse state bases are useful for selecting state components or interpolating
from a smaller state grid. An exact identity state basis reuses the original
prior covariance. Sparse bases with diagonal priors also avoid dense state
projection products; a general correlated prior can still require dense
preparation and a dense reduced covariance.
