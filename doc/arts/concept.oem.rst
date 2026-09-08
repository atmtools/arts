.. _Sec OEM:

Optimal estimation
###################

For configuration and interpretation of retrieval outputs, see
:ref:`sec-user-oem`.

The core expression of optimal estimation in ARTS
follows from :cite:t:`rodgers:00`.
He formulates the core expression of a measurement as

.. math::

  \vec{y} = F\left(\vec{x}\right) + \epsilon,

where :math:`\vec{y}` is a measurement vector
(i.e., :attr:`~pyarts3.workspace.Workspace.measurement_vec`
in the ARTS workspace),
:math:`\vec{x}` is the state of the model
(i.e., :attr:`~pyarts3.workspace.Workspace.model_state_vec`
in the ARTS workspace),
:math:`F` is the model (i.e., ARTS itself), and
:math:`\epsilon` is some measurement error that cannot
be modelled.

The forward model result from this is just

.. math::

  \vec{y}_f = F\left(\vec{x}\right).

The goal of the optimal estimation is to find an :math:`\vec{x}`
that minimizes some cost-function, which in ARTS is defined as

.. math::
  \chi^2 = \chi^2_y + \chi^2_x,

where

.. math::
  \chi^2_y = \frac{1}{m} \left(\vec{y}-\vec{y}_f\right)^\top \mathbf{S}_\epsilon^{-1} \left(\vec{y}-\vec{y}_f\right).

  \chi^2_x = \frac{1}{m} \left(\vec{x}-\vec{x}_a\right)^\top \mathbf{S}_a^{-1} \left(\vec{x}-\vec{x}_a\right),

where
:math:`m` is the number of measurements,
:math:`\vec{x}_a` is an a priori of :math:`\vec{x}`,
:math:`\mathbf{S}_\epsilon` is the covariance of the measurement noise, and
:math:`\mathbf{S}_a` is the covariance of the model state a priori.

Linearization
=============

The iterative methods use a local linear approximation to the forward model
around the current state :math:`\vec{x}_i`:

.. math::

  F(\vec{x}) \approx F(\vec{x}_i) + \mathbf{J}_i (\vec{x} - \vec{x}_i).

Here

.. math::

  \mathbf{J}_i = \left.\frac{\partial F}{\partial \vec{x}}\right|_{\vec{x}_i}

is the Jacobian matrix
(i.e., :attr:`~pyarts3.workspace.Workspace.measurement_jac`
in the ARTS workspace).

One approach to minimize the cost function is to
use a Gauss-Newton approach to update the state
of the atmosphere.  This might look like

.. math::
  \vec{x}_{i+1} = \vec{x}_a + \mathbf{S}_a
  \mathbf{J}^\top\left(\mathbf{J}\mathbf{S}_a \mathbf{J} ^\top+\mathbf{S}_\epsilon\right)^{-1}
  \left[\vec{y}-\vec{y}_f+
  \mathbf{J}\left(\vec{x}_i-\vec{x}_a\right)\right],

.. math::

  \vec{y}_f = F\left(\vec{x}_i\right)

The Jacobian and simulated measurement in this expression are evaluated at
the current state.  Iteration stops according to the convergence criterion
and iteration limit.  The state-step measures are defined below;
settings are described in :ref:`sec-user-oem`.  A decrease in cost alone does not
establish convergence.

.. _sec-oem-covariance:

Covariances and coordinates
===========================

A covariance matrix has variances on its diagonal.  Its off-diagonal entries
can be expressed in terms of standard deviations and a correlation matrix:

.. math::

   S_{ij}=\sigma_i R_{ij}\sigma_j.

Thus each covariance entry has the product of the units of its two coordinates.
The covariance is symmetric and positive semidefinite; the inverses in the
optimal-estimation expressions require positive definiteness.  A positive
diagonal alone does not imply positive definiteness.  The inverse covariance
is the precision matrix.  In general, correlations imply

.. math::

   (\mathbf{S}^{-1})_{ii} \ne \frac{1}{S_{ii}}.

For a coordinate transformation :math:`\vec{x}=f(\vec{t})`, the local
covariance transformation is

.. math::

   \mathbf{S}_x \approx \mathbf{B}\mathbf{S}_t\mathbf{B}^{\top},
   \qquad \mathbf{B}=\frac{\partial f}{\partial\vec{t}}.

This is exact for an affine transformation.  A nonlinear transformation
also changes the shape of the distribution; a Gaussian prior in native
coordinates need not remain Gaussian in transformed coordinates.

For correlated measurement errors, the quadratic measurement cost can be
written using a whitened residual.  With a Cholesky factorization
:math:`\mathbf{S}_\epsilon=\mathbf{L}\mathbf{L}^{\top}`,

.. math::

   \vec{r}_w=\mathbf{L}^{-1}\bigl(\vec{y}-F(\vec{x})\bigr),
   \qquad \chi_y^2=\frac{\vec{r}_w^{\top}\vec{r}_w}{m}.

Numerical state scaling is distinct from a change of statistical coordinates.
For a diagonal matrix of positive scales :math:`\mathbf{T}`, solving
:math:`\mathbf{H}\Delta\vec{x}=-\vec{g}` is equivalent in exact arithmetic to

.. math::

   (\mathbf{T}\mathbf{H}\mathbf{T})\vec{z}=-\mathbf{T}\vec{g},
   \qquad \Delta\vec{x}=\mathbf{T}\vec{z}.

This changes the conditioning of the linear system while retaining the
same objective and solution in the original coordinates.

.. _sec-oem-damping:

Levenberg--Marquardt damping
============================

Define the half-gradient of the unnormalized objective and its
Gauss--Newton approximation to the half-Hessian by

.. math::

   \vec{g}=\mathbf{J}^{\top}\mathbf{S}_\epsilon^{-1}
     \bigl(F(\vec{x})-\vec{y}\bigr)
     +\mathbf{S}_a^{-1}(\vec{x}-\vec{x}_a),
   \qquad
   \mathbf{H}=\mathbf{J}^{\top}\mathbf{S}_\epsilon^{-1}\mathbf{J}
     +\mathbf{S}_a^{-1}.

The damped step used by ARTS is

.. math::

   (\mathbf{H}+\gamma\mathbf{D})\Delta\vec{x}=-\vec{g},
   \qquad
   \mathbf{D}=\operatorname{diag}
               \bigl(\operatorname{diag}(\mathbf{S}_a^{-1})\bigr).

Larger :math:`\gamma` penalizes larger steps in the precision-scaled
coordinates.  At :math:`\gamma=0`, the step equals the Gauss--Newton step.
The damping penalty modifies the local step, not the objective being
minimized.

.. _sec-oem-convergence:

State-step convergence measures
================================

For :math:`n` retrieved state elements, the state-space formulation uses
the Rodgers 5.31 measure

.. math::

   d_{531}=\frac{|\Delta\vec{x}^{\top}\vec{g}|}{n},

where the half-gradient :math:`\vec{g}` is evaluated before the step.
The measurement-space formulation uses the Rodgers 5.30 measure

.. math::

   d_{530}=\frac{\Delta\vec{x}^{\top}\mathbf{H}\Delta\vec{x}}{n}.

Both quantities are dimensionless.  The measures agree for an exact
undamped Gauss--Newton step, but need not agree for damped or approximate
steps.  They are distinct from relative changes in the cost function and
from an unweighted distance between state vectors.

.. _sec-oem-uncertainty:

Gain, averaging kernel, and retrieval uncertainty
=================================================

In the linear Gaussian model, the posterior covariance, gain, and
averaging kernel are

.. math::

   \widehat{\mathbf{S}}=
     \bigl(\mathbf{S}_a^{-1}+\mathbf{J}^{\top}
     \mathbf{S}_\epsilon^{-1}\mathbf{J}\bigr)^{-1},
   \qquad
   \mathbf{G}=\widehat{\mathbf{S}}\mathbf{J}^{\top}\mathbf{S}_\epsilon^{-1},
   \qquad
   \mathbf{A}=\mathbf{G}\mathbf{J}.

The observation and smoothing contributions in state coordinates are

.. math::

   \mathbf{S}_{\rm obs}=\mathbf{G}\mathbf{S}_\epsilon\mathbf{G}^{\top},
   \qquad
   \mathbf{S}_{\rm smooth}=(\mathbf{I}-\mathbf{A})\mathbf{S}_a
                           (\mathbf{I}-\mathbf{A})^{\top}.

With the stated covariances and model assumptions,
:math:`\widehat{\mathbf{S}}=\mathbf{S}_{\rm obs}+\mathbf{S}_{\rm smooth}`.
For a nonlinear retrieval, evaluating these expressions at the retrieved
state gives a local approximation.  Their uncertainty interpretation
requires that the prior covariance represents the assumed prior uncertainty.
The observation contribution alone is not the full posterior covariance.

.. _sec-oem-information:

Information carried by the measurements
=======================================

Factor the prior and measurement covariances as
:math:`\mathbf{S}_a=\mathbf{L}_a\mathbf{L}_a^{\top}` and
:math:`\mathbf{S}_\epsilon=\mathbf{L}_\epsilon\mathbf{L}_\epsilon^{\top}`.
The dimensionless sensitivity matrix and its singular value decomposition are

.. math::

   \widetilde{\mathbf{J}}=
     \mathbf{L}_\epsilon^{-1}\mathbf{J}\mathbf{L}_a
     =\mathbf{U}\mathbf{\Sigma}\mathbf{V}^{\top}.

Each singular value :math:`s_i` measures the response of a unit-prior
state mode relative to the measurement error.  Missing singular values
when :math:`m<n` are zero.  In coordinates
:math:`\vec{z}=\mathbf{L}_a^{-1}(\vec{x}-\vec{x}_a)`, the columns of
:math:`\mathbf{V}` describe independent prior modes.  Their posterior
variances are :math:`1/(1+s_i^2)`, giving the mode variance reductions

.. math::

   a_i=\frac{s_i^2}{1+s_i^2}.

These are the eigenvalues of the averaging kernel expressed in prior
coordinates.  Their use as an information spectrum is described by
:cite:t:`nesser:21`.  A mode with :math:`s_i=1` loses half its prior
variance; an unobserved mode with :math:`s_i=0` retains its prior variance.

The physical state modes are the columns of
:math:`\mathbf{P}=\mathbf{L}_a\mathbf{V}`.  They satisfy
:math:`\mathbf{P}^{\top}\mathbf{S}_a^{-1}\mathbf{P}=\mathbf{I}`.
The posterior covariance can be reconstructed as

.. math::

   \widehat{\mathbf{S}}=
     \mathbf{P}\operatorname{diag}\left(\frac{1}{1+s_i^2}\right)
     \mathbf{P}^{\top}.

The marginal variance reduction for state element :math:`j` is
:math:`1-\widehat{S}_{jj}/S_{a,jj}`.  This refers to an individual
coordinate and is generally different from the variance reduction of
a mode combining several state elements.  A mode's sign is arbitrary;
equal singular values also allow rotations within their shared subspace.

Two scalar summaries are the degrees of freedom for signal and the
Gaussian entropy reduction in bits:

.. math::

   d_s=\operatorname{tr}(\mathbf{A})=\sum_{i=1}^{n}a_i,
   \qquad
   H=\frac{1}{2}\log_2\frac{\det\mathbf{S}_a}{\det\widehat{\mathbf{S}}}
    =\frac{1}{2}\sum_{i=1}^{n}\log_2(1+s_i^2).

These quantities depend on the Jacobian and assumed covariances, not on
the realized measurement residual.  For a nonlinear forward model they
describe only the local linear approximation at the Jacobian's state.

For a measurement and a forward prediction at the prior mean, define
the innovation and its covariance by

.. math::

   \vec{r}=\vec{y}-F(\vec{x}_a),
   \qquad
   \mathbf{S}_r=\mathbf{S}_\epsilon+
                   \mathbf{J}\mathbf{S}_a\mathbf{J}^{\top}.

Under the linear Gaussian model with independent prior and observation
errors, :math:`\vec{r}^{\top}\mathbf{S}_r^{-1}\vec{r}` follows a
chi-squared distribution with :math:`m` degrees of freedom and mean
:math:`m`.  This is a check of a prior prediction and the assumed
uncertainties.  Replacing the prior prediction with a fitted measurement
does not give this distribution.

Transforming the Jacobian matrix
================================

It is sometimes desired to transform the Jacobian matrix away from the
native units that are available in the model.  Instead of retrieving
volume-mixing-ratio, it is sometimes more stable to retrieve the
relative change in the volume-mixing-ratio, or perhaps even
the logarithm of the volume-mixing-ratio, to avoid negative
values creeping in to what is otherwise an intensive and numerical
exercise.

There are three ways perform transformations on the Jacobian
matrix in ARTS.
You can change the measurement unit (e.g., from spectral radiance
to brightness temperature), you can change the model state unit
(e.g., from volume-mixing-ratio to relative humidity), or you can
map the Jacobian matrix to another vector space (e.g., from
Cartesian coordinates to spherical coordinates).

The first type of transformations are performed by an
:class:`~pyarts3.arts.SpectralRadianceTransformOperator`
in ARTS.  This is a local operation that is almost trivial to undo.
Generally, you should just use one of the
provided enumeration values of
:class:`~pyarts3.arts.SpectralRadianceUnitType` to
set up these types of transformations.
The enumeration class also describes how the transformations
are done.  This will not be repeated here.

The second and third types of transformations are a bit more complicated
because they transform the model state vector.
The optimal estimation methods are sometimes iterative, and they must be
able to update the state of the model between iterations.
These two types of transformations are therefore required to
be able to transform the model state vector from native to non-native
units, and vice versa, as well as to transform the Jacobian matrix
from native to non-native units (but only in one-way, this operation
does not need to be reversible).

From a mathematical point of view, these two behave very similar.
But from a practical point of view they are quite different.
Transforming the model state vector requires just the model state
vector of the native units and the Jacobian matrix that corresponds
to just that model state parameter.
Mapping the model state vector to non-native units may require
knowing more than just a single model state parameter.

Core mapping/transformation expression
--------------------------------------

If we define the native units of :math:`\vec{x}` as :math:`\vec{t}`
so that

.. math::

  \vec{x} = f\left(\vec{t}\right),

there must be a reversible functions so that

.. math::

  \vec{t} = f^{-1}\left(\vec{x}\right)

for any transformation or mapping to work.  It must also be possible
to take the partial derivative of :math:`\vec{t}` with regards
to :math:`\vec{x}`.

The Jacobian :math:`\mathbf{J}` is the derivative with respect to
the retrieved coordinates :math:`\vec{x}`.  The forward model first
computes derivatives in its native coordinates :math:`\vec{t}`.
If we introduce

.. math::

  \mathbf{J}' = \frac{\partial \vec{y}}{\partial \vec{t}},

it is clear we can write

.. math::

  \mathbf{J} = \mathbf{J}' \frac{\partial}{\partial \vec{x}} f^{-1}\left(\vec{x}\right).

This step right here is what we consider the transformation
of the Jacobian matrix.
To make use of this style of transformation, we must provide
matching :math:`f` and :math:`f^{-1}`, as well as a way to compute
the partial derivative of :math:`f^{-1}` with regards to
:math:`\vec{x}`.

See :doc:`user.oem` for assigning these operators to a Jacobian target.

Relative retrievals
^^^^^^^^^^^^^^^^^^^

This is a model state vector transformation.
By relative retrievals, we mean that the value itself is not
retrieved, but instead its ratio is retrieved.

In this scenario:

.. math::

  \vec{x} = \vec{t} \oslash \vec{t}_0,

.. math::

  \vec{t} = \vec{x} \odot \vec{t}_0,

.. math::

  \mathbf{J} = \mathbf{J}' \odot \vec{t}_0,

where :math:`\oslash` and :math:`\odot`
are element-wise division and multiplication,
respectively.  :math:`\vec{t}_0` is
simply the a priori value of :math:`\vec{t}`.

.. note::

  The first iteration of a retrieval setup is going to be :math:`\vec{x} = \vec{1}`.

Logarithmic retrievals
^^^^^^^^^^^^^^^^^^^^^^

This is a model state vector transformation.
By logarithmic retrievals, we mean that the value itself is not
retrieved, but instead its logarithm is retrieved.

In this scenario:

.. math::

  \vec{x} = \log\left(\vec{t}\right),

.. math::

  \vec{t} = \exp\left(\vec{x}\right),

.. math::

  \mathbf{J} = \mathbf{J}' \odot \exp\left(\vec{x}\right),

where the exponential and logarithmic operations are element-wise.

Logarithmic relative retrievals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a model state vector transformation.
By logarithmic relative retrievals, we mean that the value itself is not
retrieved, but instead the logarithm of its relative value is retrieved.

In this scenario:

.. math::

  \vec{x} = \log\left(\vec{t} \oslash \vec{t}_0\right),

.. math::

  \vec{t} = \exp\left(\vec{x}\right) \odot \vec{t}_0,

.. math::

  \mathbf{J} = \mathbf{J}' \odot \exp\left(\vec{x}\right) \odot \vec{t}_0,

where the operations are still element-wise on the product that is created.

.. note::

  The first iteration of a retrieval setup is going to have :math:`\vec{x} = \vec{0}`.

Relative humidity retrievals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a model state vector transformation.
By relative humidity retrievals, we mean that the value itself is not
retrieved, but instead the its conversion to relative humidity is retrieved.

In this scenario:

.. math::

  \vec{x} = \vec{t} \odot \vec{p} \oslash p_{\textrm{sat}}\left(\vec{T}\right),

.. math::

  \vec{t} = \vec{x} \odot p_{\textrm{sat}}\left(\vec{T}\right) \oslash \vec{p},

.. math::

  \mathbf{J} = \mathbf{J}' \odot p_{\textrm{sat}}\left(\vec{T}\right) \oslash \vec{p},

where
:math:`\vec{p}` is the pressure at the position of :math:`\vec{t}`,
:math:`\vec{T}` is the temperature at the position of :math:`\vec{t}`, and
:math:`p_{\textrm{sat}}` is a user-provided method to compute the element-wise
saturation pressure.

.. tip::

  There is a flag that can be provided to this transformation that
  turns negative relative humidities off.

.. note::

  Be aware that the implementation in ARTS is general,
  and that while you can choose to treat temperature as, e.g.,
  relative humidity... please don't.  It makes sense only for
  some species.

Absolute field retrievals
^^^^^^^^^^^^^^^^^^^^^^^^^

This is a model state vector mapping.
By absolute field retrievals, we mean that the value itself is not
retrieved, but instead the absolute value of the field is retrieved.

In this scenario:

.. math::

  \begin{array}{rcl}
    \vec{x} &=& \sqrt{\vec{t}_u \odot \vec{t}_u + \vec{t}_v \odot \vec{t}_v + \vec{t}_w \odot \vec{t}_w},\\
    \vec{\theta} &=& \arcsin\left(\vec{t}_w \oslash \vec{x}\right),\\
    \vec{\phi} &=& \arctan\left(\vec{t}_v \oslash \vec{t}_u\right),\\
    \vec{t}_u &=& \vec{x} \odot \cos\left(\vec{\theta}\right) \odot \cos\left(\vec{\phi}\right),\\
    \vec{t}_v &=& \vec{x} \odot \cos\left(\vec{\theta}\right) \odot \sin\left(\vec{\phi}\right),\\
    \vec{t}_w &=& \vec{x} \odot \sin\left(\vec{\theta}\right),\\
  \end{array}

.. math::

  \mathbf{J} = \mathbf{J}_u' \odot \vec{t}_u \oslash \vec{x} +
               \mathbf{J}_v' \odot \vec{t}_v \oslash \vec{x} +
               \mathbf{J}_w' \odot \vec{t}_w \oslash \vec{x},

where the subscripts :math:`u`, :math:`v`, and :math:`w`
indicate the three components of the vector north, east, and up,
respectively,
and :math:`\mathbf{J}_u'`, :math:`\mathbf{J}_v'`, and :math:`\mathbf{J}_w'`
are the Jacobian matrices of these three components
in the native units of the model state vector.

.. note::

  Neither :math:`\vec{\theta}` nor :math:`\vec{\phi}` are part of the
  model state vector, but are instead derived from the model state vector
  and are used to map the Cartesian coordinates to spherical coordinates.
  They are fixed during the retrieval process
  and are not updated between iterations.
