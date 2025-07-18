.. _Sec OEM:

Optimal estimation
###################

The core expression of optimal estimation in ARTS
follows from :cite:t:`rodgers:00`.
He formulates the core expression of a measurement as

.. math::

  \vec{y} = F\left(\vec{x}\right) + \epsilon,

where :math:`\vec{y}` is a measurement vector
(i.e., :attr:`~pyarts.workspace.Workspace.measurement_vector`
in the ARTS workspace),
:math:`\vec{x}` is the state of the model
(i.e., :attr:`~pyarts.workspace.Workspace.model_state_vector`
in the ARTS workspace),
:math:`F` is the model (i.e., ARTS itself), and
:\math:`\epsilon` is some measurement error that cannot
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
:math:`\mathbf{S}_a` is the covariance of the model state apriori.

Linearization
=============

The optimal estimation methods in ARTS only work if it is possible to
linearize :math:`F\left(\vec{x}\right)` around :math:`\vec{x}`.
This is equivalent to stating that

.. math::

  \vec{y} = \vec{y}_f + \mathbf{J} \left(\vec{x} - \vec{x}_a\right)

where we can take the partial derivative to find that

.. math::

  \mathbf{J} = \frac{\partial \vec{y}}{\partial \vec{x}}

is the Jacobian matrix
(i.e., :attr:`~pyarts.workspace.Workspace.measurement_jacobian`
in the ARTS workspace).

One approach to minimize :math:`\vec{x}` is to
use a Gauss-Newton approach to update the state
of the atmosphere.  This might look like

.. math::
  \vec{x}_{i+1} = \vec{x}_a + \mathbf{S}_a
  \mathbf{J}^\top\left(\mathbf{J}\mathbf{S}_a \mathbf{J} ^\top+\mathbf{S}_\epsilon\right)^{-1}
  \left[\vec{y}-\vec{y}_f+
  \mathbf{J}\left(\vec{x}_i-\vec{x}_a\right)\right],

.. math::

  \vec{y}_f = F\left(\vec{x}_i\right)

for :math:`i` starting at :math:`a` or :math:`0` and increasing the
count until the cost functions above are not decreasing anymore, or
too slowly to warrant more iterations.

Transforming the Jacobian matrix
================================

It is often desired to transform the Jacobian matrix away from the
native units that are available in the model.  Instead of retrieving
volume-mixing-ratio, it is often more stable to retrieve the
relative change in the volume-mixing-ratio, or perhaps even
the logarithm of the volume-mixing-ratio, to avoid negative
values creeping in to what is otherwise an intensive and numerical
exercise.

There are two ways perform transformations on the Jacobian
matrix in ARTS.
You can change the measurement unit (e.g., from spectral radiance
to brightness temperature) or you can change the model state unit
(e.g., from volume-mixing-ratio to relative humidity).

The first type of transformations are performed by an
:class:`~pyarts.arts.SpectralRadianceTransformOperator`
in ARTS.  Generally, you should just use one of the
provided enumeration values of
:class:`~pyarts.arts.SpectralRadianceUnitType` to
set up these types of transformations.
The enumeration class also describes how the transformations
are done.

The second type of transformation are a bit more complicated
because it's generally required to be able to perform
the transformation both ways.  This is the case because
optimal estimation methods are iterative, and they must be
able to update the state of the model between iterations.

.. note::
  A third type of transformation is imaginable as a target
  selection transformation.
  That is, for example, the wind-field using all three components
  as targets but transforming the Jacobian matrix to work on
  only the absolute wind-speed.  This type of transformation is
  not available in ARTS at this time but would be a welcome
  addition.  This comment is thus ignored by the rest of this
  text.

Core transformation expression
------------------------------

If we define the native units of :math:`\vec{x}` as :math:`\vec{t}`
so that

.. math::

  \vec{x} = f\left(\vec{t}\right),

there must be a reversible functions so that

.. math::

  \vec{t} = f^{-1}\left(\vec{x}\right)

for any transformation of this sort to work.  It must also be possible
to take the partial derivative of :math:`\vec{t}` with regards
to :math:`\vec{x}`.

If we put this in the form of the linearized forward simulation,

.. math::

  \vec{y}_f = \mathbf{J} \vec{x} = \mathbf{J} f\left(\vec{t}\right).

Here :math:`\mathbf{J}` is still the partial derivative with regards to
:math:`\vec{x}`.  However, all partial derivatives will have been
computed in terms of :math:`\vec{t}`, since this is the native unit.
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

We provide several such solutions built-in to ARTS
as listed below but it is possible to specify these
directly from python by simply providing the three
operators above.

Relative retrievals
^^^^^^^^^^^^^^^^^^^

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
