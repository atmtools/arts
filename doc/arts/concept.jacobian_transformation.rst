.. _Sec Jacobian Transformation:

Jacobian transformation
#######################

There are two ways perform transformations on the Jacobian
matrix.
You can change the measurement unit (e.g., from spectral radiance
to brightness temperature) or you can change the model state unit
(e.g., from volume-mixing-ratio to relative humidity).

The first type of transformations are performed by an
:class:`pyarts.arts.SpectralRadianceTransformOperator`
in ARTS.  Generally, you should just use one of the
provided enumeration values of
:class:`pyarts.arts.SpectralRadianceUnitType` to
set up these types of transformations.
The enumeration class also describes how the transformations
are done.

The second type of transformation are a bit more complicated
because it's generally required to be able to perform
the transformation both ways.  This is the case because
optimal estimation methods are iterative, and they must be
able to update the state of the model between iterations.

This section will describe how to use existing transformations
and extend to use custom transformations.

Core transformation
===================

The core expression for a measurement and the Jacobian
matrix are from

.. math::

  \vec{y} = \mathbf{J} \vec{x},

  \mathbf{J} = \frac{\partial \vec{y}}{\partial \vec{x}},

where :math:`\vec{x}` is the model state vector,
:math:`\vec{y}` is the measurement vector, and
:math:`\mathbf{J}` is the Jacobian matrix.

If we define the native units of :math:`\vec{x}` as :math:`\vec{t}`
so that

.. math::

  \vec{x} = f\left(\vec{t}\right),

there must be a reversible functions so that

.. math::

  \vec{t} = f^{-1}\left(\vec{x}\right)

for any transformation of this sort to work.

If we put this in the form of the first equation above,

.. math::

  \vec{y} = \mathbf{J} \vec{x} = \mathbf{J} f\left(\vec{t}\right).

Here :math:`\mathbf{J}` is still the partial derivative with regards to
:math:`\vec{x}`.  However, all partial derivatives will have been
computed in terms of :math:`\vec{t}`, since this is the native unit.
If we introduce

.. math::

  \mathbf{K} = \frac{\partial \vec{y}}{\partial \vec{t}},

it is clear we can write

.. math::

  \mathbf{J} = \mathbf{K} \frac{\partial}{\partial \vec{x}} f^{-1}\left(\vec{x}\right).

This step right here is what we consider the transformation
of the Jacobian matrix.
To make use of this style of transformation, we must provide
matching :math:`f` and :math:`f^{-1}`, as well as the
partial derivative of :math:`f^{-1}`.

We provide several such solutions built-in to ARTS
as listed below

Relative retrievals
-------------------

By relative retrievals, we mean that the value itself is not
retrieved, but instead its ratio is retrieved.

In this scenario:

.. math::

  \vec{x} = \vec{t} \odot  \vec{t}_0,

  \vec{t} = \vec{x} \oslash \vec{t}_0,

  \mathbf{J} = \mathbf{K} \odot \vec{t}_0,

where :math:`\odot` and :math:`\oslash`
are elemtwise multiplicaton and division,
respectively.  :math:`\vec{t}_0` is
simply the a priori value of :math:`\vec{t}`.

Logarithmic retrievals
-----------------------

By logarithmic retrievals, we mean that the value itself is not
retrieved, but instead its logarithm is retrieved.

In this scenario:

.. math::

  \vec{x} = \log\left(\vec{t}\right),

  \vec{t} = \exp\left(\vec{x}\right),

  \mathbf{J} = \mathbf{K} \odot \vec{t},

where the exponential and logarithmic operations are elementwise.

Logarithmic relative retrievals
-------------------------------

By logarithmic relative retrievals, we mean that the value itself is not
retrieved, but instead the logarithm of its relative value is retrieved.

In this scenario:

.. math::

  \vec{x} = \log\left(\vec{t} \odot \vec{t}_0\right),

  \vec{t} = \exp\left(\vec{x} \oslash \vec{t}_0\right),

  \mathbf{J} = \mathbf{K} \odot \vec{t}_0 \odot \exp\left(\vec{x}\right),

where the operations chain.
