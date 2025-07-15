.. _Sec OEM:

Optimal estimation
###################

The core of optimal estimation in ARTS follows from :cite:t:`rodgers:00`.
He formulates the core expression

.. math::

  \vec{y} = F(\vec{x}),

where :\math:`\vec{y}` is a measurement vector
(i.e., :attr:`~pyarts.workspace.Workspace.measurement_vector` in the ARTS workspace),
:math:`\vec{x}` is the state of the model
(i.e., :attr:`~pyarts.workspace.Workspace.model_state_vector` in the ARTS workspace),
and :math:`F` is the model (i.e., ARTS).

Optimal estimation comes from linearization of the model itself.  Think instead of the
above as

.. math::

  \vec{y} = \mathbf{J} \vec{x},

where :math:`\mathbf{J}` is this linearization state, so that a naive way to get
the model state from a measurement may be written like

.. math::

  \vec{x} = \mathbf{J}^{-1} \vec{y}.

This :math:`\mathbf{J}` is the Jacobian matrix
(i.e., :attr:`~pyarts.workspace.Workspace.measurement_jacobian` in the ARTS workspace)

.. math::

  \mathbf{J} = \frac{\partial \vec{y}}{\partial \vec{x}}.


