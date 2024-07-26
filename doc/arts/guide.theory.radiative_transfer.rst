Radiative Transfer
==================

The step-by-step equation of ARTS radiative transfer calculation is described in this section.
The key variables to look out for are :attr:`~pyarts.Workspace.spectral_radiance` and :attr:`~pyarts.Workspace.spectral_radiance_jacobian`.

A single step
-------------

A single step of the radiative transfer calculation in ARTS solvers assume to follow from

.. math::

  I_{1} = J_{0} + T_{0} \left(I_{0} - J_{0}\right)

where
:math:`I_{1}` is the spectral radiance after the step,
:math:`I_{0}` is the spectral radiance before the step,
:math:`J_{0}` is the spectral radiance source at the step, and
:math:`T_{0}` is the transmittance through the current step.

To get more of an intuitive understanding of this equation,
a "step" could be considered as a single layer of the atmosphere.
Imagine that the radiation is coming into the layer at some level is :math:`I_{0}`,
and that the radiation is leaving at another level is :math:`I_{1}`.  The radiation that is leaving
is the sum of the radiation that is coming in, and the transmitted incoming radiation
differences to the local source of the radiation :math:`T_{0} \left(I_{0} - J_{0}\right)`,
and the radiation that is emitted by the layer itself :math:`J_{0}`

The way to compute :math:`J_{0}` and :math:`T_{0}` is not discussed here.
What is important is that they are computed based on the state of the model
in the surrounding levels of the model.  Noting that the radiation coming into the layer
may also have been reflected and can thus depend on all model states, we can write this as

.. math::

  J_{0} &=& J(x_0, x_1) \\
  T_{0} &=& T(x_0, x_1) \\
  I_{0} &=& I_0(\vec{x})

where :math:`x_0`, :math:`x_1` the levels above and below the layer, and :math:`\vec{x}` is the state of the model everywhere.
Note that both :math:`x_0` and :math:`x_1` are covered by :math:`\vec{x}`.  This will matter when we compute the Jacobian.

Multiple steps
--------------

Now we are not always interested in the radiation at level 0, but quite often in the radiation
at the end of the atmospheric path.  This, and the radiation at intermittent levels, can be computed
as

.. math::

  I_{2}   &=& J_{1} + T_{1} \left(I_{1} - J_{1}\right) \\
  I_{3}   &=& J_{2} + T_{2} \left(I_{2} - J_{2}\right) \\
  I_{4}   &=& J_{3} + T_{3} \left(I_{3} - J_{3}\right) \\
  I_{5}   &=& J_{4} + T_{4} \left(I_{4} - J_{4}\right) \\
  &\cdots& \\
  I_{N} &=& J_{N-1} + T_{N-1} \left(I_{N-1} - J_{N-1}\right)

Note that :math:`J_i` and :math:`T_i` depends on values at levels :math:`x_i` and :math:`x_{i+1}`,
even though these have been omitted for brevity.
There are here :math:`N-1` steps or layers, and :math:`N` levels.  Note that these steps are
non-commutative, since the step through a layer polarizes the radiation.

If we expand these expressions and collect them by transmittance, we get

.. math::

  I_1 &=& J_{0}   &+&                   T_{0} \left(I_{0} - J_{0}\right)  \\
  I_2 &=& J_{1}   &+&             T_{1} T_{0} \left(I_{0} - J_{0}\right)      &+&             T_{1} \left(J_{0} - J_{1}\right)  \\
  I_3 &=& J_{2}   &+&       T_{2} T_{1} T_{0} \left(I_{0} - J_{0}\right)      &+&       T_{2} T_{1} \left(J_{0} - J_{1}\right) &+&       T_{2} \left(J_{1} - J_{2}\right)  \\
  I_4 &=& J_{3}   &+& T_{3} T_{2} T_{1} T_{0} \left(I_{0} - J_{0}\right)      &+& T_{3} T_{2} T_{1} \left(J_{0} - J_{1}\right) &+& T_{3} T_{2} \left(J_{1} - J_{2}\right) &+& T_{3} \left(J_{2} - J_{3}\right) \\
  &\cdots& \\
  I_N &=& J_{N-1} &+& T_{N-1}T_{N-2} \cdots T_1T_0 \left(I_{0} - J_{0}\right) &+& \cdots                                       &+& \cdots                                 &+& T_{N-1} \left(J_{N-2} - J_{N-1}\right)

which we can rewrite as a simple sum of terms:

.. math::
  
  \vec{I} = \vec{J}' + \left[
  \begin{array}{rrcrr}
    T_0                          & 0                            & \cdots & 0              & 0 \\
    T_1T_0                       & T_1                          & \cdots & 0              & 0 \\
    \cdots                       & \cdots                       & \cdots & 0              & 0 \\
    T_{N-2}T_{N-3} \cdots T_1T_0 & T_{N-2}T_{N-3} \cdots T_2T_1 & \cdots & T_{N-2}        & 0 \\
    T_{N-1}T_{N-2} \cdots T_1T_0 & T_{N-1}T_{N-2} \cdots T_2T_1 & \cdots & T_{N-1}T_{N-2} & T_{N-1}
  \end{array}
  \right] \left[
  \begin{array}{lcl}
    I_0 &-& J_0         \\
    J_0 &-& J_1         \\
    &\cdots&            \\
    J_{N-3} &-& J_{N-2} \\
    J_{N-2} &-& J_{N-1}
  \end{array} \right]

and so on.  Here, :math:`\vec{J}'` is the source vector but we use a prime-tick because it is
0-indexed whereas :math:`\vec{I}` is 1-indexed.

In practice, it is often more convenient to just use
the step-by-step equations, but the matrix above
may help offer important insights about potential optimizations in the future.

Partial derivatives
-------------------

The Jacobian after the last step is important for the retrieval of the state of the atmosphere.
This can be computed following

.. math::

  \frac{\partial I_{N}}{\partial\vec{x}} =
  \frac{\partial J_{N-1}}{\partial\vec{x}} +
  \frac{\partial T_{N-1}}{\partial\vec{x}} \left(I_{N-1} - J_{N-1} \right) +
  T_{N-1} \left( \frac{\partial I_{N-1}}{\partial\vec{x}} - \frac{\partial J_{N-1}}{\partial\vec{x}}\right)

The only complicated term hese is the last one, which depends on the entire state of the atmosphere.
If we propagate the derivatives through the steps, we need 
:math:`\frac{\partial I_{0}}{\partial\vec{x}}` to start the process.
This can be seen if we expand the last term in the Jacobian, starting from the first step:

.. math::

  \begin{array}{lclclclclcr}
    \frac{\partial I_{1}}{\partial\vec{x}} &=&\frac{\partial J_{0}}{\partial\vec{x}} &+&\frac{\partial T_{0}}{\partial\vec{x}} \left(I_{0} - J_{0}\right) &-&T_{0} \frac{\partial J_{0}}{\partial\vec{x}} &&&+&T_{0} \frac{\partial I_{0}}{\partial\vec{x}} \\
    \frac{\partial I_{2}}{\partial\vec{x}} &=&\frac{\partial J_{1}}{\partial\vec{x}} &+&\frac{\partial T_{1}}{\partial\vec{x}} \left(I_{1} - J_{1}\right) &-&T_{1} \frac{\partial J_{1}}{\partial\vec{x}} &&&+&T_{1} \frac{\partial I_{1}}{\partial\vec{x}} \\
    &=&\frac{\partial J_{1}}{\partial\vec{x}} &+&\frac{\partial T_{1}}{\partial\vec{x}} \left(I_{1} - J_{1}\right) &-&T_{1} \frac{\partial J_{1}}{\partial\vec{x}} &+&T_{1} \left(\frac{\partial J_{0}}{\partial\vec{x}} +\frac{\partial T_{0}}{\partial\vec{x}} \left(I_{0} - J_{0}\right) -T_{0} \frac{\partial J_{0}}{\partial\vec{x}}\right) &+&T_{1}T_{0} \frac{\partial I_{0}}{\partial\vec{x}}\\
    \frac{\partial I_{3}}{\partial\vec{x}} &=&\frac{\partial J_{2}}{\partial\vec{x}} &+&\frac{\partial T_{2}}{\partial\vec{x}} \left(I_{2} - J_{2}\right) &-&T_{2} \frac{\partial J_{2}}{\partial\vec{x}} &&&+&T_{2} \frac{\partial I_{2}}{\partial\vec{x}} \\
    &=&\frac{\partial J_{2}}{\partial\vec{x}} &+&\frac{\partial T_{2}}{\partial\vec{x}} \left(I_{2} - J_{2}\right) &-&T_{2} \frac{\partial J_{2}}{\partial\vec{x}} &+&T_{2} \left(\frac{\partial J_{1}}{\partial\vec{x}} +\frac{\partial T_{1}}{\partial\vec{x}} \left(I_{1} - J_{1}\right) +T_{1} \left(\frac{\partial J_{0}}{\partial\vec{x}} +\frac{\partial T_{0}}{\partial\vec{x}} \left(I_{0} - J_{0}\right) -T_{0} \frac{\partial J_{0}}{\partial\vec{x}} -  \frac{\partial J_{1}}{\partial\vec{x}}\right)\right) &+&T_{2}T_{1}T_{0} \frac{\partial I_{0}}{\partial\vec{x}} \\
    &\cdots\\\frac{\partial I_{N}}{\partial\vec{x}} &=&\frac{\partial J_{N-1}}{\partial\vec{x}} &+&\frac{\partial T_{N-1}}{\partial\vec{x}} \left(I_{N-1} - J_{N-1}\right) &-&T_{N-1} \frac{\partial J_{N-1}}{\partial\vec{x}} &&&+&T_{N-1} \frac{\partial I_{N-1}}{\partial\vec{x}} \\
    &=&\frac{\partial J_{N-1}}{\partial\vec{x}} &+&\frac{\partial T_{N-1}}{\partial\vec{x}} \left(I_{N-1} - J_{N-1}\right) &-&T_{N-1} \frac{\partial J_{N-1}}{\partial\vec{x}} &+&T_{N-1} \left(\cdots\right) &+&T_{N-1}T_{N-2} \cdots T_{2}T_{1}T_{0} \frac{\partial I_{0}}{\partial\vec{x}}
  \end{array}

As can be seen, as long as we can determine
:math:`\frac{\partial I_{0}}{\partial\vec{x}}` at the start of the solver,
we can compute the Jacobian at the end of the solver.

Defining the partial *partial* derivatives as

.. math::

  \begin{array}{lclclcl}
    \frac{\partial I_{1}'}{\partial\vec{x}} &=&
    \frac{\partial J_{0}}{\partial\vec{x}}  &+&
    \frac{\partial T_{0}}{\partial\vec{x}} \left(I_{0} - J_{0}\right) &-&
    T_{0} \frac{\partial J_{0}}{\partial\vec{x}} \\
    \frac{\partial I_{2}'}{\partial\vec{x}} &=&
    \frac{\partial J_{1}}{\partial\vec{x}} &+&
    \frac{\partial T_{1}}{\partial\vec{x}} \left(I_{1} - J_{1}\right) &+&
    T_{1} \left(\frac{\partial I_{1}'}{\partial\vec{x}} - \frac{\partial J_{1}}{\partial\vec{x}}\right) \\
    \frac{\partial I_{3}'}{\partial\vec{x}} &=&
    \frac{\partial J_{2}}{\partial\vec{x}} &+&
    \frac{\partial T_{2}}{\partial\vec{x}} \left(I_{2} - J_{2}\right) &+&
    T_{2} \left(\frac{\partial I_{2}'}{\partial\vec{x}} - \frac{\partial J_{2}}{\partial\vec{x}}\right) \\
    &\cdots\\
    \frac{\partial I_{N}'}{\partial\vec{x}} &=&
    \frac{\partial J_{N-1}}{\partial\vec{x}} &+&
    \frac{\partial T_{N-1}}{\partial\vec{x}} \left(I_{N-1} - J_{N-1}\right) &+&
    T_{N-1} \left(\frac{\partial I_{N-1}'}{\partial\vec{x}} - \frac{\partial J_{N-1}}{\partial\vec{x}}\right)
  \end{array}


The only other terms are :math:`J_i` and :math:`T_i`.  These can be computed from the state of atmosphere or from the state of the model.
As mentioned before, we have :math:`J_i = J(x_i, x_{i+1})` and :math:`T_i = T(x_i, x_{i+1})`, where the :math:`x_i` are levels in the model.
Note that :math:`x_i` does not necessarily need to be a part of :math:`\vec{x}`, but that there must be a way to map :math:`x_i` back to :math:`\vec{x}`.

This significantly simplifies the expressions as we just need to compute 
:math:`\frac{\partial I_{N}'}{\partial x_i}`, and since
:math:`x_0` only affects :math:`J_0` and :math:`T_0`,
:math:`x_1` only affects :math:`J_0`, :math:`T_0`, :math:`J_1`, and :math:`T_1`,
and so on, we can compute the Jacobian by just computing the partial derivatives
of the source and transmittance functions in the local coordinate system.

In short:

.. math::

  \begin{array}{ll}
    \frac{\partial I_0}{\partial x_0} &=\left(1 - T_{0}\right) \frac{\partial}{\partial x_{0}} J_{0} - \frac{\partial}{\partial x_{0}} T_{0} J_{0}\\
    \frac{\partial I_1}{\partial x_0} &=T_{1} \frac{\partial}{\partial x_{0}} I_{0}\\
    \frac{\partial I_2}{\partial x_0} &=T_{2} \frac{\partial}{\partial x_{0}} I_{1}\\
    \frac{\partial I_3}{\partial x_0} &=T_{3} \frac{\partial}{\partial x_{0}} I_{2}\\
    \cdots\\
    \frac{\partial I_0}{\partial x_1} &=\left(1 - T_{0}\right) \frac{\partial}{\partial x_{1}} J_{0} - \frac{\partial}{\partial x_{1}} T_{0} J_{0}\\
    \frac{\partial I_1}{\partial x_1} &=T_{1} \frac{\partial}{\partial x_{1}} I_{0} - T_{1} \frac{\partial}{\partial x_{1}} J_{1} + \frac{\partial}{\partial x_{1}} J_{1} + \frac{\partial}{\partial x_{1}} T_{1} I_{0} - \frac{\partial}{\partial x_{1}} T_{1} J_{1}\\
    \frac{\partial I_2}{\partial x_1} &=T_{2} \frac{\partial}{\partial x_{1}} I_{1}\\
    \frac{\partial I_3}{\partial x_1} &=T_{3} \frac{\partial}{\partial x_{1}} I_{2}\\
    \cdots\\
    \frac{\partial I_0}{\partial x_2} &=0\\
    \frac{\partial I_1}{\partial x_2} &=- T_{1} \frac{\partial}{\partial x_{2}} J_{1} + \frac{\partial}{\partial x_{2}} J_{1} + \frac{\partial}{\partial x_{2}} T_{1} J_{0} - \frac{\partial}{\partial x_{2}} T_{1} J_{1} - \frac{\partial}{\partial x_{2}} T_{1} T_{0} J_{0}\\
    \frac{\partial I_2}{\partial x_2} &=T_{2} \frac{\partial}{\partial x_{2}} I_{1} - T_{2} \frac{\partial}{\partial x_{2}} J_{2} + \frac{\partial}{\partial x_{2}} J_{2} + \frac{\partial}{\partial x_{2}} T_{2} I_{1} - \frac{\partial}{\partial x_{2}} T_{2} J_{2}\\
    \frac{\partial I_3}{\partial x_2} &=T_{3} \frac{\partial}{\partial x_{2}} I_{2}\\
    \cdots\\
    \frac{\partial I_0}{\partial x_3} &=0\\
    \frac{\partial I_1}{\partial x_3} &=0\\
    \frac{\partial I_2}{\partial x_3} &=- T_{2} \frac{\partial}{\partial x_{3}} J_{2} + \frac{\partial}{\partial x_{3}} J_{2} + \frac{\partial}{\partial x_{3}} T_{2} J_{1} - \frac{\partial}{\partial x_{3}} T_{2} J_{2} + \frac{\partial}{\partial x_{3}} T_{2} T_{1} J_{0} - \frac{\partial}{\partial x_{3}} T_{2} T_{1} J_{1} - \frac{\partial}{\partial x_{3}} T_{2} T_{1} T_{0} J_{0}\\
    \frac{\partial I_3}{\partial x_3} &=T_{3} \frac{\partial}{\partial x_{3}} I_{2} - T_{3} \frac{\partial}{\partial x_{3}} J_{3} + \frac{\partial}{\partial x_{3}} J_{3} + \frac{\partial}{\partial x_{3}} T_{3} I_{2} - \frac{\partial}{\partial x_{3}} T_{3} J_{3}\\
    \cdots\\
    \frac{\partial I_0}{\partial x_4} &=0\\
    \frac{\partial I_1}{\partial x_4} &=0\\
    \frac{\partial I_2}{\partial x_4} &=0\\
    \frac{\partial I_3}{\partial x_4} &=- T_{3} \frac{\partial}{\partial x_{4}} J_{3} + \frac{\partial}{\partial x_{4}} J_{3} + \frac{\partial}{\partial x_{4}} T_{3} J_{2} - \frac{\partial}{\partial x_{4}} T_{3} J_{3} + \frac{\partial}{\partial x_{4}} T_{3} T_{2} J_{1} - \frac{\partial}{\partial x_{4}} T_{3} T_{2} J_{2} + \frac{\partial}{\partial x_{4}} T_{3} T_{2} T_{1} J_{0} - \frac{\partial}{\partial x_{4}} T_{3} T_{2} T_{1} J_{1} - \frac{\partial}{\partial x_{4}} T_{3} T_{2} T_{1} T_{0} J_{0}\\
  \end{array}
