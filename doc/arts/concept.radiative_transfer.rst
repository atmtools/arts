.. _Sec Radiative Transfer:

Radiative Transfer
##################

The step-by-step equation of ARTS radiative transfer calculation is described in this section.
The key variables to look out for are :attr:`~pyarts3.workspace.Workspace.spectral_radiance` and :attr:`~pyarts3.workspace.Workspace.spectral_radiance_jacobian`.

A single step
*************

A single step of the radiative transfer calculation in ARTS solvers assume to follow from

.. math::

  I_{1} = J_{0} + T_{0} \left(I_{0} - J_{0}\right),

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
  I_{0} &=& I_0(\vec{x}),

where :math:`x_0`, :math:`x_1` the levels above and below the layer, and :math:`\vec{x}` is the state of the model everywhere.
Note that both :math:`x_0` and :math:`x_1` are covered by :math:`\vec{x}`.  This will matter when we compute the Jacobian.

Multiple steps
**************

We are not always interested in the radiation at level 1, but quite often in the radiation
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
  \end{array} \right],

and so on.  Here, :math:`\vec{J}'` is the source vector, but we use a prime-tick because it is
0-indexed whereas :math:`\vec{I}` is 1-indexed.

In practice, it is often more convenient to just use
the step-by-step equations, but the matrix above
may help offer important insights about potential optimizations in the future.

Partial derivatives
*******************

The Jacobian after the last step is important for the retrieval of the state of the atmosphere.
This can be computed following

.. math::

  \frac{\partial I_{N}}{\partial\vec{x}} =
  \frac{\partial J_{N-1}}{\partial\vec{x}} +
  \frac{\partial T_{N-1}}{\partial\vec{x}} \left(I_{N-1} - J_{N-1} \right) +
  T_{N-1} \left( \frac{\partial I_{N-1}}{\partial\vec{x}} - \frac{\partial J_{N-1}}{\partial\vec{x}}\right)

As we know that the solution can be written as a sum of terms, we can expand the last term as

.. math::

  \begin{array}{llrl}
    \frac{\partial I_N}{\partial \vec{x}} &= T_{N-1}T_{N-2} \cdots T_1T_0\frac{\partial I_0}{\partial \vec{x}} + \frac{\partial J_{N-1}}{\partial \vec{x}}& \\
            &+& \left[
    \begin{array}{rrcrr}
      \frac{\partial T_{N-1}T_{N-2} \cdots T_1T_0}{\partial \vec{x}} &
      \frac{\partial T_{N-1}T_{N-2} \cdots T_2T_1}{\partial \vec{x}} &
      \cdots &
      \frac{\partial T_{N-1}T_{N-2}}{\partial \vec{x}} &
      \frac{\partial T_{N-1}}{\partial \vec{x}}
    \end{array}
     \right] &\left[
    \begin{array}{lcl}
      I_0 &-& J_0         \\
      J_0 &-& J_1         \\
      &\cdots&            \\
      J_{N-3} &-& J_{N-2} \\
      J_{N-2} &-& J_{N-1}
    \end{array} \right] \\
    &+&  \left[
    \begin{array}{rrcrr}
      T_{N-1}T_{N-2} \cdots T_1T_0 & T_{N-1}T_{N-2} \cdots T_2T_1 & \cdots & T_{N-1}T_{N-2} & T_{N-1}
    \end{array}
    \right] &\left[
    \begin{array}{lcl}
      &-& \frac{\partial J_0 }{\partial \vec{x}}        \\
      \frac{\partial J_0}{\partial \vec{x}} &-& \frac{\partial J_1}{\partial \vec{x}}         \\
      &\cdots&            \\
      \frac{\partial J_{N-3}}{\partial \vec{x}} &-& \frac{\partial J_{N-2}}{\partial \vec{x}} \\
      \frac{\partial J_{N-2}}{\partial \vec{x}} &-& \frac{\partial J_{N-1}}{\partial \vec{x}}
    \end{array} \right]
  \end{array}

Remember we already defined that
:math:`J_i` and :math:`T_i` depends on values at levels :math:`x_i` and :math:`x_{i+1}`.
These are part of :math:`\vec{x}` only via mapping, :math:`\vec{x}` covers both :math:`x_i` and :math:`x_{i+1}`.
Thus, the :math:`\frac{\partial I_{0}}{\partial\vec{x}}`-term has been lifted from the above expression,
as the other terms may be computed on a separate grid before being mapped back to :math:`\vec{x}`.
This mapping is not discussed here in details.  Introducing an alternative notation to make the expressions below
more compact,

.. math::
  
  \Pi_{n}^{m} = \left\{
  \begin{array}{ll}
    \prod_{i=n}^m T_i = T_n T_{n-1} \cdots T_{m+1} T_m & n \geq m \\
    1 & n < m
  \end{array}\right.,

where it is important to note that all :math:`T_i` are only functions of :math:`x_i` and :math:`x_{i+1}`.

For sake of keeping the expressions short, we add :math:`T_N=1` and assume :math:`N>>0` below.
We can extract the two dot products and expand them to see what the Jacobian looks like.
For the source partial derivatives:

.. math::

  \begin{array}{llll}
    \frac{\partial I_N^{(1)}}{\partial x_0} &=& &&
    \left(\Pi_N^1 - \Pi_N^0\right) &\frac{\partial J_0}{\partial x_0}
    \\
    \frac{\partial I_N^{(1)}}{\partial x_1} &=&
    \left(\Pi_N^1 - \Pi_N^0\right) &\frac{\partial J_0}{\partial x_1} &+&
    \left(\Pi_N^2 - \Pi_N^1\right) &\frac{\partial J_1}{\partial x_1}
    \\
    \frac{\partial I_N^{(1)}}{\partial x_2} &=&
    \left(\Pi_N^2 - \Pi_N^1\right) &\frac{\partial J_1}{\partial x_2} &+&
    \left(\Pi_N^3 - \Pi_N^2\right) &\frac{\partial J_2}{\partial x_2}
    \\
    \cdots
    \\
    \frac{\partial I_N^{(1)}}{\partial x_{N-2}} &=&
    \left(\Pi_N^{N-2} - \Pi_N^{N-3}\right) &\frac{\partial J_{N-3}}{\partial x_{N-2}} &+&
    \left(\Pi_N^{N-1} - \Pi_N^{N-2}\right) &\frac{\partial J_{N-2}}{\partial x_{N-2}}
    \\
    \frac{\partial I_N^{(1)}}{\partial x_{N-1}} &=&
    \left(\Pi_N^{N-1} - \Pi_N^{N-2}\right) &\frac{\partial J_{N-2}}{\partial x_{N-1}} &-&
    \Pi_N^{N-1}&\frac{\partial J_{N-1}}{\partial x_{N-1}}
  \end{array}

and for the transmittance partial derivatives:

.. math::

  \begin{array}{rrrrrrrrrrrrrrrrrr}
    \frac{\partial I_N^{(2)}}{\partial x_0} &=&
    \Pi_{N}^{1} &\Bigl[& \frac{\partial T_{0}}{\partial x_{0}} \left(I_0 - J_0\right) &\Bigr]
    \\
    \frac{\partial I_N^{(2)}}{\partial x_1} &=&
    \Pi_{N}^{2}
    &\Bigl[&
    \left(\frac{\partial T_1}{\partial x_1}T_{0} + T_{1}\frac{\partial T_{0}}{\partial x_{1}}\right) \left(I_0 - J_0\right) &+&
    \frac{\partial T_{1}}{\partial x_{1}} \left(J_0 - J_1\right)&\Bigr]
    \\
    \cdots
    \\
    \frac{\partial I_N^{(2)}}{\partial x_{N-1}} &=&
    T_{N}
    &\Bigl[&
    \left(\frac{\partial T_{N-1}}{\partial x_{N-1}}T_{N-2} + T_{N-1}\frac{\partial T_{N-2}}{\partial x_{N-1}}\right) \Pi_{N-3}^{0} \left(I_0 - J_0\right) &+&
    \left(\frac{\partial T_{N-1}}{\partial x_{N-1}}T_{N-2} + T_{N-1}\frac{\partial T_{N-2}}{\partial x_{N-1}}\right) \Pi_{N-3}^{1} \left(J_0 - J_1\right) &+&
    \cdots &+&
    \left(\frac{\partial T_{N-1}}{\partial x_{N-1}}T_{N-2} + T_{N-1}\frac{\partial T_{N-2}}{\partial x_{N-1}}\right) \Pi_{N-3}^{N-2} \left(J_{N-3} - J_{N-2}\right) &+&
    \frac{\partial T_{N-1}}{\partial x_{N-1}} \left(J_{N-2} - J_{N-1}\right)&\Bigr]
  \end{array}

The expression in the grid of :math:`\vec{x}` is then the following:

.. math::

  \frac{\partial I_{N}}{\partial \vec{x}} =
  \Pi_N^0\frac{\partial I_0}{\partial \vec{x}} +
  f
  \left(
  \frac{\partial I_N^{(1)}}{\partial \vec{x}_i} +
  \frac{\partial I_N^{(2)}}{\partial \vec{x}_i} +
  \frac{\partial J_{N-1}}{\partial \vec{x}_i}
  \right),

where the last term is 0 for all but :math:`i=N` and :math:`i=N-1`
and where the function :math:`f` is defined as
the map from :math:`\vec{x}_i\rightarrow\vec{x}`.
