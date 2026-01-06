.. _Sec Radiative Transfer:

Radiative Transfer
##################

The step-by-step equation of ARTS radiative transfer calculation is described in this section.
The key variables to look out for are :attr:`~pyarts3.workspace.Workspace.spectral_rad` and :attr:`~pyarts3.workspace.Workspace.spectral_rad_jac`.

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

  J_{0} &=& J(\vec{x}_0, \vec{x}_1) \\
  T_{0} &=& T(\vec{x}_0, \vec{x}_1) \\
  I_{0} &=& I_0(\mathbf{x}),

where :math:`\vec{x}_0`, :math:`\vec{x}_1` is the state of the levels above and below the layer, and :math:`\mathbf{x}` is the state of the model everywhere.
Both :math:`\vec{x}_0` and :math:`\vec{x}_1` are covered by :math:`\mathbf{x}`, though they are not necessarily in :math:`\mathbf{x}`,
but can be derived from it, e.g., via interpolation.  This will matter when we compute the Jacobian - the Jacobian is on :math:`\mathbf{x}` but radiative transfer quantities
are on :math:`\vec{x}_0` and :math:`\vec{x}_1`, etc, so there must always be a map back from :math:`\vec{x}_0` and :math:`\vec{x}_1`, etc, to :math:`\mathbf{x}`.

.. note::

  The notation of ARTS often calls the content of :math:`\vec{x}_i` as targets of the Jacobian and
  a partially flattened version of :math:`\mathbf{x}` as the state vector.

  The unflattened dimension being the measurement dimension, e.g., frequency or channels.

  The concept of targets comes from the fact that the Jacobian is only computed for
  user-selected variables in :math:`\mathbf{x}`.  These selected variables are thus the
  users' "target".

Multiple steps
**************

We are not always interested in the radiation at level 1, but quite often in the radiation
at the end of the atmospheric path.  This, and the radiation at intermittent levels, can be computed
by simply taking multiple steps of the above equation.  Thus, going from level 0 to level N, we get

.. math::

  I_{1}   &=& J_{0} + T_{0} \left(I_{0} - J_{0}\right) \\
  I_{2}   &=& J_{1} + T_{1} \left(I_{1} - J_{1}\right) \\
  I_{3}   &=& J_{2} + T_{2} \left(I_{2} - J_{2}\right) \\
  I_{4}   &=& J_{3} + T_{3} \left(I_{3} - J_{3}\right) \\
  I_{5}   &=& J_{4} + T_{4} \left(I_{4} - J_{4}\right) \\
  &\cdots& \\
  I_{N} &=& J_{N-1} + T_{N-1} \left(I_{N-1} - J_{N-1}\right)

Both :math:`J_i` and :math:`T_i` depend on values at levels :math:`\vec{x}_i` and :math:`\vec{x}_{i+1}`,
even though these have been omitted for brevity.
There are here :math:`N-1` steps or layers, and :math:`N` levels.  Note that these steps are
non-commutative, since the step through a layer polarizes the radiation.

.. note::

  :math:`I_0` is often called the background radiation, as it is the radiation
  entering the atmosphere from space, from the ground, or elsewhere.
  The background, it other words.

If we expand the above expressions and collect them by transmittance, we get

.. math::

  I_1 &=& J_{0}   &+&                   T_{0} \left(I_{0} - J_{0}\right)  \\
  I_2 &=& J_{1}   &+&             T_{1} T_{0} \left(I_{0} - J_{0}\right)      &+&             T_{1} \left(J_{0} - J_{1}\right)  \\
  I_3 &=& J_{2}   &+&       T_{2} T_{1} T_{0} \left(I_{0} - J_{0}\right)      &+&       T_{2} T_{1} \left(J_{0} - J_{1}\right) &+&       T_{2} \left(J_{1} - J_{2}\right)  \\
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
0-indexed whereas the resulting :math:`\vec{I}` is 1-indexed.

In practice, it is often more convenient to just use
the step-by-step equations for forward calculations, but the matrix above
gives insight on optimizations for the Jacobian of the "measured" :math:`I_N`.

Partial derivatives
*******************

The Jacobian after the last step is important for the retrieval of the state of the atmosphere.
It is in this formalism defined as

.. math::

  \frac{\partial I_{N}}{\partial\mathbf{x}} =
  \frac{\partial J_{N-1}}{\partial\mathbf{x}} +
  \frac{\partial T_{N-1}}{\partial\mathbf{x}} \left(I_{N-1} - J_{N-1} \right) +
  T_{N-1} \left( \frac{\partial I_{N-1}}{\partial\mathbf{x}} - \frac{\partial J_{N-1}}{\partial\mathbf{x}}\right)

We now introduce an alternative notation to make the expressions below
more compact,

.. math::
  
  \Pi_{n}^{m} = 
    \prod_{i=n}^m T_i = T_n T_{n-1} \cdots T_{m+1} T_m,

where for sake of keeping the expressions less cluttered, we add :math:`T_N=1`.

We can summarize the full Jacobian after N steps as

.. math::

  \begin{array}{llrl}
    \frac{\partial I_N}{\partial \mathbf{x}} &=& 
    \Pi_{N}^{0}\frac{\partial I_0}{\partial \mathbf{x}} + \frac{\partial J_{N-1}}{\partial \mathbf{x}} & \\
    &+& \left[
    \begin{array}{rrcrr}
      \frac{\partial \Pi_{N}^{0}}{\partial \mathbf{x}} &
      \frac{\partial \Pi_{N}^{1}}{\partial \mathbf{x}} &
      \cdots &
      \frac{\partial \Pi_{N}^{N-2}}{\partial \mathbf{x}} &
      \frac{\partial \Pi_{N}^{N-1}}{\partial \mathbf{x}} 
    \end{array}
     \right] & \left[
    \begin{array}{lcl}
      I_0 &-& J_0         \\
      J_0 &-& J_1         \\
      &\cdots&            \\
      J_{N-3} &-& J_{N-2} \\
      J_{N-2} &-& J_{N-1}
    \end{array} \right] \\
    &+&  \left[
    \begin{array}{rrcrr}
      \Pi_{N}^{0} & \Pi_{N}^{1} & \cdots & \Pi_{N}^{N-2} & \Pi_{N}^{N-1}
    \end{array}
    \right] & \left[
    \begin{array}{lcl}
      &-& \frac{\partial J_0 }{\partial \mathbf{x}}        \\
      \frac{\partial J_0}{\partial \mathbf{x}} &-& \frac{\partial J_1}{\partial \mathbf{x}}         \\
      &\cdots&            \\
      \frac{\partial J_{N-3}}{\partial \mathbf{x}} &-& \frac{\partial J_{N-2}}{\partial \mathbf{x}} \\
      \frac{\partial J_{N-2}}{\partial \mathbf{x}} &-& \frac{\partial J_{N-1}}{\partial \mathbf{x}}
    \end{array} \right]
  \end{array}

The first term is the background radiation contribution, which is lifted from the
rest of the expression because we can then simplify the rest of the expression.
Since we know that we can map from :math:`\vec{x}_i` to :math:`\mathbf{x}`, we can
write the partial derivatives in terms of :math:`\vec{x}_i` instead of :math:`\mathbf{x}`.
The source term then becomes

.. math::

  \begin{array}{llllllll}
    \frac{\partial I_N^{(1)}}{\partial \vec{x}_0} &=& &&&&
    \Pi_N^1 & \left(1 - T_0\right)& \frac{\partial J_0}{\partial \vec{x}_0}
    \\
    \frac{\partial I_N^{(1)}}{\partial \vec{x}_1} &=&
    \Pi_N^1 & \left(1 - T_0\right)& \frac{\partial J_0}{\partial \vec{x}_1} &+&
    \Pi_N^2 & \left(1 - T_1\right)& \frac{\partial J_1}{\partial \vec{x}_1}
    \\
    \frac{\partial I_N^{(1)}}{\partial \vec{x}_2} &=&
    \Pi_N^2 & \left(1 - T_1\right)& \frac{\partial J_1}{\partial \vec{x}_2} &+&
    \Pi_N^3 & \left(1 - T_2\right)& \frac{\partial J_2}{\partial \vec{x}_2}
    \\
    &\cdots&
    \\
    \frac{\partial I_N^{(1)}}{\partial \vec{x}_{N-1}} &=&
    \Pi_N^{N-1} & \left(1 - T_{N-2}\right)& \frac{\partial J_{N-2}}{\partial \vec{x}_{N-1}} &+&
    \Pi_N^{N} & \left(1 - T_{N-1}\right)& \frac{\partial J_{N-1}}{\partial \vec{x}_{N-1}}
    \\
    \frac{\partial I_N^{(1)}}{\partial \vec{x}_{N}} &=&
    \Pi_N^{N} & \left(1 - T_{N-1}\right)& \frac{\partial J_{N-1}}{\partial \vec{x}_{N}} 
  \end{array}

and for the transmittance, the term becomes

.. math::

  \begin{array}{llllllllll}
    \frac{\partial I_N^{(2)}}{\partial \vec{x}_0} &=&
    &&&&
    \Pi_{N}^{1} & \frac{\partial T_{0}}{\partial \vec{x}_{0}}& \left(I_0 - J_0\right) 
    \\
    \frac{\partial I_N^{(2)}}{\partial \vec{x}_1} &=&
    \Pi_{N}^{1} & \frac{\partial T_0}{\partial\vec{x}_1} & \left(I_0 - J_0\right) &+&
    \Pi_{N}^{2} & \frac{\partial T_1}{\partial\vec{x}_1} & \left(I_1 - J_1\right)
    \\
    \frac{\partial I_N^{(2)}}{\partial \vec{x}_2} &=&
    \Pi_{N}^{2} & \frac{\partial T_1}{\partial\vec{x}_2} & \left(I_1 - J_1\right) &+&
    \Pi_{N}^{3} & \frac{\partial T_2}{\partial\vec{x}_2} & \left(I_2 - J_2\right)
    \\
    &\cdots&
    \\
    \frac{\partial I_N^{(2)}}{\partial \vec{x}_{N-1}} &=&
    \Pi_{N}^{N-1} & \frac{\partial T_{N-2}}{\partial\vec{x}_{N-1}} & \left(I_{N-2} - J_{N-2}\right) &+&
    \Pi_{N}^{N} & \frac{\partial T_{N-1}}{\partial\vec{x}_{N-1}} & \left(I_{N-1} - J_{N-1}\right)
    \\
    \frac{\partial I_N^{(2)}}{\partial \vec{x}_{N}} &=&
    \Pi_{N}^{N} & \frac{\partial T_{N-1}}{\partial\vec{x}_{N}} & \left(I_{N-1} - J_{N-1}\right) &&&
  \end{array}

The expression is then mapped back to :math:`\mathbf{x}` in the following:

.. math::

  \frac{\partial I_{N}}{\partial \mathbf{x}} =
  \Pi_N^0\frac{\partial I_0}{\partial \mathbf{x}} +
  \sum_i^N f
  \left(
  \frac{\partial I_N^{(1)}}{\partial \vec{x}_i} +
  \frac{\partial I_N^{(2)}}{\partial \vec{x}_i} +
  \frac{\partial J_{N-1}}{\partial \vec{x}_i}
  \right),

where the last term is 0 for all but :math:`i=N` and :math:`i=N-1`
and where the function :math:`f` is defined as
the map from :math:`\vec{x}_i\rightarrow\mathbf{x}`.
Often in ARTS, :math:`f` is just an inverse interpolation operator.

Potential improvements
**********************

Note that here, the indexing is not the same as above as
it is experimental notation.  This section exists mostly
as a reminder that we are using constant source and constant
propagation matrix above.  If these are allowed to change
within a layer, the equations below offers some approach to
achieving it.

.. note::

  These are untested and not how ARTS compute things,
  they are left as a future exercise so that we can
  implement them later.  Just the concept, for instance,
  of a Dawson function of a matrix is questionable.

Linear source function
^^^^^^^^^^^^^^^^^^^^^^

This is based on the assumption of a linearly changing source function
and constant propagation matrix through a layer.  The indexing is the same for the
transmittance and incoming background radiation as above - in layers.  However, the source
function indexing below is now shifted to represent levels - the same indexing
as we use for the spectral radiance term.

Linear in source means that we assume :math:`J = J_0 + (J_1 - J_0) \frac{r}{r_x}`,
where :math:`r` is the distance through the layer and :math:`r_x` is the total
thickness of the layer.  With constant propagation matrix :math:`K`, we define
:math:`T_0 = \exp(-K_0 r_x)` and get

.. math::
  I_1 = T_0 I_0 -
  \left[\left(\log T_0\right)^{-1}\left(1 - T_0\right) + T_0\right]J_0 +
  \left[1 + \left(\log T_0\right)^{-1}\left(1 - T_0\right)\right]J_1.

or

.. math::
  I_1 = J_1 + T_0        \left(I_0 - J_0\right) +
  \frac{1}{r_x} K_0^{-1} \left(1   - T_0\right) \left(J_0 - J_1\right).

This is a classical solution to the radiative transfer equation to improve
numerical stability when the propagation matrix :math:`K` is
large.  Going through the motion of extending this to multiple steps as for the
constant source case is then just a matter of adding 1 extra term
per layer in the dot products above.  If we define

.. math::
  \Lambda_i = \frac{1}{r_i} K_i^{-1} \left(1 - T_i\right),

we get

.. math::

  I_{1}   &=& J_{1} &+& T_{0} \left(I_{0} - J_{0}\right)       &+& \Lambda_0 \left(J_0 - J_1\right) \\
  I_{2}   &=& J_{2} &+& T_{1} \left(I_{1} - J_{1}\right)       &+& \Lambda_1 \left(J_1 - J_2\right) \\
  I_{3}   &=& J_{3} &+& T_{2} \left(I_{2} - J_{2}\right)       &+& \Lambda_2 \left(J_2 - J_3\right) \\
  I_{4}   &=& J_{4} &+& T_{3} \left(I_{3} - J_{3}\right)       &+& \Lambda_3 \left(J_3 - J_4\right) \\
  I_{5}   &=& J_{5} &+& T_{4} \left(I_{4} - J_{4}\right)       &+& \Lambda_4 \left(J_4 - J_5\right) \\
  &\cdots& \\
  I_{N}   &=& J_{N} &+& T_{N-1} \left(I_{N-1} - J_{N-1}\right) &+& \Lambda_{N-1} \left(J_{N-1} - J_{N}\right)

and again, we can expand these expressions and collect them by transmittance,
but this time with an additional term per layer

.. math::
  
  \begin{array}{rrrrrrrrrrrrr}
  I_1 &=& J_1 &+&                       T_0(I_0 - J_0) &+&         \Lambda_0(J_0 - J_1) \\
  I_2 &=& J_2 &+&                   T_1 T_0(I_0 - J_0) &+&     T_1 \Lambda_0(J_0 - J_1) &+&     \Lambda_1(J_1 - J_2) \\
  I_3 &=& J_3 &+&               T_2 T_1 T_0(I_0 - J_0) &+& T_2 T_1 \Lambda_0(J_0 - J_1) &+& T_2 \Lambda_1(J_1 - J_2) &+& \Lambda_2(J_2 - J_3) \\
  &\cdots& \\
  I_N &=& J_N &+& T_{N-1}\cdots T_2 T_1 T_0(I_0 - J_0) &+& \cdots                       &+& \cdots                   &+& \cdots &+& \Lambda_{N-1}(J_{N-1} - J_N) \\
  \end{array}

Which in shortened matrix form is:

.. math::
  
  \vec{I} = \vec{J} &+& \left[
  \begin{array}{lrrcrr}
    \Pi_0^0     &             \Lambda_0 & 0                     & \cdots & 0                             & 0 \\
    \Pi_1^0     & \Pi_1^1     \Lambda_0 &             \Lambda_1 & \cdots & 0                             & 0 \\
    \cdots      & \cdots                & \cdots                & \cdots & 0                             & 0 \\
    \Pi_{N-2}^0 & \Pi_{N-2}^1 \Lambda_0 & \Pi_{N-2}^2 \Lambda_1 & \cdots &                 \Lambda_{N-2} & 0 \\
    \Pi_{N-1}^0 & \Pi_{N-1}^1 \Lambda_0 & \Pi_{N-1}^2 \Lambda_1 & \cdots & \Pi_{N-1}^{N-1} \Lambda_{N-2} & \Lambda_{N-1}
  \end{array}
  \right] \left[
  \begin{array}{lcll}
    I_0 &-& J_0         \\
    J_0 &-& J_1         \\
    J_1 &-& J_2         \\
    &\cdots&            \\
    J_{N-2} &-& J_{N-1} \\
    J_{N-1} &-& J_{N}
  \end{array} \right]

The :math:`\vec{J}` term is now without a prime-tick as
it is 1-indexed like :math:`\vec{I}`,
though :math:`J_0` is still used in the RHS vector.

The partial derivative propagation is then derivable in the same way as above,

.. math::

  \frac{\partial I_{N}}{\partial\mathbf{x}} =
  \frac{\partial J_{N}}{\partial\mathbf{x}} +
  \frac{\partial T_{N-1}}{\partial\mathbf{x}} \left(I_{N-1} - J_{N-1} \right) +
  T_{N-1} \left( \frac{\partial I_{N-1}}{\partial\mathbf{x}} - \frac{\partial J_{N-1}}{\partial\mathbf{x}}\right) +
  \frac{\partial \Lambda_{N-1}}{\partial\mathbf{x}} \left(J_{N-1} - J_{N}\right) +
  \Lambda_{N-1} \left(\frac{\partial J_{N-1}}{\partial\mathbf{x}} - \frac{\partial J_{N}}{\partial\mathbf{x}}\right)

And again doing the expansion of the dot products and using :math:`T_N=1` for brevity, we get for the source partial derivatives at the end of the path:

.. math::
  \frac{\partial I_N}{\partial \mathbf{x}} &=& \Pi_{N}^0 \frac{\partial I_0}{\partial \mathbf{x}} + \frac{\partial J_{N}}{\partial \mathbf{x}}
  \\
  &+& \left[
  \begin{array}{rrrcrr}
    \frac{\partial \Pi_{N}^0}{\partial \mathbf{x}} &
    \frac{\partial \Pi_{N}^1 \Lambda_0}{\partial \mathbf{x}} &
    \frac{\partial}{\partial \mathbf{x}}\Pi_{N}^2 \Lambda_1 &
    \cdots &
    \frac{\partial \Pi_{N}^{N-1} \Lambda_{N-2}}{\partial \mathbf{x}} &
    \frac{\partial \Lambda_{N-1}}{\partial \mathbf{x}}
  \end{array} \right] & \left[
  \begin{array}{lcll}
    I_0 &-& J_0         \\
    J_0 &-& J_1         \\
    J_1 &-& J_2         \\
    &\cdots&            \\
    J_{N-2} &-& J_{N-1} \\
    J_{N-1} &-& J_{N}
  \end{array} \right]
  \\
  &+& \left[
  \begin{array}{rrrcrr}
    \Pi_{N}^0 & \Pi_{N}^1 \Lambda_0 & \Pi_{N}^2 \Lambda_1 & \cdots & \Pi_{N}^{N-1} \Lambda_{N-2} & \Pi_{N}^{N} \Lambda_{N-1}
  \end{array} \right] & \left[
  \begin{array}{lcll}
                                                 &-& \frac{\partial J_0}    {\partial \mathbf{x}} \\
    \frac{\partial J_0}    {\partial \mathbf{x}} &-& \frac{\partial J_1}    {\partial \mathbf{x}} \\
    \frac{\partial J_1}    {\partial \mathbf{x}} &-& \frac{\partial J_2}    {\partial \mathbf{x}} \\
    &\cdots&                                                                                      \\
    \frac{\partial J_{N-2}}{\partial \mathbf{x}} &-& \frac{\partial J_{N-1}}{\partial \mathbf{x}} \\
    \frac{\partial J_{N-1}}{\partial \mathbf{x}} &-& \frac{\partial J_{N}}{\partial \mathbf{x}}
  \end{array} \right]

Following the same procedure as above, we can rewrite this in terms of :math:`\vec{x}_i` instead of :math:`\mathbf{x}`.
Here, a key difference from before is that the source term now only depends on :math:`\vec{x}_i` and not on :math:`\vec{x}_{i+1}`.

The derivative contribution is then given by

.. math::
  \frac{\partial I_N} {\partial \vec{x}_0} &=& \Pi_N^1 &\Bigl[&
    \frac{\partial T_0} {\partial \vec{x}_0} \left(I_0 - J_0\right) &+&
    \frac{\partial \Lambda_0} {\partial \vec{x}_0} \left(J_0 - J_1\right) &+&
    \left(\Lambda_0 - T_0\right) \frac{\partial J_0} {\partial \vec{x}_0}
    &\Bigr]
  \\
  \frac{\partial I_N} {\partial \vec{x}_1} &=& \Pi_N^2 &\Bigl[&
    \frac{\partial T_1} {\partial \vec{x}_1} \left(I_1 - J_1\right) &+&
    \frac{\partial \Lambda_1} {\partial \vec{x}_1} \left(J_1 - J_2\right) &+&
    \left(\Lambda_1 - T_1\right) \frac{\partial J_1} {\partial \vec{x}_1}
    &\Bigr] &+&
    \Pi_N^1 &\Bigl[&
      \left(1 - \Lambda_0\right) \frac{\partial J_1} {\partial \vec{x}_1}
      &+& \frac{\partial T_0} {\partial \vec{x}_1} \left(I_0 - J_0\right)
      &+& \frac{\partial \Lambda_0} {\partial \vec{x}_1} \left(J_0 - J_1\right)
    &\Bigr]
  \\
  \frac{\partial I_N} {\partial \vec{x}_2} &=& \Pi_N^3 &\Bigl[&
    \frac{\partial T_2} {\partial \vec{x}_2} \left(I_2 - J_2\right) &+&
    \frac{\partial \Lambda_2} {\partial \vec{x}_2} \left(J_2 - J_3\right) &+&
    \left(\Lambda_2 - T_2\right) \frac{\partial J_2} {\partial \vec{x}_2}
    &\Bigr] &+&
    \Pi_N^2 &\Bigl[&
      \left(1 - \Lambda_1\right) \frac{\partial J_2} {\partial \vec{x}_2}
      &+& \frac{\partial T_1} {\partial \vec{x}_2} \left(I_1 - J_1\right)
      &+& \frac{\partial \Lambda_1} {\partial \vec{x}_2} \left(J_1 - J_2\right)
    &\Bigr]
  \\
  &\cdots&
  \\
  \frac{\partial I_N} {\partial \vec{x}_{N-1}} &=& \Pi_N^{N} &\Bigl[&
    \frac{\partial T_{N-1}} {\partial \vec{x}_{N-1}} \left(I_{N-1} - J_{N-1}\right) &+&
    \frac{\partial \Lambda_{N-1}} {\partial \vec{x}_{N-1}} \left(J_{N-1} - J_N\right) &+&
    \left(\Lambda_{N-1}-T_{N-1}\right) \frac{\partial J_{N-1}} {\partial \vec{x}_{N-1}}
    &\Bigr] &+&
    \Pi_N^{N-1} &\Bigl[&
      \left(1 - \Lambda_{N-2}\right) \frac{\partial J_{N-1}} {\partial \vec{x}_{N-1}}
      &+& \frac{\partial T_{N-2}} {\partial \vec{x}_{N-1}} \left(I_{N-2} - J_{N-2}\right)
      &+& \frac{\partial \Lambda_{N-2}} {\partial \vec{x}_{N-1}} \left(J_{N-2} - J_{N-1}\right)
    &\Bigr]
  \\
  \frac{\partial I_N} {\partial \vec{x}_{N}} &=& 
    &&&&&&&&&
    \Pi_N^{N} &\Bigl[&
      \left(1 - \Lambda_{N-1}\right) \frac{\partial J_{N}} {\partial \vec{x}_{N}}
      &+& \frac{\partial T_{N-1}} {\partial \vec{x}_{N}} \left(I_{N-1} - J_{N-1}\right)
      &+& \frac{\partial \Lambda_{N-1}} {\partial \vec{x}_{N}} \left(J_{N-1} - J_{N}\right)
    &\Bigr]

Linear propagation matrix and source function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This relaxes the assumption of constant propagation matrix within
the layer but is otherwise similar to linear in source.
The indexing is the same for the
transmittance and incoming background radiation as above - in layers.  However, the source
function and propagation matrix indexing below is now shifted to represent levels - the same indexing
as we use for the spectral radiance term.

Linear in source means that we assume :math:`J = J_0 + (J_1 - J_0) \frac{r}{r_x}`,
where :math:`r` is the distance through the layer and :math:`r_x` is the total
thickness of the layer.  With a linearly varying propagation matrix :math:`K = K_0 + (K_1 - K_0) \frac{r}{r_x}`, we define
:math:`T_0 = \exp\left(- \overline{K_{0, 1}} r_x\right)`, with the overline indicate
the arithmetic mean, and get

.. math::
  I_1 = J_1 + T_0 \left(I_0 - J_0\right) + \left(\Lambda_0^{(0)} + \Lambda_0^{(1)}\right) \left(J_0 - J_1\right).

This is a classical solution to the radiative transfer equation to improve
numerical stability when the propagation matrix :math:`K` is
large.  Going through the motion of extending this to multiple steps as for the
constant source case is then just a matter of adding 1 extra term
per layer in the dot products above.  If we define

.. math::
  \Lambda_i = \Lambda_i^{(0)} + \Lambda_i^{(1)},

where the first term is identical to the linear source function case,
and the second term is a second order correction for the linear
variation of the propagation matrix, we get


.. math::

  \Lambda_i^{(0)} &=& \frac{1}{r_i} \overline{K_{i, i+1}}^{-1} &\left(1 - T_i\right), \\
  \Lambda_i^{(1)} &=& \frac{1}{r_i^2} \overline{K_{i, i+1}}^{-1} &D\left(\left[K_{i+1} - K_{i}\right] r_i \middle/ 2\right),

where :math:`D` is the Dawson function of a matrix.

The rest of the steps that follows are the same as for the linear source function
case above, just replacing :math:`\Lambda_i` with the new definition.

.. caution::

  The Dawson function of a matrix is not implemented in ARTS yet.  Instead,
  the element-wise Dawson function is used as an approximation.  This means
  that it might be better to use the linear source function method.

  Indeed, we do not know how well this approximation performs in practice.
