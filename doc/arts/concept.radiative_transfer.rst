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

And again doing the expansion of the dot products and
using :math:`T_N=1` for brevity, we get for the source
partial derivatives at the end of the path:

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

Following the same procedure as above, we can rewrite this in
terms of :math:`\vec{x}_i` instead of :math:`\mathbf{x}`.
Here, a key difference from before is that the source term
now only depends on :math:`\vec{x}_i` and not on :math:`\vec{x}_{i+1}`.

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

This partly relaxes the assumption of constant propagation matrix within
the layer but is otherwise similar to linear in source.
The indexing is the same for the
transmittance and incoming background radiation as above - in layers.
However, the source
function and propagation matrix indexing below is now shifted
to represent levels - the same indexing
as we use for the spectral radiance term.

Linear in source means that we assume
:math:`J = J_0 + (J_1 - J_0) \frac{r}{r_x}`,
where :math:`r` is the distance through the layer and :math:`r_x` is the total
thickness of the layer.
With a linearly changing propagation matrix
:math:`K = K_0 + (K_1 - K_0) \frac{r}{r_x}`, we define
:math:`T_0 = \exp\left(- \overline{K_{0, 1}} r_x\right)`, and again get

.. math::
  I_1 = J_1 + T_0 \left(I_0 - J_0\right) + \Lambda_0 \left(J_0 - J_1\right),

where

.. math::
  \Lambda_i = \frac{1}{r_i} T_i \int_0^{r_i} \exp\left(K_i s + \frac{K_{i+1} - K_i}{2 r_i} s^2\right) ds.

The solution to this may be approximated with the
Matrix Dawson function :math:`\mathcal{D}`:

.. math::
  \mathbf{\Lambda}_i =
  \left\{
  \begin{array}{ll}
  \frac{1}{r_i} \mathbf{\alpha}^{-1} \left( \mathcal{D}(\mathbf{u}_1) - T_i \mathcal{D}(\mathbf{u}_0) \right) & \text{if } K_{i+1} \gg K_i \\
  \\
  \frac{1}{r_i} K_i^{-1} (1 - T_i) & \text{otherwise}
  \end{array},
  \right.

where :math:`\mathbf{\alpha} = \sqrt{(K_{i+1} - K_i)/(2 r_i)}`
(principal square root),
:math:`\mathbf{u}_0 = \frac{1}{2} \mathbf{\alpha}^{-1} K_i`,
and :math:`\mathbf{u}_1 = \frac{1}{2} \mathbf{\alpha}^{-1} K_{i+1}`.
The expression should be valid for imaginary input, as well as real,
however we find this very unstable numerically.  Thus, we only use this
approach when the propagation matrix increases significantly
along the path, i.e., :math:`K_{i+1} \gg K_i`.

The rest of the steps that follow are the same as for the
linear source function
case above, just replacing :math:`\Lambda_i` with this new definition.

.. admonition:: Caveats
  :class: warning

  The approach here uses a matrix Dawson function, which is not
  what is currently implemented in ARTS.  Instead, an element-wise
  Dawson function is used, which is only correct if the propagation
  matrix is diagonal in the same basis at both ends of the layer, or
  potentially if the non-commutative part is very small.

  The switch from using linear propagation matrix to constant propagation
  matrix when the propagation matrix decreases along the path is also
  a crude approximation that should be improved in the future.  Currently,
  the argument to the square root must be greater that :math:`10^{-6}` for
  the algorithms to not revert to the constant propagation matrix case.

  Also be aware that the square root of a matrix is not unique, leading to
  further potential numerical issues.

Deriving the expressions
************************

This is not an important section for users of ARTS, but is included
for completeness to show how the expressions above are derived.
It may also be useful if you want to implement
other variations of the radiative transfer equation solution.

The matrix case is the same, provided you can treat :math:`K` and :math:`J` as
*mostly* commuting.  This only affect the last case.  The path
variables are in :math:`s \in [0,r]` for the distance through
a layer of understood
physics, with :math:`I(0)=I_0` and :math:`I(r)=I_1`.

The differential equation in all cases is

.. math::

  \frac{dI}{ds} &= -K(s)\bigl(I(s) - J(s)\bigr)
  \quad\Longleftrightarrow

  I'(s) + K(s) I(s) &= K(s) J(s).

The integrating factor is

.. math::

  \mu(s) = \exp\Bigl(\int_0^s K(u)du\Bigr),

and then

.. math::

  \bigl(\mu I\bigr)' &= \mu K J
  \quad\Rightarrow\quad

  I(r) &= T(r)I_0 + T(r)\int_0^r \mu(s)K(s)J(s)ds

with

.. math::

  T(r) &= \mu(r)^{-1} \\
       &= \exp\Bigl(-\int_0^r K(u)du\Bigr).

All three cases are just different assumptions
for :math:`K(s)` and :math:`J(s)`.

1. Constant :math:`K`, constant :math:`J`.
   
   Take

   .. math::

      K(s)&=K_0 \\ J(s)&=J_0 \\ T_0 &= e^{-K_0 r}.

   Then

   .. math::

      \mu(s) &= e^{K_0 s} \\
      T(r)   &=e^{-K_0 r} \\
             &=T_0.

   The solution:

   .. math::

      I(r)
        &= T_0 I_0 + T_0 \int_0^r e^{K_0 s}K_0J_0ds \\
        &= T_0 I_0 + T_0 \left(\int_0^r e^{K_0 s}ds\right) K_0 J_0  \\
        &= T_0 I_0 + T_0 \left(e^{K_0 r} - 1\right) K_0^{-1} K_0 J_0  \\
        &= T_0 I_0 + T_0 \left(e^{K_0 r} - 1\right) J_0 \\
        &= T_0 I_0 + \left(1 - T_0\right) J_0.

   Writing :math:`I_1 = I(r)`, we get exactly:

   .. math::

      I_1 = J_0 + T_0 (I_0 - J_0),

   which is the form used above.  Note that no parts of this derivation
   depends on whether or not the expression is scalar or matrix.

2. Constant :math:`K`, linearly varying :math:`J`. Assume

   .. math::

      K(s)&=K_0 \\ J(s)&=J_0 + \frac{s}{r}(J_1 - J_0).

   Again :math:`T_0 = e^{-K_0 r}`, and :math:`\mu(s)=e^{K_0 s}`.

   Use the general solution:

   .. math::

      I(r) = T_0 I_0 + T_0 \int_0^r e^{K_0 s} K_0 J(s)ds.

   Insert :math:`J(s)` and split:

   .. math::

      I(r)
      &= T_0 I_0 + T_0 \int_0^r e^{K_0 s} K_0
      \Bigl[J_0 + \frac{s}{r}(J_1-J_0)\Bigr] ds \\
      &= T_0 I_0 +
      T_0 \left(\int_0^r e^{K_0 s} ds\right) K_0 J_0 +
      T_0 \left(\int_0^r s e^{K_0 s} ds\right) K_0 \frac{J_1-J_0}{r}.

   Compute the two integrals:

   .. math::

      \int_0^r e^{K_0 s} ds &= \left(e^{K_0 r}-1\right) K_0^{-1} \\
      \int_0^r s e^{K_0 s} ds
      &= \Bigl[se^{K_0 s}K_0^{-1}\Bigr]_0^r - \int_0^r e^{K_0 s} K_0^{-1} ds\\
      &= r e^{K_0 r} K_0^{-1} - \left(e^{K_0 r}-1\right)K_0^{-2}.

   Insert:

   .. math::

      I(r)
      = T_0 I_0
      +
      T_0 \left(e^{K_0 r}-1\right)K_0^{-1} K_0 J_0  +
      T_0 \Bigl[r e^{K_0 r}K_0^{-1} - \left(e^{K_0 r}-1\right){K_0^{-2}}\Bigr]
      K_0 \frac{J_1-J_0}{r}.

   Simplify each term, using :math:`T_0=e^{-K_0 r}` and cancelling
   where possible.

   First integral term:

   .. math::

      T_0 \left(e^{K_0 r}-1\right)K_0^{-1} K_0 J_0
      &= T_0 \left(e^{K_0 r}-1\right) J_0 \\
      &= \left( 1-T_0 \right) J_0.

   Second integral term:

   .. math::

      T_0
      \Bigl[r e^{K_0 r}K_0^{-1} - \left(e^{K_0 r}-1\right)K_0^{-2}\Bigr]
      K_0 \frac{J_1-J_0}{r}
      &=
      \left[
      T_0 r e^{K_0 r} +
      T_0 \left(e^{K_0 r}-1\right)K_0^{-1}
      \right] \frac{J_1-J_0}{r} \\
      &=
      \left[
      r -
      \left(1-T_0\right)K_0^{-1}
      \right] \frac{J_1-J_0}{r}.

   So

   .. math::

      I(r) = T_0 I_0 + \left( 1-T_0 \right) J_0
      +\left[r - \left( 1-T_0 \right)K_0^{-1}\right] \frac{J_1-J_0}{r}.

   Rearrange as

   .. math::

      I_1 = J_1 + T_0 (I_0 - J_0)
      + \underbrace{\frac{1}{r} K_0^{-1}(1-T_0)}_{\displaystyle \Lambda_0}(J_0 - J_1),

   which gives the form used above


   .. math::

      I_1 &= J_1 + T_0 (I_0 - J_0) + \Lambda_0 (J_0 - J_1),\\
      \Lambda_0 &= \frac{1}{r} K_0^{-1} \left( 1-T_0 \right) \\ &=-\log T_0 \left( 1-T_0 \right).

   Again, no assumptions were made about scalar vs. matrix.
   The only *scary* step of matrix notation non-commutativity is that you have to
   be aware that you can write, e.g., :math:`A^3=A A^2=A^2 A`,
   which means with :math:`\exp(A) = \sum_{n=0}^\infty \frac{A^n}{n!}`, you can do
   :math:`\exp(A) A = \left(\sum_{n=0}^\infty \frac{A^{n}}{n!}\right) A = \sum_{n=0}^\infty \frac{A^{n+1}}{n!} = A \sum_{n=0}^\infty \frac{A^{n}}{n!} = A \exp(A)`.

   .. admonition:: Implementation note
      :class: tip

      It is very important to implement the expression for :math:`\Lambda_0`
      in a numerically stable way.  The expression above is not
      stable for small :math:`K_0 r` (i.e., large :math:`T_0`) *as written*.
      The key instability stems from the subtraction of two nearly equal
      terms in :math:`1 - T_0`.  The IEEE floating point standard
      provides a function :code:`expm1(x)`, which computes :math:`e^x - 1`
      in a numerically stable way for small :math:`x`.

      Likewise, the matrix expansion of :math:`1 - T_0` might be unstable.

      So we use a special solution implementing our own version of the reduced
      Cayley-Hamilton theorem to compute :math:`\Lambda_0` in a numerically
      stable way
      for matrices that conform to the propagation matrix notation in ARTS.
      This makes use of the inversion of :math:`K` to remove components from
      the expansion of the matrix exponential that would otherwise
      cause numerical instability.

3. Linear :math:`K(s)` and linear :math:`J(s)`.

   Now set both :math:`K` and :math:`J` to
   be linear in :math:`s` across the layer:

   .. math::

      K(s) &= K_0 + \frac{s}{r}(K_1-K_0),\\
      J(s) &= J_0 + \frac{s}{r}(J_1-J_0).

   Let

   .. math::

      \alpha &= \frac{K_1-K_0}{r},\\
      \beta &= \frac{J_1-J_0}{r},\\
      K(s) &= K_0 + \alpha s,\\
      J(s) &= J_0 + \beta s.\\

   3.1. Integrating factor and transmittance.

     The integrating factor and transmittance become

     .. math::

        \mu(s)
        &= \exp\Bigl(\int_0^s K(u)du\Bigr)\\
        &= \exp\Bigl(\int_0^s (K_0 + \alpha u)du\Bigr) \\
        &= \exp\bigl(K_0 s + \tfrac12 \alpha s^2\bigr)\\
        T(r) &= \mu(r)^{-1}\\
        &= \exp\bigl(-K_0 r - \tfrac12 \alpha r^2\bigr)\\
        &= T_0.

     That's the same :math:`T_0` as always, assuming that
     we only evaluate the expression at :math:`0` and :math:`r`,
     and that the layer value of :math:`K` above is the average of
     the endpoints (:math:`K_0` and :math:`K_1`).

   3.2. General solution. The general solution is still

     .. math::

        I(r) = T_0 I_0 + T_0 \int_0^r \mu(s) K(s) J(s)\,ds.

     Insert the forms:

     .. math::

        \mu(s) &= e^{K_0 s + \frac12 \alpha s^2},\\
        K(s)   &= K_0 + \alpha s,\\
        J(s)   &= J_0 + \beta s\\
               &= J_1 + (J_0 - J_1)\Bigl(1 - \frac{s}{r}\Bigr),

     where in the last step we used that :math:`J(r)=J_1`.

     Then

     .. math::

        I(r)
        = T_0 I_0
        + T_0 \int_0^r \mu(s) K(s)
          \Bigl[J_1 + (J_0-J_1)\Bigl(1 - \frac{s}{r}\Bigr)\Bigr]\,ds.

     Split the integral into the part proportional to :math:`J_1` and
     the part proportional to :math:`(J_0-J_1)`:

     .. math::

        I(r)
        = T_0 I_0
        + T_0 \left(\int_0^r \mu(s) K(s)\,ds\right) J_1
        + T_0 \left( \int_0^r \mu(s) K(s)
                     \Bigl(1 - \frac{s}{r}\Bigr)\,ds \right) (J_0-J_1).

     The first integral is the same as for a *constant* source
     :math:`J(s)\equiv J_1`.  For any :math:`K(s)`, that problem has the
     known solution

     .. math::

        I(r) = J_1 + T_0 (I_0 - J_1),

     so comparing with the integral form gives

     .. math::

        T_0 \left(\int_0^r \mu(s) K(s)\,ds\right) J_1 &= (1 - T_0) J_1
        \quad\Rightarrow\quad\\
        \int_0^r \mu(s) K(s)\,ds &= T_0^{-1} - 1.

     Insert this back:

     .. math::

        I(r)
        = T_0 I_0 + (1-T_0) J_1
          + T_0 \left( \int_0^r \mu(s) K(s)
                       \Bigl(1 - \frac{s}{r}\Bigr)\,ds \right) (J_0-J_1).

     Rearranging the first two terms:

     .. math::

        T_0 I_0 + (1-T_0) J_1
        = J_1 + T_0 (I_0 - J_1)
        = J_1 + T_0 (I_0 - J_0) + T_0 (J_0 - J_1),

     so

     .. math::

        I(r)
        = J_1 + T_0 (I_0 - J_0)
        + \Biggl[
              T_0
              + T_0 \int_0^r \mu(s) K(s)\Bigl(1 - \frac{s}{r}\Bigr)\,ds
            \Biggr] (J_0-J_1).

     This shows that the solution has the affine form

     .. math::

        I_1 = J_1 + T_0 (I_0 - J_0) + \Lambda_0 (J_0 - J_1),

     with

     .. math::

        \Lambda_0
        = T_0
          \Biggl[
            1 + \int_0^r \mu(s) K(s)\Bigl(1 - \frac{s}{r}\Bigr)\,ds
          \Biggr].

     To obtain a more compact expression for :math:`\Lambda_0`,
     use integration by parts on the remaining integral.  Note that

     .. math::

        \bigl[\mu(s)\bigl(1 - \tfrac{s}{r}\bigr)\bigr]'
        = \mu'(s)\Bigl(1 - \frac{s}{r}\Bigr) - \frac{1}{r}\mu(s)
        = \mu(s)K(s)\Bigl(1 - \frac{s}{r}\Bigr) - \frac{1}{r}\mu(s),

     hence

     .. math::

        \mu(s)K(s)\Bigl(1 - \frac{s}{r}\Bigr)
        = \bigl[\mu(s)\bigl(1 - \tfrac{s}{r}\bigr)\bigr]' + \frac{1}{r}\mu(s).

     Integrating this from :math:`0` to :math:`r` gives

     .. math::

        \int_0^r \mu(s)K(s)\Bigl(1 - \frac{s}{r}\Bigr) ds
        &= \Bigl[\mu(s)\Bigl(1 - \frac{s}{r}\Bigr)\Bigr]_0^r
           + \frac{1}{r}\int_0^r \mu(s)\,ds \\
        &= -1 + \frac{1}{r}\int_0^r \mu(s)\,ds,

     because :math:`\mu(0)=1` and :math:`1-r/r=0`.  Inserting this into
     the expression for :math:`\Lambda_0` gives

     .. math::

        \Lambda_0
        = T_0\Biggl[1 - 1 + \frac{1}{r}\int_0^r \mu(s)\,ds\Biggr]
        = \frac{1}{r} T_0 \int_0^r \mu(s)\,ds.

     Finally, for the linear :math:`K(s)` assumed above,

     .. math::

        \mu(s)
        = \exp\Bigl(\int_0^s K(u)du\Bigr)
        = \exp\Bigl(K_0 s + \frac{K_1-K_0}{2r} s^2\Bigr),

     so we obtain exactly

     .. math::

        \Lambda_0 = \frac{1}{r} T_0 \int_0^r
        \exp\Bigl(K_0 s + \frac{K_1-K_0}{2r} s^2\Bigr)\,ds,

     which is the expression used above (with index :math:`i`):

     .. math::

        \Lambda_i = \frac{1}{r_i} T_i \int_0^{r_i}
        \exp\left(K_i s + \frac{K_{i+1}-K_i}{2 r_i} s^2\right) ds.

     In other words, :math:`\Lambda_i` is precisely the coefficient
     multiplying :math:`(J_i - J_{i+1})` that comes from solving the
     ODE with both :math:`K` and :math:`J` linear in :math:`s`.

   3.3. Dawson function connection. The starting point is

      .. math::

        \Lambda_0 = \frac{1}{r} T_0 \int_0^{r}
        \exp\Bigl(K_0 s + \tfrac{K_1-K_0}{2r} s^2\Bigr)\,ds,

      which is valid both for scalar and matrix-valued :math:`K_0, K_1`. In the
      remainder of this subsection we first treat the **scalar** case (or,
      equivalently, the case where :math:`K_0` and :math:`K_1` commute and can be
      simultaneously diagonalized), because only then can we complete the square and
      express the integral in terms of the (scalar) Dawson function.

      Using :math:`\alpha = (K_1-K_0)/r` as above, the integral

      .. math::

        \int_0^{r} \exp\Bigl(K_0 s + \tfrac12 \,\alpha s^2\Bigr)\,ds

      is a Gaussian-type integral with a linear term in the exponent.
      Completing the square gives

      .. math::

        K_0 s + \tfrac12 \,\alpha s^2
        &= \tfrac12 \,\alpha\Bigl(s^2 + 2\frac{K_0}{\alpha} s\Bigr) \\
        &= \tfrac12 \,\alpha\Bigl\{\bigl(s + \tfrac{K_0}{\alpha}\bigr)^2
        - \bigl(\tfrac{K_0}{\alpha}\bigr)^2\Bigr\}.

      Thus

      .. math::

        \int_0^{r} \exp\Bigl(K_0 s + \tfrac12 \,\alpha s^2\Bigr)ds
        = e^{-\frac{K_0^2}{2\alpha}}
        \int_0^{r}
        \exp\Bigl(\tfrac12 \,\alpha\bigl(s + \tfrac{K_0}{\alpha}\bigr)^2\Bigr)ds.

      With the substitution

      .. math::

        t = \sqrt{\tfrac{\alpha}{2}}\Bigl(s + \tfrac{K_0}{\alpha}\Bigr)
        \quad\Rightarrow\quad
        ds = \sqrt{\tfrac{2}{\alpha}}\,dt,

      the integral becomes

      .. math::

        \sqrt{\tfrac{2}{\alpha}}\, e^{-\frac{K_0^2}{2\alpha}}
        \int_{t_0}^{t_1} e^{t^2}\,dt,

      which can be expressed in terms of the Dawson function

      .. math::

        D(z) = e^{-z^2}\int_0^z e^{u^2}\,du.

      Carrying this through in the scalar case gives the representation

      .. math::

        \Lambda_0
        = \frac{1}{r}\,\sigma^{-1}\bigl(D(u_1) - T_0\,D(u_0)\bigr),

      where

      .. math::

        \sigma &= \sqrt{\frac{K_1-K_0}{2r}}, \\
        u_0 &= \tfrac12 \sigma^{-1} K_0, \\
        u_1 &= \tfrac12 \sigma^{-1} K_1.

      **Extension to matrices.** For *matrix*-valued :math:`K_0,K_1`, the
      exact expression for :math:`\Lambda_0` is still

      .. math::

        \Lambda_0
        = \frac{1}{r}\,T_0 \int_0^{r}
        \exp\Bigl(K_0 s + \tfrac{K_1-K_0}{2r} s^2\Bigr)\,ds,

      but the “complete the square” steps above only go through without
      change if all matrices involved commute (for example, when
      :math:`K_0` and :math:`K_1` are simultaneously diagonalizable).
      Under this *commuting* assumption we may treat the matrices as if
      they were scalars and define matrix analogues

      .. math::

        \boldsymbol{\sigma} = \sqrt{\frac{K_{i+1} - K_i}{2 r_i}}
        \quad\text{(principal matrix square root)},\\
        \mathbf{u}_0 = \tfrac12 \boldsymbol{\sigma}^{-1} K_i, \qquad
        \mathbf{u}_1 = \tfrac12 \boldsymbol{\sigma}^{-1} K_{i+1},

      and a matrix Dawson function :math:`\mathcal{D}` by functional
      calculus applied to these commuting matrices. This leads to the
      formal matrix expression

      .. math::

        \mathbf{\Lambda}_i =
        \frac{1}{r_i}\,\boldsymbol{\sigma}^{-1}
        \Bigl(\mathcal{D}(\mathbf{u}_1) - T_i\,\mathcal{D}(\mathbf{u}_0)\Bigr),

      where all factors on the right-hand side commute because they are
      analytic functions of :math:`K_i` and :math:`K_{i+1}`.

      In practice, ARTS does *not* evaluate a full matrix Dawson
      function. Instead, it applies the scalar Dawson function element-wise
      in a fixed basis, which is only strictly valid if the propagation
      matrix is diagonal (or nearly diagonal) in that basis at both ends of
      the layer. This is why the documentation refers to the Dawson-based
      expression as an approximation in the general polarized case.

   .. admonition:: Implementation note
      :class: tip

      This method requires computing the matrix square root of
      :math:`(K_{i+1}-K_i)/(2 r_i)`.  It also requires using a matrix
      Dawson function.  Neither of these operations are numerically
      stable.  It is therefore not recommended to use this method
      unless you have a good reason to do so.  It is mainly included
      here for completeness, and as an indicator for future improvements
      of the radiative transfer solver in ARTS.
