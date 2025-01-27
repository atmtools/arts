Line-by-line Absorption
#######################

This section describes the physical process of absorption lines of 
different molecules absorbing and emitting spectral radiance
in the atmosphere.

These are the types of line-by-line absorption considered here:


- :ref:`lbl-plain`, where each absorption line is considered separately.

  - Without Zeeman effect
  - With Zeeman effect

- :ref:`lbl-ecs`,
  where the absorption lines of similar energies of a molecule
  are mixed together.

.. _lbl-plain:

Line-by-line Absorption Overview
********************************

The absorption in plain line-by-line absorption is simply the sum of all 
absorption by each absorption line.  The absorption of a single absorption line
is described by the following equations:

.. math::

  \alpha = S(\cdots) F(\cdots),

where
:math:`\alpha` is the absorption coefficient,
:math:`S` is the :ref:`lbl-line-strength` operator, and
:math:`F` is the :ref:`lbl-line-shape` operator.

Both :math:`S` and :math:`F` change slightly if Zeeman effect is considered.
The main way that Zeeman effect changes the calculations is via the polarization
it introduces to the propagation matrix summation.

Without Zeeman effect
=====================

The contribution to :ref:`propagation matrix <prop-mat>` from all non Zeeman-split absorption lines is simply

.. math::

  K_{A, lbl} = \mathrm{Re} \sum_i \alpha_{i},

and from this the full matrix is

.. math::

  \mathbf{K}_{lbl} = \left[ \begin{array}{llll} K_{A, lbl}&0&0&0\\ 0&K_{A, lbl}&0&0\\0&0&K_{A, lbl}&0\\0&0&0&K_{A, lbl} \end{array} \right],

where :math:`i` is the pseudo-index of the absorption line and ``lbl`` is the pseudo-index of the plain line-by-line absorption for the sake of :ref:`summing up absorption <eq-prop-mat-sumup>`.

.. note::

  Plain line-by-line absorption only contribute towards the diagonal of the :ref:`propagation matrix <prop-mat>`.

With Zeeman effect
==================

When Zeeman effect is considered, there are effectively 3 separate kinds of polarized absorption added to the :ref:`propagation matrix <prop-mat>`

.. math::

  K_{\sigma_\pm, z} &= \sum_i \alpha_{i, \sigma_\pm} \\
  K_{\pi, z} &= \sum_i \alpha_{i, \pi}

and from this, the full matrix contribution is

.. math::

  \mathbf{K}_{z} =
  \sum_\pm \left(
  \mathrm{Re} K_{\sigma_\pm,z} \left[
  \begin{array}{rrrr}
  1 + \cos^2\theta_m &
  \sin^2\theta_m\cos 2\eta_m &
  \sin^2\theta_m\sin 2\eta_m &
  \mp 2\cos\theta_m \\
  \sin^2\theta_m\cos 2\eta_m  &
  1 + \cos^2\theta_m  &
  0 &
  0 \\
  \sin^2\theta_m\sin 2\eta_m &
  0 &
  1 + \cos^2\theta_m &
  0 \\
  \mp 2\cos\theta_m &
  0 &
  0 &
  1 + \cos^2\theta_m 
  \end{array}  \right] +
  \mathrm{Im} K_{\sigma_\pm,z} \left[
  \begin{array}{rrrr}
  0 &
  0 &
  0 &
  0 \\
  0 &
  0 &
  \mp 4\cos\theta_m &
  2\sin^2\theta_m\sin 2\eta_m  \\
  0 &
  \pm 4\cos\theta_m  &
  0 &
  - 2 \sin^2\theta_m\cos 2\eta_m \\
  0 &
  - 2\sin^2\theta_m\sin 2\eta_m &
  2 \sin^2\theta_m\cos 2\eta_m &
  0
  \end{array}  \right]
  \right) +
  \\
  \mathrm{Re} K_{\pi,z} \left[
  \begin{array}{rrrr}
  \sin^2\theta_m &
  - \sin^2\theta_m\cos 2\eta_m  &
  - \sin^2\theta_m\sin 2\eta_m  &
  0 \\
  - \sin^2\theta_m\cos 2\eta_m  &
  \sin^2\theta_m  &
  0 &
  0 \\
  - \sin^2\theta_m\sin 2\eta_m  &
  0 &
  \sin^2\theta_m  &
  0 \\
  0 &
  0 &
  0 &
  \sin^2\theta_m 
  \end{array}  \right] +
  \mathrm{Im} K_{\pi,z} \left[
  \begin{array}{rrrr}
  0 &
  0 &
  0 &
  0 \\
  0 &
  0 &
  0 &
  - 2\sin^2\theta_m\sin 2\eta_m \\
  0 &
  0 &
  0 &
  2 \sin^2\theta_m\cos 2\eta_m \\
  0 &
  2\sin^2\theta_m\sin 2\eta_m &
  - 2 \sin^2\theta_m\cos 2\eta_m &
  0
  \end{array}  \right],

where the somewhat weird :math:`\pm`-sum is over the sigma components.
Here the angles :math:`\theta_m` and :math:`\eta_m` are the angles with regards to the magnetic field.

Given a spherical coordinate observation system with zenith angle :math:`\theta_z` and aziumuth angle :math:`\eta_a` and a
local magnetic field with upwards facing strength :math:`B_w`, eastward facing strength :math:`B_u` and northward facing strength :math:`B_v`,
these angles are given by

.. math::

  \theta_m = \arccos\left(\frac{B_v \cos\eta_a \sin\theta_z + B_u \sin\eta_a \sin\theta_z + B_w \cos\theta_z}{ \sqrt{B_w^2 + B_u^2 + B_v^2} } \right)
  \\
  \eta_m = \mathrm{atan2}\left(B_u \cos\eta_a - B_v \sin\eta_a,\; B_w \cos\eta_a\cos\theta_z + B_u\sin\eta_a\cos\theta_z - B_w\sin\theta_z \right)

.. _lbl-line-shape:

Line Shapes
===========

Line shapes should distribute absorption as a function of frequency.
By convention, the line shape is normalized to have an integral of 1.

Voigt Line Shape
----------------

.. math::

  F = \frac{1 + G_{lm} - iY_{lm}}{\sqrt{\pi}G_D} w(z),

where

.. math::

  z = \frac{\nu - \nu_0 - \Delta\nu_{lm} - \Delta_nu_Z - \Delta\nu_{P,0} + iG_{P,0}}{G_D},

where

.. list-table::
  :header-rows: 1

  * - Parameter
    - Description
  * - :math:`\nu`
    - The sampling frequency
  * - :math:`\nu_0`
    - The line center frequency
  * - :math:`G_D`
    - The scaled Doppler broadening half-width half-maximum
  * - :math:`\Delta\nu_Z`
    - The Zeeman shift
  * - :math:`G_{P,0}`
    - The pressure broadening - half width half maximum in the Lorentz profile
  * - :math:`\Delta\nu_{P,0}`
    - The pressure shift
  * - :math:`Y_{lm}`
    - The 1st order Line-mixing parameter
  * - :math:`G_{lm}`
    - The 2nd-order strength modifying line mixing parameter
  * - :math:`\Delta\nu_{lm}`
    - The 2nd-order line-mixing shift
  * - :math:`w(z)`
    - The Faddeeva function.

For more information about how :math:`G_{P,0}`, :math:`\Delta\nu_{P,0}`, :math:`Y_{lm}`, :math:`G_{lm}`, and :math:`\Delta\nu_{lm}` are computed see :ref:`lbl-line-shape-params`.
The scaled Doppler broadening half width half maximum is given by

.. math::

  G_D = \sqrt{\frac{2000 R T}{mc^2}} \nu_0,

where

.. list-table::
  :header-rows: 1

  * - Parameter
    - Description
  * - :math:`R`
    - The ideal gas constant in Joules per mole per Kelvin,
  * - :math:`T`
    - The temperature in Kelvin,
  * - :math:`m`
    - The mass of the molecule in grams per mole, and
  * - :math:`c`
    - The speed of light in meters per second.

The factor 2000 is to convert to SI units.

The Faddeeva function is in ARTS
computed using the MIT-licensed `Faddeeva package <http://ab-initio.mit.edu/faddeeva/>`_,
which is based in large parts on the work by :cite:t:`zaghloul12:_algorithm916_acm`.

The Zeeman line-shift is derived from the magnetic field strength and the magnetic quantum number of the transition.
A linear Zeeman effect is assumed such that

.. math::

  \Delta\nu_Z = \frac{e} {4 \pi m_e} \left(M_l g_{l,z} - M_u g_{u,z} \right),

where :math:`e` is the elementary charge, :math:`m_e` is the mass of an electron, 
:math:`M_l` and :math:`M_u` are the projection of the lower and upper states, respectively, 
of the angular momentum quantum number on the magnetic field, and :math:`g_{l,z}`
and :math:`g_{u,z}` are the lower and upper state Landé g-factors, respectively.
The latter are generally computed ahead of time, e.g., as by :cite:t:`larsson19:_updated_jqsrt,larsson20:_zeeman_jqsrt`.

.. note::

  It is important to not confuse the line-mixing parameters used here with full line mixing as described below.
  The line-mixing paramters here are still plain line-by-line absorption, but it is important that there are no
  cut lines and that the data for *all* line-paramters are derived toghether. 

.. _lbl-line-shape-params:

Line Shape Parameters
---------------------

The line shape parameters supported by ARTS are

.. list-table::
  :header-rows: 1

  * - Parameter
    - Description
    - Pressure Dependency
  * - :math:`G_{P,0}`
    - Pressure broadening half width half maximum, collision-independent.
    - :math:`P`
  * - :math:`G_{P,2}`
    - Pressure broadening half width half maximum, collision-dependent.
    - :math:`P`
  * - :math:`\Delta\nu_{P,0}`
    - Pressure shift, collision-independent.
    - :math:`P`
  * - :math:`\Delta\nu_{P,2}`
    - Pressure shift, collision-dependent.
    - :math:`P`
  * - :math:`\nu_{VC}`
    - Velocity changing frequency.
    - :math:`P`
  * - :math:`\eta`
    - Correlation parameter.
    - :math:`-`
  * - :math:`Y_{lm}`
    - 1st order Line-mixing parameter.
    - :math:`P`
  * - :math:`G_{lm}`
    - 2nd-order strength modifying line mixing parameter.
    - :math:`P^2`
  * - :math:`\Delta\nu_{lm}`
    - 2nd-order line-mixing shift.
    - :math:`P^2`

These parameters are all computed species-by-species before being volume-mixing ratio weighted and summed up.
In equation form:

.. math::

  L = \frac{\sum_i x_i L_i}{\sum_i x_i},

where :math:`L` is a placeholder for any of the line shape parameters, and :math:`x` is the volume mixing ratio, and :math:`i` is a species index.
The normalization is there to allow fewer than all species to contribute to the line shape parameters.

The temperature dependencies of the individual :math:`L_i` are computed based on avaiable data.
There is no general form avaiable, so instead the temperature dependencies are computed based on the data avaiable for each species as:

.. list-table::
  :header-rows: 1

  * - Name
    - Equation
    - Description
  * - ``T0``
    - :math:`L_i(T) = X_0`
    - Constant regardless of temperature
  * - ``T1``
    - :math:`L_i(T) = X_0 \left(\frac{T_0}{T}\right)^{X_1}`
    - Simple power law
  * - ``T2``
    - :math:`L_i(T) = X_0 \left(\frac{T_0}{T}\right) ^ {X_1} \left[1 + X_2 \log\left(\frac{T_0}{T}\right)\right]`
    - Power law with compensation.
  * - ``T3``
    - :math:`L_i(T) = X_0 + X_1 \left(T - T_0\right)`
    - Linear in temperature
  * - ``T4``
    - :math:`L_i(T) = \left[X_0 + X_1 \left(\frac{T_0}{T} - 1\right)\right] \left(\frac{T_0}{T}\right)^{X_2}`
    - Power law with compensation.  Used for line mixing.
  * - ``T5``
    - :math:`L_i(T) = X_0 \left(\frac{T_0}{T}\right)^{\frac{1}{4} + \frac{3}{2}X_1}`
    - Power law with offset.
  * - ``AER``
    - :math:`L_i(200) = X_0`, :math:`L_i(250) = X_1`, :math:`L_i(296) = X_2`, :math:`L_i(340) = X_3`.  Linear interpolation inbetween.
    - Inspired by the way `AER <http://rtweb.aer.com/lblrtm.html>`_ deals with linemixing.
  * - ``DPL``
    - :math:`L_i(T) = X_0 \left(\frac{T_0}{T}\right) ^ {X_1} + X_2 \left(\frac{T_0}{T}\right) ^ {X_3}`
    - Double power law.
  * - ``POLY``
    - :math:`L_i(T) = X_0 + X_1 T + X_2 T ^ 2 + X_3 T ^ 3 + \cdots`
    - Polynomial in temperature.  Used internal in ARTS when training our own linemixing.

here, :math:`X_0` ... :math:`X_N` are all model supplied constants whereas :math:`T` is the temperature in Kelvin and :math:`T_0` is the reference temperature of the model parameters.

.. _lbl-line-strength:

Line Strength
=============

.. _lbl-lte:

Local thermodynamic equilibrium
-------------------------------

For local thermodynamic equilibrium (LTE), the line strength is given by

.. math::

  S_{LTE} = \rho \frac{c^2\nu}{8\pi} \left[1 - \exp\left(-\frac{h\nu}{kT}\right)\right]
  \frac{g_u\exp\left(-\frac{E_l}{kT}\right)}{Q(T)} \frac{A_{lu}}{\nu_0^3}

.. _lbl-nlte:

Non-local thermodynamic equilibrium
-----------------------------------

For non-LTE, the line strength is given by

.. math::

  S_{NLTE} = \rho \frac{c^2\nu}{8\pi} \left(r_l \frac{g_u}{g_l} - r_u\right) \frac{A_{lu}} {\nu_0^3},

and the added emissions are given by

.. math::

  K_{NLTE} = \rho \frac{c^2\nu}{8\pi} \left\{r_u\left[
  1 - \exp\left(\frac{h\nu_0}{kT}\right)\right] - \left(r_l \frac{g_u}{g_l} - r_u\right)
  \right\} \frac{ A_{lu}}{\nu_0^3},

where :math:`r_l` and :math:`r_u` are the ratios of the populations of the lower and upper states, respectively.
Note that :math:`K_{LTE} = 0`, as it represents "additional" emission due to non-LTE conditions.
Also note that :math:`K_{NLTE}` may be negative.

To ensure ourselves that this can be turned into the expression for LTE,
we can rewrite the above for the expression that :math:`r_l` and :math:`r_u`
would have in LTE according to the Boltzmann distribution:

.. math::

  r_l = \frac{g_l\exp\left(-\frac{E_l}{kT}\right)}{Q(T)}

and

.. math::

  r_u = \frac{g_u\exp\left(-\frac{E_u}{kT}\right)}{Q(T)}

Putting this into the ratio-expression for :math:`S_{NLTE}` with the following simplification steps:

Expansion:

.. math::

  \left(r_l \frac{g_u}{g_l} - r_u\right) =
  \frac{g_u}{Q(T)}\left[\exp\left(-\frac{E_l}{kT}\right) - \exp\left(-\frac{E_u}{kT}\right)\right].

Extract lower state energies:

.. math::

  \frac{g_u}{Q(T)}\left[\exp\left(-\frac{E_l}{kT}\right) - \exp\left(-\frac{E_u}{kT}\right)\right]
  \frac{\exp\left(-\frac{E_l}{kT}\right)}{\exp\left(-\frac{E_l}{kT}\right)} \rightarrow
  \left[1 - \exp\left(-\frac{h\nu_0}{kT}\right)\right]\frac{g_u\exp\left(-\frac{E_l}{kT}\right)}{Q(T)},

where this last step is possible because we estimate that :math:`E_u-E_l = h\nu_0`.  Note how the
expression for :math:`K_{NLTE}` is 0 under LTE conditions. As it should be.
This is seen by putting the above RHS and the expression for :math:`r_u` into the expression for :math:`K_{NLTE}`:

.. math::

  K_{NLTE} = \rho \frac{c^2\nu}{8\pi} \left\{\frac{g_u\exp\left(-\frac{E_u}{kT}\right)}{Q(T)}\left[
    1 - \exp\left(\frac{h\nu_0}{kT}\right)\right] - \left[1 - \exp\left(-\frac{h\nu_0}{kT}\right)\right]\frac{g_u\exp\left(-\frac{E_l}{kT}\right)}{Q(T)}
    \right\} \frac{ A_{lu}}{\nu_0^3} = 0.

The ratio between LTE and non-LTE line strength remaining is:

.. math::

  \frac{S_{NLTE}}{S_{LTE}} = \frac{1 - \exp\left(-\frac{h\nu_0}{kT}\right)}{1 - \exp\left(-\frac{h\nu}{kT}\right)}.

It is clear that the non-LTE expression is the one that is incorrect here.
The energy of the emitted photon is not :math:`h\nu_0` but :math:`h\nu`, and
as such the actual energy of the transition is :math:`E'_u-E'_l = h\nu`, but
this should be relatively close in cases where we actually care about non-LTE
(which is low density, low collision atmospheres).

Zeeman effect
-------------

If Zeeman effect is considered, the emission and absorption terms above are modified by quantum number state distribution.
For O :sub:`2`, for example, this introduces a factor of

.. math::

  S_z = f(\Delta M) \left( \begin{array}{ccc} J_l & 1 & J_u \\ M_l & \Delta M & - M_u \end{array} \right)^2,

where :math:`\Delta M \in \left[-1,\;0,\;1\right]` is the change in quantum number for angular rotational momentum projection along the magnetic field
for :math:`\sigma_-`, :math:`\pi`, and :math:`\sigma_+`, respectively,
:math:`f(\Delta M)` is the 0.75 for :math:`\sigma_\pm` and 1.5 for :math:`\pi`,
and :math:`J_l` and :math:`J_u` are the lower and upper total angular rotational momentum quantum number.
The :math:`(:::)` construct is the Wigner 3-j symbol.
It can be `computed using software <https://fy.chalmers.se/subatom/wigxjpf/>`_ such as that by :cite:t:`johansson2016`.

.. _lbl-ecs:

Line-mixing using Error-corrected Sudden
****************************************

TBD
