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
     \mathbf{K}_{z} =\sum_\pm\left(\\
     \mathrm{Re} K_{\sigma_\pm,z} \left[\begin{array}{llll}
          1 + \cos^2\theta_m           &  \sin^2\theta_m \cos 2 \eta_m & -\sin^2 \theta_m \sin 2 \eta_m & \mp 2 \cos \theta_m  \\
          \sin^2\theta_m \cos 2 \eta_m &  1 + \cos^2 \theta_m          &  0                             &     0                \\
         -\sin^2\theta_m \sin 2 \eta_m &  0                            &  1 + \cos^2 \theta_m           &     0                \\
         \mp 2 \cos\theta_m            &  0                            &  0                             &     1 + \cos^2 \theta_m
      \end{array}\right] +
       \mathrm{Im} K_{\sigma_\pm,z} \left[\begin{array}{llll}
          0 &  0                           &      0                      &  0                             \\
          0 &  0                           &  \pm 2 \cos \theta_m        & -\sin^2 \theta_m \sin 2 \eta_m \\
          0 & \mp 2 \cos \theta_m          &      0                      & -\sin^2 \theta_m \cos 2 \eta_m \\
          0 & \sin^2\theta_m \sin 2 \eta_m &  \sin^2\theta_m \cos 2 \eta_m &  0
      \end{array}\right] \right)  +\\
     \mathrm{Re} K_{\pi,z} \left[\begin{array}{llll}
         \sin^2\theta_m               &  -\sin^2\theta_m \cos 2 \eta_m &  \sin^2 \theta_m \sin 2 \eta_m &   0 \\
        -\sin^2\theta_m \cos 2 \eta_m &   \sin^2\theta_m               &  0                             &   0 \\
         \sin^2\theta_m \sin 2 \eta_m &   0                            &  \sin^2\theta_m                &   0 \\
         0                            &   0                            &  0                             &  \sin^2\theta_m
      \end{array}\right] +
       \mathrm{Im} K_{\pi,z} \left[ \begin{array}{llll}
          0 &  0                            &      0                         &  0                           \\
          0 &  0                            &      0                         & \sin^2 \theta_m \sin 2 \eta_m \\
          0 &  0                            &      0                         & \sin^2 \theta_m \cos 2 \eta_m \\
          0 & -\sin^2\theta_m \sin 2 \eta_m &  -\sin^2\theta_m \cos 2 \eta_m &  0
      \end{array}\right],

where the somewhat weird :math:`\pm`-sum is over the sigma components.
Here the angles :math:`\theta_m` and :math:`\eta_m` are the angles with regards to the magnetic field.

Given a spherical coordinate observation system with zenith angle :math:`\theta_z` and aziumuth angle :math:`\eta_a` and a
local magnetic field with upwards facing strength :math:`B_w`, eastward facing strength :math:`B_u` and northward facing strength :math:`B_v`,
these angles are given by

.. math::

  \theta_m = \arccos\left(\frac{B_v \cos\eta_a \sin\theta_z + B_u \sin\eta_a \sin\theta_z + B_w \cos\theta_z}{ \sqrt{B_w^2 + B_u^2 + B_v^2} } \right)
  \\
  \eta_m = -\mathrm{atan2}\left(B_u \cos \eta_a-B_v \sin\eta_a,\; B_u \cos\theta_z \sin\eta_a + B_v \cos\theta_z\cos\eta_a - B_w\sin\theta_z \right)

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

where :math:`L` is a placeholder for any of the line shape parameters, and :math:`x` is the volume-mixing ratio, and :math:`i` is a species index.
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
  \frac{g_u\exp\left(-\frac{E_l}{kT}\right)}{Q(T)} \frac{A_{lu}}{\nu_0^3},

where :math:`\rho` is the number density of the absorbing species,

.. math::

  \rho = \mathrm{VMR}\frac{P}{kT},

where VMR is the volume-mixing ratio of the absorbing species.

.. list-table::
  :header-rows: 1

  * - Parameter
    - Description
  * - :math:`\rho`
    - Number density of the absorbing isotopologue, :math:`\mathrm{VMR} \cdot P / (kT)`
  * - :math:`\mathrm{VMR}`
    - Volume-mixing ratio of the absorbing species
  * - :math:`P`
    - Atmospheric pressure
  * - :math:`c`
    - Speed of light
  * - :math:`\nu`
    - Sampling frequency
  * - :math:`\nu_0`
    - Line centre frequency
  * - :math:`h`
    - Planck constant
  * - :math:`k`
    - Boltzmann constant
  * - :math:`T`
    - Temperature
  * - :math:`g_u`
    - Degeneracy of the upper state
  * - :math:`E_l`
    - Energy of the lower state
  * - :math:`Q(T)`
    - Partition function at temperature :math:`T`
  * - :math:`A_{lu}`
    - Einstein A coefficient for spontaneous emission

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
For O\ :sub:`2`, for example, this introduces a factor of

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

When the atmosphere is at sufficient pressure, collisions occur frequently enough that
absorption lines of a vibrational-rotational band can no longer be treated
independently.  Molecules undergoing collisions may exchange rotational angular
momentum, transferring population between rotational levels.  At intermediate
pressures this introduces off-diagonal couplings between lines in a spectral band,
leading to the phenomenon of *line mixing*.  At high pressures the lines collapse
towards a pressure-broadened Q-branch.

The Error-corrected Sudden (ECS) approximation provides a rigorous, quantum-mechanical
framework for computing this mixing.  The "Sudden" part refers to the
Infinite-Order-Sudden (IOS) approximation, in which the collision time is assumed short
compared with the rotational period.  The "Error-corrected" part refers to the
subsequent rescaling of the relaxation matrix elements to satisfy the first-order
optical sum rule exactly, removing a systematic bias that arises from the sudden
approximation.

.. _lbl-ecs-lineshape:

ECS Line Shape
==============

For a band of :math:`n` interacting absorption lines, the ECS complex absorption shape
for a single broadening species (see :ref:`lbl-ecs-multispecies` for the full
expression) is written in terms of the complex relaxation matrix :math:`\mathbf{W}` as

.. math::

  \chi(\nu) \propto \mathrm{Im}\left[\mathbf{d}^T \left(\nu \mathbf{I} - \mathbf{W}\right)^{-1} \mathbf{p}\, \mathbf{d}\right],

where :math:`\mathbf{d}` is the vector of reduced dipole matrix elements (one entry per
line), :math:`\mathbf{p}` is the diagonal matrix of lower-state thermal populations
(:math:`p_j = g_j\exp(-E_l^{(j)}/kT)/Q(T)`, the Boltzmann fractional population of
the lower state of line :math:`j`, following the :ref:`lbl-lte` notation where
:math:`g_j` is the lower-state degeneracy and :math:`Q(T)` the partition function),
and :math:`\mathbf{W}` is the full complex :math:`n \times n` relaxation matrix whose
diagonal and off-diagonal elements are described in :ref:`lbl-ecs-relaxmat`.

The relaxation matrix is diagonalised as :math:`\mathbf{W} = \mathbf{V} \tilde{\boldsymbol{\nu}} \mathbf{V}^{-1}`,
where :math:`\tilde{\boldsymbol{\nu}}` is the diagonal matrix of complex
*equivalent line* positions.  Each equivalent line :math:`k` has a complex
frequency :math:`\tilde{\nu}_k` (real part: position, imaginary part: pressure
broadening) and a complex *equivalent strength* :math:`\tilde{S}_k`.

Explicitly, the equivalent strength for line :math:`k` is

.. math::

  \tilde{S}_k = \left(\sum_j d_j V_{jk}\right) \left(\sum_j p_j d_j V^{-1}_{kj}\right),

and the ECS line shape function is

.. math::

  F_{ECS}(\nu) = \frac{\sqrt{\ln 2}}{\sqrt{\pi}} \sum_k \tilde{S}_k \frac{w(z_k)}{G_{D,k}},

where :math:`w` is the Faddeeva function and

.. math::

  z_k = \frac{\left(\tilde{\nu}_k - \nu\right)\sqrt{\ln 2}}{G_{D,k}}, \qquad
  G_{D,k} = G_D^{fac} \cdot \mathrm{Re}\!\left[\tilde{\nu}_k\right],

with the Doppler scale factor

.. math::

  G_D^{fac} = \sqrt{\frac{2000 R T}{m c^2}},

where :math:`R` is the ideal gas constant in J mol\ :sup:`-1` K\ :sup:`-1`,
:math:`m` the molar mass in g mol\ :sup:`-1`, :math:`c` the speed of light, and
:math:`T` the temperature (same symbols as in the plain LBL definition in
:ref:`lbl-line-shape`).

Note that :math:`G_{D,k}` is formed by multiplying :math:`G_D^{fac}` by
:math:`\mathrm{Re}[\tilde{\nu}_k]` — the real part of the :math:`k`-th eigenvalue
— rather than by any original line centre :math:`\nu_{0,j}`.  This is necessary
because eigenvalue decomposition does not in general return eigenvalues in the same
order as the input lines, so there is no well-defined mapping from equivalent line
:math:`k` to a single physical line :math:`j`.  Using :math:`\mathrm{Re}[\tilde{\nu}_k]`
keeps the Doppler width self-consistent with the actual position of each equivalent
line.

The contribution to the :ref:`propagation matrix <prop-mat>` from the entire band is
then

.. math::

  K_{A, ecs} = N \nu \left(1 - \exp\!\left(-\frac{h\nu}{kT}\right)\right) \mathrm{Re}\!\left[F_{ECS}(\nu)\right],

where :math:`N` is the total number density of the absorbing species.

.. note::

  Zeeman splitting within a band is not currently supported together with ECS
  line mixing.

.. _lbl-ecs-relaxmat:

Relaxation Matrix
=================

The complex relaxation matrix :math:`\mathbf{W}` has dimensions
:math:`n \times n` (lines in the band) and is constructed as follows.

The *real* part of the diagonal elements carries the (pressure-shifted) line
centre frequencies:

.. math::

  \mathrm{Re}\, W_{ii} = \nu_{0,i} + \Delta\nu_{P,0,i},

where :math:`\nu_{0,i}` is the vacuum line centre and :math:`\Delta\nu_{P,0,i}` is the
pressure shift of line :math:`i`.

The *imaginary* part of the diagonal elements carries the pressure broadening:

.. math::

  \mathrm{Im}\, W_{ii} = G_{P,0,i},

where :math:`G_{P,0,i}` is the pressure-broadening half-width half-maximum.

The off-diagonal elements :math:`W_{ij}` (:math:`i \neq j`) encode the rate of
transfer from line :math:`j` to line :math:`i` via collisions.  Their construction
from the ECS basis rates is described below.

.. _lbl-ecs-rates:

ECS Basis Rates
===============

The ECS approach introduces two species-dependent functions of the integer angular
momentum transfer channel :math:`L`:

**Basic rate** :math:`Q(L)`:
  This encodes the intrinsic probability of a collision transferring :math:`L` units of
  angular momentum.  Its temperature dependence is parameterised as

  .. math::

    Q(L, T) = s(T) \cdot \frac{e^{-\beta(T)\, E_L / kT}}{[L(L+1)]^{\lambda(T)}},

  where :math:`E_L` is the rotational energy of level :math:`L`,
  and :math:`s(T)`, :math:`\beta(T)`, and :math:`\lambda(T)` are
  temperature-dependent model parameters stored per broadening species
  (see :ref:`lbl-line-shape-params` for the available temperature dependence forms).

**Adiabatic factor** :math:`\Omega(L)`:
  The IOS approximation becomes inaccurate when the rotational period approaches
  the collision duration.  The adiabatic factor corrects for this using the
  coupling model:

  .. math::

    \Omega(L) = \frac{1}{\left[1 + \dfrac{\omega_{L,L-2}^2\, \tau_c^2}{24}\right]^2},

  where :math:`\omega_{L,L-2} = (E_L - E_{L-2})/\hbar` is the angular frequency of the
  :math:`L \to L-2` rotational transition and

  .. math::

    \tau_c = \frac{\sigma_c(T)}{\bar{v}}, \qquad
    \bar{v} = \sqrt{\frac{8kT}{\pi\mu}},

  with :math:`\sigma_c(T)` the temperature-dependent mean collisional diameter (a
  per-species model parameter, distinct from the reduced dipole :math:`d_j`),
  :math:`\mu` the reduced mass of the colliding pair, and :math:`\bar{v}` the
  mean relative thermal speed.

.. _lbl-ecs-offdiag:

Off-diagonal Elements
=====================

The off-diagonal elements of :math:`\mathbf{W}` are computed species-by-species.
All variants follow the formal IOS structure: they are written as a sum over
even angular momentum transfer channels :math:`L`, weighted by the ratio
:math:`Q(L)/\Omega(L)` and by Wigner 3-j and 6-j coupling coefficients.
After the IOS computation an error-correction (sum-rule rescaling) is applied
(see :ref:`lbl-ecs-sumrule`).

Detailed balance is enforced throughout: the rate of transfer from a lower strength
line :math:`j` to a higher strength line :math:`i` is obtained from the downward
rate :math:`W_{ij}` via

.. math::

  W_{ji} = W_{ij} \exp\!\left(\frac{E_j - E_i}{kT}\right),

where :math:`E_i` is the energy of the lower rotational state of line :math:`i`
(using the same :math:`E_l` convention as in :ref:`lbl-lte`).

The four variants implemented in ARTS, corresponding to the four line shape model
types ``VP_ECS_HARTMANN``, ``VP_ECS_MAKAROV``, ``VP_ECS_STOTOP``, and ``VP_ECS_SPHTOP``,
are described below.

Linear Molecules — Hartmann (CO\ :sub:`2`)
------------------------------------------

For linear molecules (e.g. CO\ :sub:`2`) the off-diagonal rate from line :math:`j`
(upper/lower rotational quantum numbers :math:`J'_i, J'_f`, vibrational angular
momentum :math:`l`) to line :math:`i` (:math:`J_i, J_f`) is
:cite:p:`NIRO2004483`

.. math::

  W_{ij} =
    \Omega(J_i)\, (2J'_i+1)\sqrt{(2J_f+1)(2J'_f+1)}
    \sum_L (2L+1)
    \begin{pmatrix} J_i & J'_i & L \\ l & -l & 0 \end{pmatrix}
    \begin{pmatrix} J_f & J'_f & L \\ l & -l & 0 \end{pmatrix}
    \begin{Bmatrix} J_i & J_f & 1 \\ J'_f & J'_i & L \end{Bmatrix}
    \frac{Q(L)}{\Omega(L)},

where the sum runs over even :math:`L \geq \max(|J_i-J'_i|, |J_f-J'_f|)`,
:math:`(\,\cdots)` denotes a Wigner 3-j symbol,
and :math:`\{\,\cdots\}` a Wigner 6-j symbol.

The corresponding reduced dipole element used in the equivalent-strength
calculation is

.. math::

  d(J_f, J_i) = (-1)^{J_f + l_f + 1} \sqrt{2J_f+1}\;
    \begin{pmatrix} J_f & 1 & J_i \\ l_i & l_f - l_i & -l_f \end{pmatrix}.

The rotational energy entering :math:`Q` and :math:`\Omega` is the rigid-rotor
expression :math:`E_J = B_0 J(J+1)`, where :math:`B_0` is the effective
ground-state rotational constant for the species.  Energy levels provided by
quantum-chemical calculations are used directly where available; the rigid-rotor
expression serves to extrapolate to levels not covered by those calculations.

Symmetric Tops with Electron Spin — Makarov (O\ :sub:`2`)
----------------------------------------------------------

Molecular oxygen (O\ :sub:`2`) has an unpaired electron spin :math:`S = 1`, so each
rotational quantum number :math:`N` gives rise to a triplet :math:`J = N-1, N, N+1`.
The off-diagonal coupling between lines :math:`(i: N_l J_l \to N_u J_u)` and
:math:`(j: N'_l J'_l \to N'_u J'_u)` is (using :math:`l`/:math:`u` for the
lower/upper state of each transition, matching the convention of :ref:`lbl-lte`)
:cite:p:`Makarov2020`

.. math::

  W_{ij} =&
    (-1)^{J'_l + J_l + 1}\,
    [N_l][N_u][N'_u][N'_l][J_u][J'_u][J_l][J'_l] \Omega(N_l) \\ &
    \begin{array}{llll}
      \sum_L (2L+1) &
      \begin{pmatrix} N'_l & N_l & L \\ 0 & 0 & 0 \end{pmatrix} &
      \begin{pmatrix} N'_u & N_u & L \\ 0 & 0 & 0 \end{pmatrix} \\ &
      \begin{Bmatrix} L & J_l & J'_l \\ S & N'_l & N_l \end{Bmatrix} &
      \begin{Bmatrix} L & J_u & J'_u \\ S & N'_u & N_u \end{Bmatrix} &
      \begin{Bmatrix} L & J_l & J'_l \\ 1 & J'_u & J_u \end{Bmatrix}
      \frac{Q(L)}{\Omega(L)},
    \end{array}

where :math:`[X] \equiv \sqrt{2X+1}`, :math:`(\cdots)` denotes a Wigner 3-j symbol,
and :math:`\{\cdots\}` a Wigner 6-j symbol.

The reduced dipole is

.. math::

  d(J_u, J_l, N) = (-1)^{J_l + N}
    \sqrt{6(2J_l+1)(2J_u+1)}
    \begin{Bmatrix} 1 & 1 & 1 \\ J_l & J_u & N \end{Bmatrix}.

The rotational energy for the O\ :sub:`2` microwave band is computed from the
full ground-state Hamiltonian including spin–rotation coupling and magnetic
interactions.

Symmetric Tops (NH\ :sub:`3`, PH\ :sub:`3`)
-------------------------------------------

This part is mostly untested and may be incorrect.
It has been generated by AI and is available in ARTS
only for experimentation to see if it produces reasonable results.

For symmetric top molecules (e.g. NH\ :sub:`3`, PH\ :sub:`3`) with :math:`\Delta K = 0`
collisions, lines within the same :math:`K` sub-band are coupled identically to the
Hartmann linear-molecule formula with the vibrational angular momentum :math:`l`
replaced by :math:`K`: :cite:p:`Hadded2002`

.. math::

  W_{ij} =
    \Omega(J_i)\, (2J'_i+1)\sqrt{(2J_f+1)(2J'_f+1)}
    \sum_L (2L+1)
    \begin{pmatrix} J_i & J'_i & L \\ K & -K & 0 \end{pmatrix}
    \begin{pmatrix} J_f & J'_f & L \\ K & -K & 0 \end{pmatrix}
    \begin{Bmatrix} J_i & J_f & 1 \\ J'_f & J'_i & L \end{Bmatrix}
    \frac{Q(L)}{\Omega(L)}.

Lines with different :math:`K` are not coupled. The reduced dipole is

.. math::

  d(J_f, J_i, K) = (-1)^{J_f + K + 1}\sqrt{2J_f+1}\;
    \begin{pmatrix} J_f & 1 & J_i \\ K & 0 & -K \end{pmatrix}.

Rotational energy levels provided by quantum-chemical calculations are used
directly where available; levels beyond those are extrapolated using the
rigid-rotor expression :math:`E_J = B_0 J(J+1)` with the species-specific
ground-state rotational constant :math:`B_0`.

Spherical Tops (CH\ :sub:`4`)
-----------------------------

This part is mostly untested and may be incorrect.
It has been generated by AI and is available in ARTS
only for experimentation to see if it produces reasonable results.

For spherical top molecules (e.g. CH\ :sub:`4`) the coupling reduces to the
:math:`l = 0` limit of the Hartmann formula: :cite:p:`Pieroni1999`

.. math::

  W_{ij} =
    \Omega(J_i)\, (2J'_i+1)\sqrt{(2J_f+1)(2J'_f+1)}
    \sum_L (2L+1)
    \begin{pmatrix} J_i & J'_i & L \\ 0 & 0 & 0 \end{pmatrix}
    \begin{pmatrix} J_f & J'_f & L \\ 0 & 0 & 0 \end{pmatrix}
    \begin{Bmatrix} J_i & J_f & 1 \\ J'_f & J'_i & L \end{Bmatrix}
    \frac{Q(L)}{\Omega(L)},

and the reduced dipole is

.. math::

  d(J_f, J_i) = (-1)^{J_f+1}\sqrt{2J_f+1}\;
    \begin{pmatrix} J_f & 1 & J_i \\ 0 & 0 & 0 \end{pmatrix}.

Rotational energy levels provided by quantum-chemical calculations are used
directly where available; levels beyond those are extrapolated using the
rigid-rotor expression :math:`E_J = B_0 J(J+1)` with the species-specific
ground-state rotational constant :math:`B_0`.

.. _lbl-ecs-sumrule:

Sum-rule Correction
===================

The pure IOS matrix elements computed above do not, in general, satisfy the
first-order optical sum rule exactly due to the finite range of the :math:`L` sum and
the approximate nature of the adiabatic factor.  The *error-correction* step
rescales each column of off-diagonal elements to enforce

.. math::

  \sum_j d_j\, W_{ji} = 0 \quad \forall\, i.

This is done as follows.  For each line :math:`i`, partition the off-diagonal
elements into those coupling to lines with lower intensity-weighted frequency
(summed into :math:`s_\downarrow`) and those coupling to higher-frequency lines
(:math:`s_\uparrow`):

.. math::

  s_\downarrow = \sum_{j > i} d_j\, W_{ji}, \qquad
  s_\uparrow   = \sum_{j < i} d_j\, W_{ji}.

All downward-coupling elements are then rescaled by :math:`-s_\uparrow / s_\downarrow`,
and the corresponding upward-coupling elements are updated by detailed balance.

This rescaling constitutes the "error-corrected" part of the ECS method and ensures
the resulting relaxation matrix produces physically consistent absorption profiles
that recover the correct integrated line intensity at all pressures.

.. _lbl-ecs-multispecies:

Multiple Broadening Species
===========================

When multiple broadening species are present, the ARTS implementation offers two
modes:

In the **single-W mode** (used by default when calling ``calculate``), the per-species
relaxation matrices are first volume-mixing-ratio weighted and summed into a single
effective :math:`\mathbf{W}`:

.. math::

  \mathbf{W}_{eff} = \sum_s x_s\, \mathbf{W}^{(s)},

and a single diagonalisation is performed.

In the **multi-W mode** (used by ``equivalent_values`` for pre-computing equivalent
lines at multiple temperatures), the diagonalisation is performed separately per
species and the resulting absorption contributions are VMR-weighted and summed.

.. _lbl-ecs-rosenkranz:

Rosenkranz Approximation
========================

The full ECS calculation requires the diagonalisation of an :math:`n \times n`
complex matrix at every temperature and pressure of interest, together with a
VMR-weighted sum over broadening species.  For many practical applications a
simpler representation is desirable: the *Rosenkranz approximation* retains the
ordinary Voigt line shape of each line but adds pressure-dependent first- and
second-order correction terms that encode the effect of line mixing to a given
order in pressure.

The corrected Voigt line shape for line :math:`i` is exactly the :ref:`Voigt
profile <lbl-line-shape>` already described,

.. math::

  F_i = \frac{1 + G_{lm,i} - iY_{lm,i}}{\sqrt{\pi}\,G_D}\,w(z_i),

where :math:`z_i` contains :math:`\Delta\nu_{lm,i}` as an additional shift, and
the three correction parameters are:

.. list-table::
  :header-rows: 1

  * - Parameter
    - Physical meaning
    - Pressure scaling
  * - :math:`Y_{lm,i}`
    - First-order line-mixing: asymmetric intensity redistribution between nearby lines.
    - :math:`P`
  * - :math:`G_{lm,i}`
    - Second-order strength correction: quadratic-in-pressure modification to the
      integrated area.
    - :math:`P^2`
  * - :math:`\Delta\nu_{lm,i}`
    - Second-order frequency shift: quadratic-in-pressure displacement of the line
      centre due to the mixing.
    - :math:`P^2`

The Rosenkranz parameters are not fitted to measured spectra directly; instead they
are derived from the ECS equivalent lines or, equivalently, from the relaxation matrix
itself via perturbation theory — both approaches are described below.

.. _lbl-ecs-rosenkranz-W:

Perturbation Theory from the Relaxation Matrix
-----------------------------------------------

When the off-diagonal elements of :math:`\mathbf{W}` are small compared with the
spacings between line centres (the "weak coupling" limit, valid for resolved lines or
moderate pressures), the Rosenkranz parameters can be obtained analytically by
expanding the resolvent :math:`(\nu\mathbf{I} - \mathbf{W})^{-1}` in powers of the
off-diagonal part.  This is the original approach of :cite:t:`rosenkranz:75`.

Write :math:`\mathbf{W} = \mathbf{D} + \mathbf{V}`, where :math:`\mathbf{D}` is the
diagonal part (line centres plus pressure broadening) and :math:`\mathbf{V}_{ij} =
W_{ij}` for :math:`i \neq j` (the off-diagonal relaxation rates, purely imaginary in
the ARTS convention: :math:`V_{ij} = i R_{ij}` with :math:`R_{ij}` real and
proportional to :math:`P`).  Let :math:`g_i(\nu) = [\nu - W_{ii}]^{-1}` be the
unperturbed resolvent for line :math:`i`.  The Neumann expansion then gives

.. math::

  (\nu\mathbf{I} - \mathbf{W})^{-1} =
    \mathbf{G}_0 + \mathbf{G}_0 \mathbf{V} \mathbf{G}_0
    + \mathbf{G}_0 \mathbf{V} \mathbf{G}_0 \mathbf{V} \mathbf{G}_0 + \cdots,

where :math:`\mathbf{G}_0 = \mathrm{diag}(g_i(\nu))`.

Collecting all contributions to the absorption of line :math:`i` through first and
second order, and evaluating the slowly varying factors involving other lines :math:`j`
at :math:`\nu = \nu_i`, yields the three Rosenkranz parameters for line :math:`i`:

**First-order mixing parameter** (:math:`Y_i \sim P`):

.. math::

  Y_i = \frac{2}{S_i} \sum_{j \neq i} S_j \frac{R_{ij}}{\nu_i - \nu_j},

where :math:`S_i = p_i d_i^2` is proportional to the LBL line strength of line
:math:`i` (see :ref:`lbl-ecs-lineshape` for the definition of :math:`p_i`),
and
:math:`R_{ij} = \mathrm{Im}[W_{ij}] / P` is the pressure-normalised off-diagonal
relaxation rate (transfer from line :math:`j` to line :math:`i`; note
:math:`\mathrm{Re}[W_{ij}] = 0` for :math:`i \neq j`), and the sum is over
all other lines :math:`j` in the band.

**Second-order strength correction** (:math:`G_i \sim P^2`):

From the squared first-order cross terms, the fractional modification to the
integrated area of line :math:`i` is

.. math::

  G_i = -\frac{1}{S_i} \sum_{j \neq i} S_j \left(\frac{R_{ij}}{\nu_i - \nu_j}\right)^2.

**Second-order line-centre shift** (:math:`\Delta\nu_i \sim P^2`):

The diagonal self-energy correction (virtual transition :math:`i \to j \to i`) gives

.. math::

  \Delta\nu_i = -\frac{1}{S_i} \sum_{j \neq i} S_j \frac{R_{ij}^2}{\nu_i - \nu_j},

where detailed balance (:math:`S_i R_{ji} = S_j R_{ij}`) has been used to express
everything in terms of the downward rate :math:`R_{ij}`.

.. note::

  Note that :math:`\Delta\nu_i = G_i (\nu_i - \nu_j)` only for a single interfering
  line.  In general they have different frequency denominators (:math:`\nu_i - \nu_j`
  vs :math:`(\nu_i - \nu_j)^2`) and thus differ quantitatively when multiple lines
  contribute.  The two parameters are both needed to correctly reproduce the
  second-order pressure dependence of the band profile.

  The perturbation theory expressions above assume the off-diagonal elements are
  small relative to the line spacing.  They break down for overlapping lines
  (e.g., at very high pressures or for lines very close in frequency).  In that
  regime the full ECS calculation should be used instead.

.. _lbl-ecs-rosenkranz-fitting:

Fitting from Equivalent Lines
------------------------------

The perturbation theory expressions above are analytically exact in the weak-coupling
limit, but in practice it is often more accurate to extract the Rosenkranz parameters
*numerically* from the ECS equivalent lines, because the equivalent-line calculation
already includes the full resummation of the relaxation matrix (not just the first few
terms of the Neumann series).  The two approaches agree at low pressure but the
equivalent-line fit is preferred at higher pressures where the perturbation series
converges slowly.

Given the complex equivalent lines :math:`(\tilde{S}_{k,s}, \tilde{\nu}_{k,s})`
(indexed by :math:`k` in eigenvalue-decomposition order, which carries no physical
meaning) computed by ECS for broadening species :math:`s` at pressure :math:`P_0`
and a grid of temperatures :math:`T_1, \ldots, T_M`, the Rosenkranz coefficients are
extracted by ``abs_bandsLineMixingAdaptation`` as follows.
Here, :math:`i` will denote the index of the physical LBL lines (the rows/columns of
:math:`\mathbf{W}`) and :math:`k` the index of the equivalent lines.

**Step 1 — Sort and match.**
The :math:`n` equivalent lines are sorted by :math:`\mathrm{Re}[\tilde{\nu}_{k,s}]`
and the :math:`n` physical LBL lines are sorted by :math:`\nu_{0,i}`.  Equivalent line
at sorted position :math:`n` is then identified with the physical line at sorted
position :math:`n`, giving a bijection :math:`k \leftrightarrow i` between the two
index sets.  This is a heuristic matching: because eigenvalue decomposition does not
guarantee any particular ordering of eigenvalues, sorting is the only way to establish a
correspondence.  The identification is reliable as long as the second-order frequency
shifts and pressure shifts remain small compared with the separations between adjacent
line centres; it can fail if either effect is comparable in magnitude to the line spacing.

**Step 2 — Form normalised differences.**
For each matched pair :math:`(k, i)`, the LTE line strength of physical line :math:`i`
at temperature :math:`T` is :math:`S_{LTE,i}(T)` as defined in :ref:`lbl-lte`
(without the number density factor :math:`\rho`; the equivalent per-molecule strength is
:math:`s_i(T) = S_{LTE,i}(T)/\rho`).  The unperturbed complex frequency of physical
line :math:`i` is

.. math::

  \nu_i^{LBL}(T) = \nu_{0,i} + \Delta\nu_{P,0,i}(T,P_0) + i\,G_{P,0,i}(T,P_0).

The residual strength ratio and frequency residual are then formed:

.. math::

  r_{i,s}(T) &= \frac{\tilde{S}_{k,s}(T)}{S_{LTE,i}(T)}, \\
  \delta\nu_{i,s}(T) &= \tilde{\nu}_{k,s}(T) - \nu_i^{LBL}(T),

where :math:`k` is the equivalent line matched to physical line :math:`i` in Step 1.

**Step 3 — Extract pressure-normalised Rosenkranz coefficients.**
The three coefficients for physical line :math:`i` at the reference pressure
:math:`P_0` are read off as:

.. math::

  Y_{lm,i,s}(T) &= \frac{\mathrm{Im}\!\left[r_{i,s}(T)\right]}{P_0}, \\[4pt]
  G_{lm,i,s}(T) &= \frac{\mathrm{Re}\!\left[r_{i,s}(T)\right] - 1}{P_0^2}, \\[4pt]
  \Delta\nu_{lm,i,s}(T) &= \frac{\mathrm{Re}\!\left[\delta\nu_{i,s}(T)\right]}{P_0^2}.

The imaginary part of :math:`\delta\nu_{i,s}(T)` — which represents the
correction to the pressure-broadening half-width — is divided by
:math:`P_0^3` but is not retained as a separate Rosenkranz coefficient
(it is already captured by the diagonal of :math:`\mathbf{W}` at first
order in :math:`P`).

**Step 4 — Polynomial fit in temperature.**
Each of the three coefficients is fitted as a polynomial in temperature of
configurable degree :math:`d`:

.. math::

  Y_{lm,i,s}(T)         &\approx \sum_{n=0}^{d} a_n^{(Y)}\,T^n, \\
  G_{lm,i,s}(T)         &\approx \sum_{n=0}^{d} a_n^{(G)}\,T^n, \\
  \Delta\nu_{lm,i,s}(T) &\approx \sum_{n=0}^{d} a_n^{(DV)}\,T^n.

These polynomials are stored using the ``POLY`` temperature model (see
:ref:`lbl-line-shape-params`) for each broadening species separately and
are then evaluated at runtime using the ordinary VMR-weighted sum of the
:ref:`line shape parameter <lbl-line-shape-params>` framework, with the
coefficients attached to physical line :math:`i`.

.. note::

  Setting ``rosenkranz_fit_order = 1`` retains only :math:`Y_{lm}` and is
  appropriate for moderate pressures where the quadratic-in-pressure corrections
  are negligible.  Setting it to 2 also fits :math:`G_{lm}` and
  :math:`\Delta\nu_{lm}`, which is necessary at higher pressures or when
  second-order effects are important (e.g. near the Q-branch of O :sub:`2`
  at tens of GHz).

  The polynomial fits implicitly assume that the reference pressure :math:`P_0`
  is fixed; the resulting coefficients must be used at the same pressure
  normalisation.  This is handled automatically when the output of
  ``abs_bandsLineMixingAdaptation`` is fed back into the ordinary line-by-line
  calculation.

