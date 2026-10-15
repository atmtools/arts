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

where :math:`i` is the pseudo-index of the absorption line and :math:`\mathrm{lbl}` is the pseudo-index of the plain line-by-line absorption for the sake of :ref:`summing up absorption <eq-prop-mat-sumup>`.

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

Given a spherical coordinate observation system with zenith angle :math:`\theta_z` and azimuth angle :math:`\eta_a` and a
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

  z = \frac{\nu - \nu_0 - \Delta\nu_{lm} - \Delta\nu_{Z} - \Delta\nu_{P,0} + iG_{P,0}}{G_D},

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

Line-mixing using Energy-Corrected Sudden
********************************************

When the atmosphere is at sufficient pressure, collisions occur frequently enough that
absorption lines of a vibrational-rotational band can no longer be treated
independently.  Molecules undergoing collisions may exchange rotational angular
momentum, transferring population between rotational levels.  At intermediate
pressures this introduces off-diagonal couplings between lines in a spectral band,
leading to the phenomenon of *line mixing*.  This redistributes absorption across
the band and can strongly change its shape.

The Energy-Corrected Sudden (ECS) approximation provides a semi-empirical
framework for computing this mixing.  The "Sudden" part refers to the
Infinite-Order-Sudden (IOS) approximation, in which the collision time is assumed short
compared with the rotational period.  The energy correction introduces an
adiabatic factor accounting for a finite collision duration.  The additional
sum-rule rescaling used below is a separate step; it does not make the model exact.

.. _lbl-ecs-lineshape:

ECS Line Shape
==============

For a band of :math:`n` interacting absorption lines, the ECS complex absorption shape
for a single broadening species (see :ref:`lbl-ecs-multispecies` for the full
expression) is written in terms of the complex relaxation matrix
:math:`\mathbf{W}` of :ref:`lbl-ecs-relaxmat` as

.. math::

  \chi(\nu) \propto \mathrm{Im}\left[\mathbf{d}^T \left(\nu \mathbf{I} - \mathbf{W}\right)^{-1} \mathbf{p}\, \mathbf{d}\right].

The CO\ :sub:`2` and O\ :sub:`2` models use the diagonal weights
:math:`p_j=g_{u,j}\exp(-E_l^{(j)}/kT)/Q(T)` and signed amplitudes

.. math::

  d_j = \mathrm{sgn}(d_{r,j})\,\frac{c}{2}
        \sqrt{\frac{A_j}{2\pi\nu_{0,j}^3}}.

Here :math:`g_{u,j}` is the upper-state statistical weight, :math:`A_j` the
Einstein A coefficient, and :math:`d_{r,j}` the angular reduced dipole used by
the sum-rule correction.  The weights :math:`p_j` therefore should not be
interpreted as the total lower-state populations. A change of normalisation
requires a consistent transformation of the amplitudes, weights, and angular
kernels.

The matrix is diagonalised as :math:`\mathbf{W} = \mathbf{V} \tilde{\boldsymbol{\nu}} \mathbf{V}^{-1}`,
where :math:`\tilde{\boldsymbol{\nu}}` is the diagonal matrix of complex
*equivalent line* positions.  Each equivalent line :math:`k` has a complex
frequency :math:`\tilde{\nu}_k` (real part: position, imaginary part: pressure
broadening) and a complex *equivalent strength* :math:`\tilde{S}_k`.

Explicitly, the equivalent strength for line :math:`k` is

.. math::

  \tilde{S}_k = \left(\sum_j d_j V_{jk}\right) \left(\sum_j p_j d_j V^{-1}_{kj}\right).

Equivalently, the second factor is :math:`a_k`, where
:math:`\mathbf{V}\mathbf{a}=\mathbf{p}\mathbf{d}`. This representation requires
a complete eigenvector basis; nonnegative damping is necessary for physical
absorption profiles.

The ECS line shape function is

.. math::

  F_{ECS}(\nu) = \frac{1}{\sqrt{\pi}} \sum_k \tilde{S}_k \frac{w(z_k)}{G_{D,k}},

where :math:`w` is the Faddeeva function and

.. math::

  z_k = \frac{\tilde{\nu}_k - \nu}{G_{D,k}}, \qquad
  G_{D,k} = G_D^{fac} \cdot \mathrm{Re}\!\left[\tilde{\nu}_k\right],

with the Doppler scale factor

.. math::

  G_D^{fac} = \sqrt{\frac{2000 R T}{m c^2}},

where :math:`R` is the ideal gas constant in J mol\ :sup:`-1` K\ :sup:`-1`,
:math:`m` the molar mass in g mol\ :sup:`-1`, :math:`c` the speed of light, and
:math:`T` the temperature.  Here :math:`G_{D,k}` is the Gaussian :math:`1/e`
half-width, as in :ref:`lbl-line-shape`, rather than its half-width at half-maximum.

Assigning each equivalent line the Doppler width of
:math:`\mathrm{Re}[\tilde{\nu}_k]` is an approximation.  It recovers isolated
Voigt profiles when coupling vanishes, but does not generally equal the velocity
average of the coupled resolvent when the physical lines have different Doppler
widths.  Diagonalising a collision matrix and averaging over velocities need not
commute.  Bands spanning a wide frequency range require particular care.

The contribution to the :ref:`propagation matrix <prop-mat>` from the entire band is
then

.. math::

  K_{A, ecs} = N \nu \left(1 - \exp\!\left(-\frac{h\nu}{kT}\right)\right) \mathrm{Re}\!\left[F_{ECS}(\nu)\right],

where :math:`N` is the total number density of the absorbing species.

ECS Jacobians
-------------

Derivatives of the absorption spectrum follow from the dependence of the
collision matrix, optical populations, dipole amplitudes, and Doppler widths on
the perturbed quantity. Temperature, pressure, composition, isotopic abundance,
line frequencies, lower-state energies, Einstein coefficients, and collisional
width and shift parameters can all contribute. A perturbation of a spectroscopic
lower-state energy changes the optical population; the collision energies
remain determined by the chosen rotational-state model unless that model is
also perturbed.

Derivatives of individual equivalent lines require resolved, distinct modes
and a well-conditioned eigenvector basis. Degenerate modes require a coupled
subspace treatment. Perturbations that leave the centered spectral operator
fixed change the populations, amplitudes, or common frequency offset without
changing its eigenvectors.

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

The off-diagonal elements are purely imaginary.  Write
:math:`R_{ij}=\mathrm{Im}\,W_{ij}` for the real collisional coupling
coefficients, where :math:`j` is the state the collision transfers population
from and :math:`i` the state it transfers to. Their construction from the ECS
basis rates is described below.

.. _lbl-ecs-rates:

ECS Basis Rates
===============

The CO\ :sub:`2` and O\ :sub:`2` ECS models introduce two species-dependent
functions of the integer angular momentum transfer channel :math:`L`:

**Basic rate** :math:`Q(L)`:
  This encodes the intrinsic probability of a collision transferring :math:`L` units of
  angular momentum.  Its temperature dependence is parameterised as

  .. math::

    Q(L, T) = s(T) \cdot \frac{e^{-\beta(T)\, E_L / kT}}{[L(L+1)]^{\lambda(T)}},

  where :math:`E_L` is the rotational energy of level :math:`L`,
  and :math:`s(T)`, :math:`\beta(T)`, and :math:`\lambda(T)` are
  temperature-dependent model parameters specified per broadening species
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

The CO\ :sub:`2` and O\ :sub:`2` off-diagonal elements follow the formal IOS
structure: both models are written as a sum over even angular momentum transfer
channels :math:`L`, weighted by the ratio
:math:`Q(L)/\Omega(L)` and by Wigner 3-j and 6-j coupling coefficients.
After the IOS computation a sum-rule rescaling is applied
(see :ref:`lbl-ecs-sumrule`).

For these two models, the reverse coupling is assigned using the following
Boltzmann relation:

.. math::

  R_{ij} = R_{ji} \exp\!\left(\frac{E_j - E_i}{kT}\right),

where :math:`E_i` is the energy of the lower rotational state of line :math:`i`
(using the same :math:`E_l` convention as in :ref:`lbl-lte`).

The Hartmann CO\ :sub:`2` band model and Makarov O\ :sub:`2` microwave band
model have different angular couplings and molecular energy models, described
below.

Linear Molecules — Hartmann (CO\ :sub:`2`)
------------------------------------------

For linear molecules (e.g. CO\ :sub:`2`) the angular coupling for line :math:`j`
(upper/lower rotational quantum numbers :math:`J'_i, J'_f`, vibrational angular
momentum :math:`l`) to line :math:`i` (:math:`J_i, J_f`) is
:cite:p:`NIRO2004483`

.. math::

  R_{ji} =
    \Omega(J_i)\, (2J'_i+1)\sqrt{(2J_f+1)(2J'_f+1)}
    \sum_L (2L+1)
    \begin{pmatrix} J_i & J'_i & L \\ l_i & -l_i & 0 \end{pmatrix}
    \begin{pmatrix} J_f & J'_f & L \\ l_f & -l_f & 0 \end{pmatrix}
    \begin{Bmatrix} J_i & J_f & 1 \\ J'_f & J'_i & L \end{Bmatrix}
    \frac{Q(L)}{\Omega(L)},

where the sum runs over allowed even :math:`L\geq 2` satisfying the angular
momentum triangle bounds,
:math:`(\,\cdots)` denotes a Wigner 3-j symbol,
and :math:`\{\,\cdots\}` a Wigner 6-j symbol.

The upper and lower angular labels in this expression are interchanged when
:math:`l_i>l_f`.  The angular reduced dipole used for the sign and
sum-rule correction is, in upper/lower state notation,

.. math::

  d_r(J_u,J_l,l_u,l_l) = (-1)^{J_u+l_u+1}\sqrt{2J_u+1}\;
    \begin{pmatrix} J_u & 1 & J_l \\ l_u & l_l-l_u & -l_l \end{pmatrix}.

The rotational energy entering :math:`Q` and :math:`\Omega` is the rigid-rotor
expression :math:`E_J = hcB_0 J(J+1)`, with :math:`B_0=0.39021\,\mathrm{cm}^{-1}`
for CO\ :sub:`2`-626.  This is currently the only isotopologue with a rotational
energy model considered here. Other isotopologues require their own consistent
rotational energies as well as collision parameters.

Linear Molecules with Electron Spin — Makarov (O\ :sub:`2`)
-----------------------------------------------------------

Molecular oxygen (O\ :sub:`2`) has an unpaired electron spin :math:`S = 1`, so each
rotational quantum number :math:`N` gives rise to a triplet :math:`J = N-1, N, N+1`.
The off-diagonal coupling between lines :math:`(i: N_l J_l \to N_u J_u)` and
:math:`(j: N'_l J'_l \to N'_u J'_u)` is (using :math:`l`/:math:`u` for the
lower/upper state of each transition, matching the convention of :ref:`lbl-lte`)
:cite:p:`Makarov2020`

.. math::

  R_{ij} =&
    (-1)^{J'_u + J_u + 1}\,
    [N_l][N_u][N'_u][N'_l][J_u][J'_u][J_l][J'_l] \Omega(N_u) \\ &
    \begin{array}{llll}
      \sum_L (2L+1) &
      \begin{pmatrix} N'_l & N_l & L \\ 0 & 0 & 0 \end{pmatrix} &
      \begin{pmatrix} N'_u & N_u & L \\ 0 & 0 & 0 \end{pmatrix} \\ &
      \begin{Bmatrix} L & J_l & J'_l \\ S & N'_l & N_l \end{Bmatrix} &
      \begin{Bmatrix} L & J_u & J'_u \\ S & N'_u & N_u \end{Bmatrix} &
      \begin{Bmatrix} L & J_u & J'_u \\ 1 & J'_l & J_l \end{Bmatrix}
      \frac{Q(L)}{\Omega(L)},
    \end{array}

where :math:`[X] \equiv \sqrt{2X+1}`, :math:`(\cdots)` denotes a Wigner 3-j symbol,
and :math:`\{\cdots\}` a Wigner 6-j symbol.

The reduced dipole is

.. math::

  d_r(J_u, J_l, N) = (-1)^{J_l + N}
    \sqrt{6(2J_l+1)(2J_u+1)}
    \begin{Bmatrix} 1 & 1 & 1 \\ J_l & J_u & N \end{Bmatrix}.

The kernel is specific to O\ :sub:`2`-66.  A single set of molecular constants
:cite:p:`tretyakov05:_60-ghz_jms` supplies two energy functions:
a reference-rotor energy depending on :math:`N`, and a resolved spin-triplet
energy depending on both :math:`N` and :math:`J`.  Both use the :math:`N=1, J=0` ground
state as their energy zero.

The reference rotor contains the rotational and centrifugal-distortion terms.
Its energies enter :math:`Q(L)` and the :math:`N\leftrightarrow N-2` spacings
in :math:`\Omega`.  This retains the ECS approximation of a spinless reference
rotor; the angular factors recouple it to the spin-triplet states.  The resolved
level function also includes the spin--rotation and spin--spin terms.  Evaluated
at each line's actual lower :math:`(N_l,J_l)`, it supplies the Boltzmann factors
relating opposite matrix elements and the sum-rule correction.

The resolved energies retain approximate spin-triplet expressions.  The
:math:`N=1,J=0` special case has the spin correction :math:`-2\lambda-\gamma`,
giving a :math:`J=1\leftarrow0` splitting of 118.750334 GHz, compared with the
measured 118.750340 GHz.  Other tested :math:`N=1,3,5` branches differ from the
measured frequencies by up to about 24 MHz.  These energies therefore do not
replace catalogue line frequencies.  The common ground-state reference also
enters the absolute energy in :math:`Q`, where a reference shift does not cancel.

.. _lbl-ecs-nh3:

Symmetric Tops with Inversion — Hadded (NH\ :sub:`3`)
-----------------------------------------------------

The angular couplings and energy corrections of :cite:t:`Hadded2002`
(Eqs. 9--19) describe one parallel band: :math:`K_u=K_l`,
:math:`|J_u-J_l|\leq1`, and opposite upper/lower inversion symmetry
(:math:`a\leftarrow s` or :math:`s\leftarrow a`). This includes :math:`\nu_2`;
perpendicular bands are outside its domain.

The signed dipole and population must use the paper's common normalization:

.. math::

  \mu_i = (-1)^{J_{u,i}+K_i}\sqrt{2J_{u,i}+1}
    \begin{pmatrix}J_{u,i}&1&J_{l,i}\\K_i&0&-K_i\end{pmatrix},
  \qquad
  \rho_i = \frac{g_i(2J_{l,i}+1)}{Z(T)}
    \exp\!\left(-\frac{E_i}{kT}\right),

where :math:`Z(T)` is the partition function and :math:`g_i` is 4 for ortho and
2 for para states. Collisions do not couple the ortho and para nuclear-spin
blocks.

For the IOS angular kernel (Eq. 10), use unprimed quantum numbers for the
source line :math:`i` and primed quantum numbers for the destination line
:math:`j`. Subscripts :math:`l,u` denote each line's lower and upper state
(the paper uses :math:`i,f` for these state labels). With
:math:`R_{ij}=\langle j|\hat R|i\rangle`, the full expression is

.. math::

  R_{ij}^{\mathrm{IOS}} ={}&
    -N_l N_u N'_l N'_u (2J'_l+1)
    \sqrt{(2J'_u+1)(2J_u+1)} \\
    &\times \sum_L (-1)^{J'_u+J_u+K'_u+K'_l+1+L}
    \begin{Bmatrix}J_l&J_u&1\\J'_u&J'_l&L\end{Bmatrix}
    \mathcal{A}_L,

where the four products of Wigner 3-j symbols and dynamical factors are

.. math::

  \mathcal{A}_L ={}&
    A_l A_u
    \begin{pmatrix}J'_l&L&J_l\\K'_l&M_l^-&-K_l\end{pmatrix}
    \begin{pmatrix}J'_u&L&J_u\\K'_u&M_u^-&-K_u\end{pmatrix}
    Q(L,M_l^-,M_u^-) \\
    &+ B_l A_u
    \begin{pmatrix}J'_l&L&J_l\\-K'_l&M_l^+&-K_l\end{pmatrix}
    \begin{pmatrix}J'_u&L&J_u\\K'_u&M_u^-&-K_u\end{pmatrix}
    Q(L,M_l^+,M_u^-) \\
    &+ A_l B_u
    \begin{pmatrix}J'_l&L&J_l\\K'_l&M_l^-&-K_l\end{pmatrix}
    \begin{pmatrix}J'_u&L&J_u\\-K'_u&M_u^+&-K_u\end{pmatrix}
    Q(L,M_l^-,M_u^+) \\
    &+ B_l B_u
    \begin{pmatrix}J'_l&L&J_l\\-K'_l&M_l^+&-K_l\end{pmatrix}
    \begin{pmatrix}J'_u&L&J_u\\-K'_u&M_u^+&-K_u\end{pmatrix}
    Q(L,M_l^+,M_u^+).

For :math:`x\in\{l,u\}`, the projections and inversion factors are

.. math::

  M_x^\pm &= K_x\pm K'_x, \qquad
  P_x = (-1)^{J'_x+K'_x+J_x+K_x+L}, \\
  A_x &= 1+\epsilon_x\epsilon'_x P_x, \qquad
  B_x = \epsilon'_x+\epsilon_x P_x.

The normalization and inversion sign for each state, and analogously for
primed states, follow Eq. 5:

.. math::

  N_x =
    \begin{cases}1,&K_x=0,\\1/\sqrt{2},&K_x>0,\end{cases}
  \qquad
  \epsilon_x =
    \begin{cases}
      0,&K_x=0,\\
      (-1)^{J_x},&K_x>0\text{, antisymmetric }(a),\\
      (-1)^{J_x+1},&K_x>0\text{, symmetric }(s).
    \end{cases}

The braces denote a Wigner 6-j symbol. The integer sum spans
:math:`\max(|J_l-J'_l|,|J_u-J'_u|)\leq L\leq
\min(J_l+J'_l,J_u+J'_u)`, subject to the 3-j projection conditions and the
nonzero dynamical factors. Allowed odd and even :math:`L` are included;
there is no additional :math:`(2L+1)` multiplier. The projections are signed
multiples of 3 for NH\ :sub:`3`. Distinct :math:`Q(L,M_l,M_u)` factors cannot
in general be identified by taking absolute values of the projections.

For downward transfer from line :math:`i` to line :math:`j`, selected by
:math:`E_j\leq E_i`, the ECS correction replaces :math:`Q` by
:math:`Q'(L,M_l,M_u)=Q(L,M_l,M_u)\,\Omega(L,M_l)` inside the angular sum
and divides the result by
:math:`\Omega(J_{l,i},K_{l,i})`. Here the paper's convention is

.. math::

  \Omega = \left[1+\frac{(\tau\,\Delta E/\hbar)^2}{24}\right]^2,
  \qquad \tau=\ell_c/\bar v,
  \qquad \rho_i R_{ij}=\rho_j R_{ji}.

Detailed balance therefore includes the lower-state rotational degeneracy in
:math:`\rho`. This :math:`\Omega\geq1` is reciprocal to the
CO\ :sub:`2`/O\ :sub:`2` convention above. Setting all factors to one gives IOS
with detailed balance enforced. The energies :math:`E_i` and adiabatic gaps
must come from one consistent molecular energy model. The gap connects a
rotational level to an appropriate lower-energy level; its selection is part
of the physical collision model.

The diagonal widths are independent spectroscopic parameters in these studies.
An alternative estimate follows from the optical sum rule
:math:`\sum_j\mu_jR_{ij}=0`. Applied to a finite set of transitions, this
estimate depends on which lines are retained. It differs from the
CO\ :sub:`2`/O\ :sub:`2` rescaling described below, which adjusts off-diagonal
couplings to supplied widths.

The He study and its H\ :sub:`2`/Ar extension :cite:p:`Hadded2004` compare with
room-temperature measurements. They do not establish a general calibrated
H\ :sub:`2`/He temperature law for planetary atmospheres. Extension to other
conditions requires appropriate collision factors, widths, collision duration,
and their temperature dependence.

Interface guidance and a plotted :math:`\nu_2` example are in :doc:`user.lbl`;
implementation conventions are in :doc:`dev.lbl`.

.. _lbl-ecs-sumrule:

CO\ :sub:`2`/O\ :sub:`2` Sum-rule Correction
============================================

For CO\ :sub:`2` and O\ :sub:`2`, the approximate angular couplings and supplied
diagonal widths need not satisfy the optical sum rule on a finite set of lines.
The rescaling seeks to impose

.. math::

  \sum_j d_{r,j}\, R_{ij} = 0.

Lines are ordered by decreasing :math:`\nu_{0,i}p_i d_i^2`.  For each row
:math:`i`, define the contribution from entries still available to rescale and
the contribution from the diagonal and entries already fixed:

.. math::

  s_\downarrow = \sum_{j > i} d_{r,j}\, R_{ij}, \qquad
  s_\uparrow   = \sum_{j \leq i} d_{r,j}\, R_{ij}.

When :math:`s_\downarrow\ne 0`, the remaining entries are rescaled by
:math:`-s_\uparrow/s_\downarrow`, and their reverse couplings are updated using
the CO\ :sub:`2`/O\ :sub:`2` Boltzmann relation above.

Both the raw reverse couplings and this correction use each line's original
lower-state rotational energy. Hartmann uses the lower :math:`J`, even when
its angular labels are interchanged; Makarov uses the resolved
:math:`(N_l,J_l)` energy. The reference-rotor energies :math:`E_L` and
:math:`E_{L-2}` entering the collision basis must belong to the same molecular
energy model. These collision energies remain distinct from the spectroscopic
energies determining the optical populations.

This sequential prescription cannot generally enforce every row: the final
row has no remaining entries to adjust, and zero or nearly cancelling sums
can prevent a stable correction.  Truncating a band or a connected symmetry
block can also leave a nonzero residual. Sum-rule closure and nonnegative
absorption must be checked for the band and conditions of interest. The
rescaling alone is not a proof of either property.

.. _lbl-ecs-multispecies:

Multiple Broadening Species
===========================

For a gas mixture, the collision contributions of the individual partners are
volume-mixing-ratio weighted and summed into a single effective
:math:`\mathbf{W}` before diagonalisation.  With
:math:`\mathbf{F}=\operatorname{diag}(\nu_{0,i})`, this is

.. math::

  \mathbf{W}_{eff} = \mathbf{F}
    + \sum_s x_s\,\bigl(\mathbf{W}^{(s)}-\mathbf{F}\bigr),

and a single diagonalisation is performed.

For fractions summing to one, the expression reduces to
:math:`\mathbf{W}_{eff}=\sum_s x_s\mathbf{W}^{(s)}`. With no collisions,
only the unperturbed line frequencies remain.

Each pure partner may also be diagonalised separately, for example to obtain
its equivalent lines for a pressure-expansion fit. Summing the resulting
pure-partner spectra is not generally equivalent to diagonalising the combined
collision matrix. Combining the matrices preserves the distinct thermal
dependence and collision dynamics of each partner, including their
contributions to higher orders of line mixing.

A catalogue may provide only a mean-air diagonal width, rather than separate
O\ :sub:`2` and N\ :sub:`2` widths. One approximation is then to combine that
width with a bath model obtained by averaging the collision parameters of its
constituents. This coefficient average is an approximate bath model: the ECS
basis rates and adiabatic factors depend nonlinearly on those parameters, so
the resulting matrix generally differs from a weighted sum of separately
constructed partner matrices. A mean-air width alone does not determine the
individual partner widths.

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

For weak coupling, :math:`\mathbf{W}` can be expanded about its diagonal
part.  Write
:math:`\mathbf{W}=\mathbf{D}+\mathbf{U}`, with :math:`\mathbf{U}` containing
only off-diagonal entries, and define
:math:`\mathbf{G}_0=(\nu\mathbf{I}-\mathbf{D})^{-1}`.  Where the Neumann
series converges,

.. math::

  (\nu\mathbf{I}-\mathbf{W})^{-1} =
    \mathbf{G}_0+\mathbf{G}_0\mathbf{U}\mathbf{G}_0
    +\mathbf{G}_0\mathbf{U}\mathbf{G}_0\mathbf{U}\mathbf{G}_0+\cdots.

The absorption follows by contracting each term with :math:`\mathbf{d}^T`
and :math:`\mathbf{p}\mathbf{d}`, using the amplitude and population convention
in :ref:`lbl-ecs-lineshape`.  The first-order term contains products
:math:`d_i d_j p_j U_{ij}`.  Second-order terms contain products
:math:`d_i d_k p_k U_{ij}U_{jk}` summed over intermediate lines :math:`j`.
Consequently they cannot in general be represented by a squared single
coupling weighted only by the isolated line strengths.

For example, if
:math:`\mathbf{W}=\mathrm{diag}(\nu_{0,i})+P\mathbf{B}` and the unperturbed
frequencies are distinct, the equivalent frequency has the expansion

.. math::

  \tilde{\nu}_i = \nu_{0,i}+P B_{ii}
    +P^2\sum_{j\ne i}\frac{B_{ij}B_{ji}}{\nu_{0,i}-\nu_{0,j}}
    +O(P^3).

For purely imaginary off-diagonal collision terms this second-order product
is negative times the product of the two real coupling coefficients.  The
corresponding strength correction also depends on the left and right amplitude
projections.  Degenerate or nearly degenerate transitions require a coupled
subspace treatment rather than these nondegenerate formulas.

.. _lbl-ecs-rosenkranz-fitting:

Fitting from Equivalent Lines
------------------------------

Rosenkranz coefficients can be extracted numerically from equivalent lines at a
chosen reference pressure.  The equivalent-line calculation includes all orders
of coupling, but retaining only linear and quadratic pressure terms in the
adapted model is still an approximation.  At finite reference pressure the
extracted coefficients can include higher-order contributions; their accuracy
must be checked over the intended pressure and temperature range.

Given the complex equivalent lines :math:`(\tilde{S}_{k,s}, \tilde{\nu}_{k,s})`
(indexed by :math:`k` in eigenvalue-decomposition order, which carries no physical
meaning) computed by ECS for broadening species :math:`s` at pressure :math:`P_0`
and a grid of temperatures :math:`T_1, \ldots, T_M`, the Rosenkranz coefficients are
obtained by the following procedure.
Here, :math:`i` will denote the index of the physical LBL lines (the rows/columns of
:math:`\mathbf{W}`) and :math:`k` the index of the equivalent lines.

**Step 1 — Sort and match.**
The :math:`n` equivalent lines are sorted by :math:`\mathrm{Re}[\tilde{\nu}_{k,s}]`
and the :math:`n` physical LBL lines are sorted by :math:`\nu_{0,i}`.  Equivalent line
at sorted position :math:`n` is then identified with the physical line at sorted
position :math:`n`, giving a bijection :math:`k \leftrightarrow i` between the two
index sets.  This is a heuristic matching: because eigenvalue decomposition does not
guarantee any particular ordering of eigenvalues, sorting defines the chosen
correspondence.  The identification is reliable as long as the second-order frequency
shifts and pressure shifts remain small compared with the separations between adjacent
line centres; it can fail if either effect is comparable in magnitude to the line spacing.

**Step 2 — Form normalised differences.**
For each matched pair :math:`(k, i)`, define the unperturbed ECS strength of
physical line :math:`i` as :math:`s_i(T)=p_i(T)d_i^2`.  This excludes number
density and the frequency-dependent stimulated-emission factor.  The
unperturbed complex frequency of physical
line :math:`i` is

.. math::

  \nu_i^{LBL}(T) = \nu_{0,i} + \Delta\nu_{P,0,i}(T,P_0) + i\,G_{P,0,i}(T,P_0).

The residual strength ratio and frequency residual are then formed:

.. math::

  r_{i,s}(T) &= \frac{\tilde{S}_{k,s}(T)}{s_i(T)}, \\
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
:math:`P_0^3` but is not retained as a separate Rosenkranz coefficient.
This higher-order correction is not supplied by the first-order diagonal width.

**Step 4 — Polynomial fit in temperature.**
Each of the three coefficients is fitted as a polynomial in temperature of
configurable degree :math:`d`:

.. math::

  Y_{lm,i,s}(T)         &\approx \sum_{n=0}^{d} a_n^{(Y)}\,T^n, \\
  G_{lm,i,s}(T)         &\approx \sum_{n=0}^{d} a_n^{(G)}\,T^n, \\
  \Delta\nu_{lm,i,s}(T) &\approx \sum_{n=0}^{d} a_n^{(DV)}\,T^n.

These polynomials describe each broadening species separately. Their
contributions to a physical line :math:`i` are combined using the
volume-mixing-ratio weights of the :ref:`line shape parameter
<lbl-line-shape-params>` model.

.. note::

  A first-order expansion retains only :math:`Y_{lm}` and is appropriate for
  moderate pressures where the quadratic-in-pressure corrections are negligible.
  A second-order expansion also retains :math:`G_{lm}` and
  :math:`\Delta\nu_{lm}` within the range where a perturbation expansion remains
  accurate; increasing the order does not make it valid for arbitrarily strong
  mixing.

  Adaptation fits each collision partner separately.  For a mixture, second-order
  perturbation theory contains products :math:`x_sx_t\mathbf{R}^{(s)}
  \mathbf{R}^{(t)}` between different partners.  Averaging independently fitted
  per-partner coefficients cannot in general reproduce these cross terms.  A
  second-order adapted mixture therefore need not agree with the full ECS
  calculation that combines the matrices before diagonalisation.

  The indicated powers of pressure are restored when evaluating the expansion.
  Coefficients inferred at one finite :math:`P_0` need not reproduce the full
  ECS calculation at another pressure, even with an accurate temperature fit.
