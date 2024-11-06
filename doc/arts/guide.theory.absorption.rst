.. _Sec Absorption:

Absorption
##########

Absorption is the physical process that reduces 
the `<Spectral Radiance_>`_ of a beam of light as it passes through a medium.
It can be described be Beer's law:

.. math::
  \vec{I}(r) = \vec{I}(0) \exp(-\mathbf{K} z),

where :math:`\vec{I}(r)` is the `<Spectral Radiance_>`_ of the light at some distance
:math:`r`, and :math:`\mathbf{K}` is the `<Propagation Matrix_>`_ of the medium.

Propagation Matrix
******************

The propagation matrix conceptually describes how `<Spectral Radiance_>`_
propagates through a system. The propagation matrix is a square matrix
with strict symmetries, and it has the form

.. math::

   \mathbf{K} = \left[ \begin{array}{rrrr}
        A & B & C & D \\
        B & A & U & V \\
        C &-U & A & W \\
        D &-V &-W & A
    \end{array} \right],

where
:math:`A` describes the total power reduction,
:math:`B` describes the difference in power reduction between horizontal and vertical linear polarizations,
:math:`C` describes the difference in power reduction between plus 45 and minus 45 linear polarizations,
:math:`D` describes the difference in power reduction between right and left circular polarizations,
:math:`U` describes the phase delay between right and left circular polarizations,
:math:`V` describes the phase delay between plus 45 and minus 45 linear polarizations, and
:math:`W` describes the phase delay between horizontal and vertical linear polarizations.
The unit of all of these is m :math:`^{-1}`.

Spectral Radiance
*****************

The unit of spectral radiance is W sr :math:`^{-1}` m :math:`^{-2}` Hz :math:`^{-1}`.

Line-by-line Absorption
***********************

This section describes the physical process of absorption lines of 
different molecules absorbing and emitting `<Spectral Radiance_>`_
in the atmosphere.

These are the types of line-by-line absorption considered here:

#. `Plain line-by-line absorption <Plain Line-by-line Absorption_>`_,
   where each absorption line is considered separately.
#. Zeeman splitting, where the otherwise single absorption line is split
   into multiple lines due to the presence of a magnetic field.
#. Line-mixing by means of error-corrected sudden approximation,
   where the absorption lines of similar energies of a molecule
   are mixed together.

Plain Line-by-line Absorption
=============================

The absorption in plain line-by-line absorption is simply the sum of all 
absorption by each absorption line.  The absorption of a single absorption line
is described by the following equations:

.. math::

  \alpha = S(\cdots) F(\cdots),

where
:math:`\alpha` is the absorption coefficient,
:math:`S` is the `<Line Strength_>`_ operator, and
:math:`F` is the `<Line Shape_>`_ operator.

Line Shape
==========

Line shapes should distribute absorption as a function of frequency.
By convention, the line shape is normalized to have an integral of 1.

The following line shapes are available in ARTS 3.

Voigt Line Shape
----------------

.. math::

  F = \frac{1 + G_{lm} - iY_{lm}}{\sqrt{\pi}G_D} w(z),

where

.. math::

  z = \frac{\nu - \nu_0 - \Delta\nu_{lm} - \Delta\nu_Z - \Delta\nu_{P,0} + iG_{P,0}}{G_D},

:math:`\nu` is the frequency,
:math:`\nu_0` is the line center frequency,
:math:`\Delta\nu_{lm}` is the line mixing shift,
:math:`\Delta\nu_Z` is the Zeeman splitting,
:math:`\Delta\nu_{P,0}` is the pressure shift, 
:math:`G_{P,0}` is the pressure broadening - half width half maximum in the Lorentz profile, 
:math:`G_{lm}` is the strength modifying line mixing parameter,
:math:`Y_{lm}` is the phase-introducing line mixing parameter, and
:math:`G_D` is the scaled Doppler broadening half-width half-maximum.
:math:`\nu` is just a sampling frequency, it can be anything positive.
:math:`\nu_0` is provided by the absorption line catalog.
The :math:`\Delta\nu_{lm}`, :math:`\Delta\nu_{P,0}`, and :math:`G_{P,0}` are `<Line Shape Parameters_>`_.  :math:`\Delta\nu_Z`
is the Zeeman splitting, which depends on molecule and magnetic field strength as is described
in `Zeeman effect <Zeeman Effect_>`_.  The scaled Doppler broadening half width half maximum is given by

.. math::

  G_D = \sqrt{\frac{2000 R T}{mc^2}} \nu_0,

where
:math:`R` is the ideal gas constant in Joules per mole per Kelvin,
:math:`T` is the temperature in Kelvin,
:math:`m` is the mass of the molecule in grams per mole, and
:math:`c` is the speed of light in meters per second.
The factor 2000 is to convert to SI units.

Line Strength
=============

Local thermodynamic equilibrium
-------------------------------

For local thermodynamic equilibrium (LTE), the line strength is given by

.. math::

  S_{LTE} = \rho \frac{c^2\nu}{8\pi} \left[1 - \exp\left(-\frac{h\nu}{kT}\right)\right]
  \frac{g_u\exp\left(-\frac{E_l}{kT}\right)}{Q(T)} \frac{A_{lu}}{\nu_0^3}

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
this should be relatively close in cases where we actualy care about non-LTE
(which is low density, low collision atmospheres).

