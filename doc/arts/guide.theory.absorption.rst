Absorption
##########

Absorption is the physical process that reduces 
the spectral radiance of a beam of light as it passes through a medium.
It can be described be Beer's law:

.. math::

  \vec{I}(r) = \vec{I}(0) \exp(-\mathbf{K} r),

where :math:`\vec{I}(r)` is the spectral radiance as a :ref:`Stokes vector <stok-vec>` of the light at some distance
:math:`r`, and :math:`\mathbf{K}` is the :ref:`propagation matrix <prop-mat>` of the medium.

The unit of spectral radiance is W sr :math:`^{-1}` m :math:`^{-2}` Hz :math:`^{-1}`.  The unit of the propagation matrix is m :math:`^{-1}`.

.. _stok-vec:

The Stokes Vector
*****************

In ARTS, the spectral radiance is represented by the Stokes vector, which is a 4-long column vector.

.. math::

  \vec{I} = \left[ \begin{array}{c} I_I \\ I_Q \\ I_U \\ I_V \end{array} \right],

where 
:math:`I_I` is the total spectral radiance,
:math:`I_Q` is the difference in spectral radiance between horizontal and vertical linear polarizations,
:math:`I_U` is the difference in spectral radiance between plus 45 and minus 45 linear polarizations, and
:math:`I_V` is the difference in spectral radiance between right and left circular polarizations.
The unit of all of these is W sr :math:`^{-1}` m :math:`^{-2}` Hz :math:`^{-1}`.

This is the same definition as, e.g., :cite:t:`Mishchenko:02`.

.. tip::

  Because of the nature of circular polarization, it is not uncommon to see different :math:`I_V` definitions
  in different context.  Be sure to check the definition of :math:`I_V` in the context you are working in.
  It has caused confusion in the past, for both users and developers of ARTS.

The Stokes Vector Unit Conversion
=================================

All ARTS simulations expect the incoming spectral radiance to be in W sr :math:`^{-1}` m :math:`^{-2}` Hz :math:`^{-1}`.

However, several ways to convert the Stokes vector to other units are available.  These are listed below in no particular order.

Rayleigh Jeans Brightness Temperature
--------------------------------------

As the spectral temperature of the Rayleigh Jeans interpretation of black body spectral radiance.   The conversion is given by

.. math::
  \vec{I}_{RJBT} = \frac{c^2}{2k\nu^2} \left[ \begin{array}{c} I_I \\ I_Q \\ I_U \\ I_V \end{array} \right],

where :math:`k` is the Boltzmann constant, :math:`c` is the speed of light, and :math:`\nu` is the frequency.

The unit of all of these is K.

Planck Brightness Temperature
-----------------------------

As the spectral temperature of the black body spectral radiance.   The conversion is given by

.. math::

  \vec{I}_{Tb} = \frac{h \nu}{k} \left[ \begin{array}{lcr} 
  \log\left(1 + \frac{2h\nu^3}{c^2I_I}\right)^{-1} &&\\
  \log\left(1 + \frac{4h\nu^3}{c^2(I_I + I_Q)}\right)^{-1} &-& \log\left(1 + \frac{4h\nu^3}{c^2(I_I - I_Q)}\right)^{-1} \\
  \log\left(1 + \frac{4h\nu^3}{c^2(I_I + I_U)}\right)^{-1} &-& \log\left(1 + \frac{4h\nu^3}{c^2(I_I - I_U)}\right)^{-1} \\
  \log\left(1 + \frac{4h\nu^3}{c^2(I_I + I_V)}\right)^{-1} &-& \log\left(1 + \frac{4h\nu^3}{c^2(I_I - I_V)}\right)^{-1}
  \end{array} \right],

where :math:`k` is the Boltzmann constant, :math:`h` is the Planck constant, :math:`c` is the speed of light, and :math:`\nu` is the frequency.

The unit of all of these is K.

The Spectral Radiance per Wavelength
------------------------------------

As the spectral temperature of the Rayleigh Jeans interpretation of black body spectral radiance.   The conversion is given by

.. math::

  \vec{I}_{\lambda} = \frac{\nu^2}{c} \left[ \begin{array}{c} I_I \\ I_Q \\ I_U \\ I_V \end{array} \right],

where :math:`c` is the speed of light, and :math:`\nu` is the frequency.

The unit of all of these is W sr :math:`^{-1}` m :math:`^{-2}` m :math:`^{-1}`.

The Spectral Radiance per Kayser
--------------------------------

As the spectral temperature of the Rayleigh Jeans interpretation of black body spectral radiance.   The conversion is given by

.. math::

  \vec{I}_{cm} = \frac{1}{c} \left[ \begin{array}{c} I_I \\ I_Q \\ I_U \\ I_V \end{array} \right],

where :math:`c` is the speed of light.

The unit of all of these is W sr :math:`^{-1}` m :math:`^{-2}` m.

.. _prop-mat:

Propagation Matrix
******************

The propagation matrix conceptually describes how :ref:`the Stokes vector <stok-vec>`
propagates through a system. The propagation matrix is a square matrix
with strict symmetries, and it has the form

.. math::

   \mathbf{K} = \left[ \begin{array}{rrrr}
        K_A & K_B & K_C & K_D \\
        K_B & K_A & K_U & K_V \\
        K_C &-K_U & K_A & K_W \\
        K_D &-K_V &-K_W & K_A
    \end{array} \right],

where
:math:`K_A` describes the total power reduction,
:math:`K_B` describes the difference in power reduction between horizontal and vertical linear polarizations,
:math:`K_C` describes the difference in power reduction between plus 45 and minus 45 linear polarizations,
:math:`K_D` describes the difference in power reduction between right and left circular polarizations,
:math:`K_U` describes the phase delay between right and left circular polarizations,
:math:`K_V` describes the phase delay between plus 45 and minus 45 linear polarizations, and
:math:`K_W` describes the phase delay between horizontal and vertical linear polarizations.
The unit of all of these is m :math:`^{-1}`.

In practice, the propagation matrix is a sum of multiple physical processes:

.. _eq-prop-mat-sumup:

.. math::

  \mathbf{K} = \sum_i \mathbf{K}_i

where :math:`i` is the pseudo-index of the physical process.  A physical process here
refers to line-by-line absorption, collision-induced absorption, Faraday rotation, scattering, etc.
Think of this sum as if each contribtion the the propagation matrix is completely independent of the others.

The way to compute different :math:`\mathbf{K}_i` are described in these sections:

.. toctree::
   :maxdepth: 2
   
   guide.theory.absorption.cia
   guide.theory.absorption.faraday
   guide.theory.absorption.lbl
   guide.theory.absorption.lookup
   guide.theory.absorption.predef
   guide.theory.absorption.xsec
  