# -*- coding: utf-8 -*-
"""Collection of physical constants and conversion factors.

Physical constants
==================

=============================== ==============================================
``g``                           Earth standard gravity in :math:`\sf ms^{-1}`
``h``                           Planck constant in :math:`\sf Js`
``k``                           Boltzmann constant in :math:`\sf JK^{-1}`
``c``                           Speed of light in :math:`\sf ms^{-1}`
``N_A``                         Avogadro constant in :math:`\sf mol^{-1}`
``K``, ``zero_celsius``         Kelvin at 0 Celsius
``triple_point_water``          Triple point temperature of water :math:`\sf K`
``R``                           Universal gas constant in
                                :math:`\sf J mol^{-1}K{^-1}`
``molar_mass_dry_air``          Molar mass for dry air in
                                :math:`\sf kg\,mol^{-1}`
``molar_mass_water``            Molar mass for water vapor in
                                :math:`\sf kg\,mol^{-1}`
``gas_constant_dry_air``        Gas constant for dry air in
                                :math:`\sf J K^{-1} kg^{-1}`
``gas_constant_water_vapor``    Gas constant for water vapor in
                                :math:`\sf J K^{-1} kg^{-1}`
``isobaric_mass_heat_capacity`` Specific heat capacity in
                                :math:`\sf J kg^{-1} K^{-1}`
``heat_of_vaporization``        Heat of vaporization in :math:`\sf J kg{^-1}`
=============================== ==============================================

Mathematical constants
======================

==========  ============
``golden``  Golden ratio
==========  ============

SI prefixes
===========

=========  ================
``yotta``  :math:`10^{24}`
``zetta``  :math:`10^{21}`
``exa``    :math:`10^{18}`
``peta``   :math:`10^{15}`
``tera``   :math:`10^{12}`
``giga``   :math:`10^{9}`
``mega``   :math:`10^{6}`
``kilo``   :math:`10^{3}`
``hecto``  :math:`10^{2}`
``deka``   :math:`10^{1}`
``deci``   :math:`10^{-1}`
``centi``  :math:`10^{-2}`
``milli``  :math:`10^{-3}`
``micro``  :math:`10^{-6}`
``nano``   :math:`10^{-9}`
``pico``   :math:`10^{-12}`
``femto``  :math:`10^{-15}`
``atto``   :math:`10^{-18}`
``zepto``  :math:`10^{-21}`
=========  ================

Non-SI ratios
=============

=======  =====================================
``ppm``  :math:`10^{-6}` `parts per million`
``ppb``  :math:`10^{-9}` `parts per billion`
``ppt``  :math:`10^{-12}` `parts per trillion`
=======  =====================================

Binary prefixes
===============

=================  ==============
``kibi``, ``KiB``  :math:`2^{10}`
``mebi``, ``MiB``  :math:`2^{20}`
``gibi``           :math:`2^{30}`
``tebi``           :math:`2^{40}`
``pebi``           :math:`2^{50}`
``exbi``           :math:`2^{60}`
``zebi``           :math:`2^{70}`
``yobi``           :math:`2^{80}`
=================  ==============

=================  ==============
``KB``             :math:`10^3`
``MB``             :math:`10^6`
=================  ==============

Earth characteristics
=====================

================  =====================================
``earth_mass``    Earth mass in :math:`\sf kg`
``earth_radius``  Earth radius in :math:`\sf m`
``atm``           Standard atmosphere in :math:`\sf Pa`
================  =====================================

"""
import numpy as np
import scipy.constants as spc


# Physical constants
g = earth_standard_gravity = spc.g  # m s^-2
h = planck = spc.Planck  # J s
k = boltzmann = spc.Boltzmann  # J K^-1
c = speed_of_light = spc.speed_of_light  # m s^-1
N_A = avogadro = N = spc.Avogadro  # mol^-1
K = zero_celsius = 273.15  # Kelvin at 0 Celsius
triple_point_water = 273.16  # Triple point temperature in K
R = gas_constant = spc.gas_constant  # J mol^-1 K^-1
mu_B = spc.e * spc.Planck * (0.25 / spc.pi) / spc.m_e  # J T^-1
molar_mass_dry_air = 28.9645e-3  # kg mol^-1
molar_mass_water = 18.01528e-3  # kg mol^-1
gas_constant_dry_air = R / molar_mass_dry_air  # J K^-1 kg^-1
gas_constant_water_vapor = R / molar_mass_water  # J K^-1 kg^-1
amu = spc.m_u
stefan_boltzmann_constant = 2 * np.pi**5 * k**4 / (15 * c**2 * h**3)
isobaric_mass_heat_capacity = 1003.5  # J kg^-1 K^-1
heat_of_vaporization = 2501000  # J kg^-1

# Mathematical constants
golden = golden_ratio = (1 + np.sqrt(5)) / 2

# SI prefixes
yotta = 1e24
zetta = 1e21
exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e9
mega = 1e6
kilo = 1e3
hecto = 1e2
deka = 1e1
deci = 1e-1
centi = 1e-2
milli = 1e-3
micro = 1e-6
nano = 1e-9
pico = 1e-12
femto = 1e-15
atto = 1e-18
zepto = 1e-21

# Non-SI ratios
ppm = 1e-6  # parts per million
ppb = 1e-9  # parts per billion
ppt = 1e-12  # parts per trillion

# Binary prefixes
kibi = KiB = 2**10
mebi = MiB = 2**20
gibi = 2**30
tebi = 2**40
pebi = 2**50
exbi = 2**60
zebi = 2**70
yobi = 2**80

KB = 10**3
MB = 10**6

# Earth characteristics
earth_mass = 5.97237e24  # kg
earth_radius = 6.3781e6  # m
atm = atmosphere = 101325  # Pa
