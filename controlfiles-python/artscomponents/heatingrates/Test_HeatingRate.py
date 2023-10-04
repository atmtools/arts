#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Manfred Brath

"""

import numpy as np
from pyarts.workspace import Workspace

# =============================================================================
# paths/constants
# =============================================================================


# geographical position
geo_pos = [0., 0.]

# set frequency grid
f_grid = np.linspace(3e11, 3e13, 10)

# surface reflectivity
reflectivity = 0.0

# Number of streams
Nstreams = 6

# specific heat capacity of dry air
Cp = 1.0035e+03  # [Jkg^{-1}K^{-1}]

# Recalculate lookup table flag
# Set to True, if you ant to recalculate the lookup table
recalc = False

# numper of pressure levels
npres = 81

# =============================================================================
# open workspace
# =============================================================================

ws = Workspace()
ws.verbositySetScreen(level=2)
ws.verbositySetAgenda(level=0)

# select/define agendas
# =============================================================================


ws.PlanetSet(option="Earth")
ws.gas_scattering_agendaSet()

# define environment
# =============================================================================

# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)

# Reference ellipsoid
ws.refellipsoidEarth(ws.refellipsoid, "Sphere")

# Frequency grid
ws.f_grid = f_grid

# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
ws.VectorNLogSpace(ws.p_grid, npres, 1013e2, 1)

# set geographical position
ws.lat_true = [geo_pos[0]]
ws.lon_true = [geo_pos[1]]

ws.AtmosphereSet1D()

# define absorption
# =============================================================================

# set absorption species
ws.abs_speciesSet(species=["H2O, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350",
                           "O3",
                           "O2-*-1e12-1e99,O2-CIAfunCKDMT100",
                           "CO2, CO2-CKDMT252",
                           "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252"])

if recalc == False:
    try:
        ws.Touch(ws.abs_lines_per_species)
        ws.ReadXML(ws.abs_lookup, "Test_HeatingRate.abs_lookup.xml")
        ws.abs_lookupAdapt()
        ws.lbl_checked = 1

    except RuntimeError:
        recalc = True

if recalc == True:
    # recalc LUT
    ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")
    ws.abs_lines_per_speciesCutoff(option="ByLine", value=750e9)
    ws.abs_lines_per_speciesTurnOffLineMixing()
    ws.abs_lookupSetupWide(t_min=150., t_max=320., p_step=0.1)
    ws.propmat_clearsky_agendaAuto()
    ws.lbl_checkedCalc()
    ws.abs_lookupCalc()
    ws.WriteXML('binary', ws.abs_lookup,
                "Test_HeatingRate.abs_lookup.xml")

ws.propmat_clearsky_agendaAuto(use_abs_lookup=1)

# define atmosphere
# =============================================================================

# Atmospheric profiles
ws.AtmRawRead(basename="testdata/tropical")
# ws.AtmRawRead(basename="../../arts_dev/arts-xml-data/planets/Earth/Fascod/tropical/tropical")
ws.AtmFieldsCalc()

# Get ground altitude (z_surface) from z_field
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.)

# set surface temperature
ws.surface_skin_t = ws.t_field.value[0, 0, 0]

# set surface reflectivity
ws.surface_scalar_reflectivity = [reflectivity]

# set sun to off
ws.sunsOff()

# set cloudbox to full atmosphere
ws.cloudboxSetFullAtm()

# No jacobian calculations
ws.jacobianOff()

# Switch off gas scattering
ws.gas_scatteringOff()

# set particle scattering
# =============================================================================

# setting here particles to zero
ws.scat_data_checked = 1
ws.Touch(ws.scat_data)
ws.pnd_fieldZero()

# =============================================================================
# the radiative transfer
# =============================================================================

# Check model atmosphere
ws.scat_data_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()

ws.spectral_irradiance_fieldDisort(nstreams=Nstreams)

# calculate irradiance (flux)
ws.RadiationFieldSpectralIntegrate(
    ws.irradiance_field, ws.f_grid, ws.spectral_irradiance_field)

# =============================================================================
# Calculate heating rates accurately with varying gravity
# =============================================================================
# The WSV (mass) specific_heat_capacity can be set to vary with pressure, latitude and
# longitude as it is a Tensor3. But we set it here to constant value
ws.Tensor3SetConstant(ws.specific_heat_capacity, npres, 1, 1, Cp)

# calculate heating rate varying gravity
ws.heating_ratesFromIrradiance()

# Uncomment, if you want to save new reference data
# ws.WriteXML('ascii', ws.heating_rates, "Test_HeatingRate.heating_ratesREFERENCE.xml")

# Do the check
ws.Tensor3Create("heating_rates0")
ws.ReadXML(ws.heating_rates0, "Test_HeatingRate.heating_ratesREFERENCE.xml")
ws.Compare(ws.heating_rates, ws.heating_rates0, 1e-6)


# =============================================================================
# Calculate heating rates the simple way
# =============================================================================


# calculate heating rate varying gravity
ws.heating_ratesFromIrradianceSimple(mass_specific_heat_capacity=Cp, gravity=9.80665)

# Uncomment, if you want to save new reference data
# ws.WriteXML('ascii', ws.heating_rates,
#             "Test_HeatingRate.heating_rates_simpleREFERENCE.xml")

# Do the check
ws.ReadXML(ws.heating_rates0, "Test_HeatingRate.heating_rates_simpleREFERENCE.xml")
ws.Compare(ws.heating_rates, ws.heating_rates0, 1e-6)

