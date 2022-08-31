#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Manfred Brath

"""

import numpy as np

from pyarts.workspace import Workspace, arts_agenda
from pyarts.constants import c  # import speed of light

# =============================================================================
# paths/constants
# =============================================================================


# geographical position
geo_pos = [0., 0.]

# position where the sun is in zenith
sun_pos = [0., 0]

# Set frequency
wavelengths = np.linspace(700e-9, 600e-9, 20)
f_grid = c / wavelengths

# surface reflectivity
reflectivity = 0.2

# Recalculate lookup table flag
# Set to True, if you ant to recalculate the lookup table
recalc = False

# =============================================================================
# open workspace
# =============================================================================

ws = Workspace()
ws.verbositySetScreen(level=2)
ws.verbositySetAgenda(level=0)

# select/define agendas
# =============================================================================


ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/planet_earth.arts")


# gas scattering agenda
@arts_agenda(ws=ws, set_agenda=True)
def gas_scattering_agenda(ws):
    ws.Ignore(ws.rtp_vmr)
    ws.gas_scattering_coefAirSimple()
    ws.gas_scattering_matRayleigh()


ws.iy_main_agendaSet( option="Clearsky" )
ws.iy_space_agendaSet( option="CosmicBackground" )
ws.ppath_step_agendaSet( option="GeometricPath" )
ws.ppath_agendaSet( option="FollowSensorLosPath" )

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
ws.VectorNLogSpace(ws.p_grid, 81, 1013e2, 1)

# set geographical position
ws.lat_true = [geo_pos[0]]
ws.lon_true = [geo_pos[1]]

ws.AtmosphereSet1D()

# define absorption
# =============================================================================

# set absorption species
ws.abs_speciesSet(species=["H2O, H2O-SelfContCKDMT320, H2O-ForeignContCKDMT320"])

if recalc == False:
    try:
        ws.Touch(ws.abs_lines_per_species)
        ws.ReadXML(ws.abs_lookup, "TestDisortCalcIrradiance.abs_lookup.xml")
        ws.abs_lookupAdapt()
        ws.lbl_checked = 1

    except RuntimeError:
        recalc = True

if recalc == True:
    # recalc LUT
    ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")
    ws.abs_lines_per_speciesCutoff(option="ByLine", value=750e9)
    ws.abs_lookupSetupWide(t_min=150., t_max=320., p_step=0.1)
    ws.propmat_clearsky_agendaAuto()
    ws.lbl_checkedCalc()
    ws.abs_lookupCalc()
    ws.WriteXML('binary', ws.abs_lookup, "TestDisortCalcIrradiance.abs_lookup.xml")

ws.propmat_clearsky_agendaAuto(use_abs_lookup=1)

# define atmosphere
# =============================================================================

# Atmospheric profiles
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()

# Get ground altitude (z_surface) from z_field
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.)

# set surface temperature
ws.surface_skin_t = ws.t_field.value[0, 0, 0]

# set surface reflectivity
ws.surface_scalar_reflectivity = [reflectivity]

# set star source
# ws.Touch(ws.stars)
ws.starsAddSingleBlackbody(latitude=sun_pos[0], longitude=sun_pos[1])

# set cloudbox to full atmosphere
ws.cloudboxSetFullAtm()

# No jacobian calculations
ws.jacobianOff()

# set particle scattering
# =============================================================================

# setting here particle basically to zero
ws.scat_data_checked = 1
ws.Touch(ws.scat_data)
ws.pnd_fieldZero()

# =============================================================================
# the calculation
# =============================================================================

# Switch on/off gas scattering
ws.IndexSet(ws.gas_scattering_do, 1)

# Check model atmosphere
ws.scat_data_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()

ws.DisortCalcIrradiance(emission=0)

# calculate irradiance (flux)
ws.RadiationFieldSpectralIntegrate(ws.irradiance_field, f_grid, ws.spectral_irradiance_field)

# Uncomment, if you want to save new reference data
# ws.WriteXML('ascii', ws.irradiance_field, "TestDisortCalcIrradiance.irradiance_fieldREFERENCE.xml")

# Do the check
ws.Tensor4Create("irradiance_field0")
ws.ReadXML(ws.irradiance_field0, "TestDisortCalcIrradiance.irradiance_fieldREFERENCE.xml")
ws.Compare(ws.irradiance_field, ws.irradiance_field0, 1e-6)
