#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Manfred Brath

This script test the cloudbox_fieldInterp2Azimuth wsm.

In ARTS a 1D atmosphere cannot have a azimuth dependency, but if a 
collimated source like a sun is present even a 1D atmosphere has an 
azimuth dependency. 
To overcome this constraint when using yCalc, the user must set an
additional local sensor line of sight azimuth angle for the true
geopgraphical location of the atmosphere. For this angle the 
cloudbox_field with azimuthal dependency is interpolated to a 
cloudbox_field without azimuthal dependency so that it can be by yCalc.

"""

import numpy as np

from pyarts.workspace import Workspace, arts_agenda
from pyarts.constants import c

# constants
# =============================================================================

# geographical position
geo_pos = [0.0, 0.0]

# position where the sun is in zenith
sun_pos = [0.0, 75]

# Set frequency
wavelength = np.array([546]) * 1e-9
f_grid = c / wavelength

# surface reflectivity
reflectivity = 0.2

# local sensor line of sight azimuth angle
local_los_azimuth_angle = -180.0

# Recalculate lookup table flag
# Set to True, if you want to recalculate the lookup table
recalc = False

# Create new reference
# Set to True, if you want to recalculate a new reference
new_reference = False

# open workspace
# =============================================================================

ws = Workspace()
ws.verbositySetScreen(level=2)
ws.verbositySetAgenda(level=0)

# select/define agendas
# =============================================================================

ws.PlanetSet(option="Earth")


# gas scattering agenda
@arts_agenda(ws=ws, set_agenda=True)
def gas_scattering_agenda(ws):
    ws.Ignore(ws.rtp_vmr)
    ws.gas_scattering_coefAirSimple()
    ws.gas_scattering_matRayleigh()


ws.iy_main_agendaSet(option="Clearsky")
ws.iy_space_agendaSet(option="CosmicBackground")

# water agenda
ws.water_p_eq_agendaSet()
ws.ppath_step_agendaSet(option="GeometricPath")
ws.ppath_agendaSet(option="FollowSensorLosPath")

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
ws.abs_speciesSet(species=["H2O, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350"])

if recalc == False:
    try:
        ws.Touch(ws.abs_lines_per_species)
        ws.ReadXML(ws.abs_lookup, "Testcloudbox_fieldInterp2Azimuth.abs_lookup.xml")
        ws.abs_lookupAdapt()
        ws.lbl_checked = 1

    except RuntimeError:
        recalc = True

if recalc == True:
    ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")
    ws.abs_lines_per_speciesCutoff(option="ByLine", value=750e9)
    ws.abs_lines_per_speciesTurnOffLineMixing()
    ws.abs_lookupSetupWide(t_min=180.0, t_max=300.0, p_step=0.1)
    ws.propmat_clearsky_agendaAuto()
    ws.lbl_checkedCalc()
    ws.abs_lookupCalc()
    ws.WriteXML(
        "binary", ws.abs_lookup, "Testcloudbox_fieldInterp2Azimuth.abs_lookup.xml"
    )

ws.propmat_clearsky_agendaAuto(use_abs_lookup=1)

# set angular grid
ws.DOAngularGridsSet(N_za_grid=20, N_aa_grid=41)

# define atmosphere
# =============================================================================

# Atmospheric profiles
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()

# Get ground altitude (z_surface) from z_field
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)

# set surface temperature
ws.surface_skin_t = ws.t_field.value[0, 0, 0]

# set surface reflectivity
ws.surface_scalar_reflectivity = [reflectivity]

# set sun source
# ws.Touch(ws.sun)
ws.sunsAddSingleBlackbody(latitude=sun_pos[0], longitude=sun_pos[1])

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

# the calculation
# =============================================================================


# Switch on/off gas scattering
ws.IndexSet(ws.gas_scattering_do, 1)

# Switch on/off sun
ws.IndexSet(ws.suns_do, 1)

# Check model atmosphere
ws.scat_data_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()

ws.DisortCalc(nstreams=20, intensity_correction=0)
ws.cloudbox_fieldInterp2Azimuth(local_los_azimuth_angle=local_los_azimuth_angle)

if new_reference:
    ws.WriteXML(
        "ascii",
        ws.cloudbox_field,
        "Testcloudbox_fieldInterp2Azimuth.cloudbox_fieldREFERENCE.xml",
    )

# Do the check
ws.Tensor7Create("cloudbox_field0")
ws.ReadXML(
    ws.cloudbox_field0, "Testcloudbox_fieldInterp2Azimuth.cloudbox_fieldREFERENCE.xml"
)
ws.Compare(ws.cloudbox_field, ws.cloudbox_field0, 1e-6)
