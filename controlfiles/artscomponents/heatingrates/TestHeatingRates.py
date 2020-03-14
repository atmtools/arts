# DEFINITIONS:  -*-sh-*-

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# ==========================================================================
####  Include files ####
# ==========================================================================
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# ==========================================================================
####  Set Agendas ####
# ==========================================================================
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# Surface Agenda, see in python script
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Set propmat_clearsky_agenda to use lookup table
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# saturation pressure agenda of water vapor
@arts_agenda
def water_p_eq_agenda(ws):
    ws.water_p_eq_fieldMK05()


ws.water_p_eq_agenda = water_p_eq_agenda


@arts_agenda
def iy_cloudbox_agenda(ws):
    ws.iyInterpCloudboxField()


ws.iy_cloudbox_agenda = iy_cloudbox_agenda

# surface agenda
ws.Copy(
    ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_surface
)
# ==========================================================================
#### Basic/Sensor settings ####
# ==========================================================================
ws.StringSet(ws.iy_unit, "1")
# Dimension / type of atmosphere
ws.AtmosphereSet1D()
# Dimension of Stokes vector
ws.IndexSet(ws.stokes_dim, 1)
# Sensor Position
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 200000.0)
# set frequency grid
ws.VectorNLinSpace(ws.f_grid, 10, 300000000000.0, 30000000000000.0)
# set angular grid
# N_za_grid is the number zenith angles, which needs to be an even number.
ws.AngularGridsSetFluxCalc(N_za_grid=6, N_aa_grid=1, za_grid_type="double_gauss")
# ==========================================================================
#### Set up atmosphere ####
# ==========================================================================
##read atm data in ARTS
ws.ReadXML(ws.batch_atm_fields_compact, "testdata/garand_profiles.xml.gz")
# as the atm data has no O2, N2, and CO2 add it
ws.batch_atm_fields_compactAddConstant(
    name="abs_species-O2", value=0.2095, prepend=0, condensibles=["abs_species-H2O"]
)
ws.batch_atm_fields_compactAddConstant(
    name="abs_species-N2", value=0.7808, prepend=0, condensibles=["abs_species-H2O"]
)
ws.batch_atm_fields_compactAddConstant(
    name="abs_species-CO2",
    value=0.00039755,
    prepend=0,
    condensibles=["abs_species-H2O"],
)
# ==========================================================================
#### Set up absorption ####
# ==========================================================================
ws.abs_speciesSet(
    species=[
        "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
        "O3",
        "O2, O2-CIAfunCKDMT100",
        "CO2, CO2-CKDMT252",
        "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
        "CH4",
        "CO",
    ]
)
# Read HITRAN catalog:
# abs_linesReadFromHitran(abs_lines,
# 	"HITRAN2012.par",
# 	2e11,4e13
# )
#
# abs_lines_per_speciesCreateFromLines
# abs_lines_per_speciesCompact
#
# #Calculate absorption lookup table
# abs_lookupSetupBatch
# abs_xsec_agenda_checkedCalc
# abs_lookupCalc
#
# WriteXML( "binary", abs_lines_per_species)
# WriteXML( "binary", abs_lookup)
# read absorption lookup table
ws.ReadXML(ws.abs_lookup, "TestHeatingRates.abs_lookup.xml")
ws.abs_lookupAdapt()
# ==========================================================================
#### RT calculations ####
# ==========================================================================
# Here we just take the first garand atmosphere
ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, 0)
# Split up *atm_fields_compact* to generate p_grid, t_field, z_field, vmr_field:
ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
# Set surface altitude
ws.Extract(ws.z_surface, ws.z_field, 0)
# Set surface skin temperature
ws.Extract(ws.t_surface, ws.t_field, 0)
# set jacobian,sensor and cloudbox off
ws.jacobianOff()
ws.cloudboxOff()
ws.sensorOff()
# Consistency checks
ws.atmgeom_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.cloudbox_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
# calculate intesity field
ws.Tensor3Create("trans_field")
# calculate the full radiation field for plane parallel atmosphere
ws.spectral_radiance_fieldClearskyPlaneParallel(trans_field=ws.trans_field)
# Calculate Radiance from cloudbox_field
ws.RadiationFieldSpectralIntegrate(
    radiation_field=ws.radiance_field,
    spectral_radiation_field=ws.spectral_radiance_field,
)
# calculate flux from radiance field
ws.irradiance_fieldFromRadiance()
# set specific heat capacity and set earth gravity, which are needed as input
# for the heating rate calculations.
ws.npagesGet(ws.npages, ws.t_field)
ws.nrowsGet(ws.nrows, ws.t_field)
ws.ncolsGet(ws.ncols, ws.t_field)
ws.Tensor3SetConstant(ws.specific_heat_capacity, ws.npages, ws.nrows, ws.ncols, 1006.0)
ws.NumericSet(ws.g0, 9.80665)
# ARTS heating rates
ws.heating_ratesFromIrradiance()
ws.Tensor3Create("heating_rates1")
ws.Copy(ws.heating_rates1, ws.heating_rates)
# For testing calculate first spectral irradiance and than do spectral
# integration
ws.spectral_irradiance_fieldFromSpectralRadianceField()
ws.RadiationFieldSpectralIntegrate(
    radiation_field=ws.irradiance_field,
    spectral_radiation_field=ws.spectral_irradiance_field,
)
ws.heating_ratesFromIrradiance()
ws.Compare(ws.heating_rates, ws.heating_rates1, 1e-14)
# ==========================================================================
#### Compare results ####
# ==========================================================================
ws.Tensor3Create("heating_rates0")
ws.ReadXML(ws.heating_rates0, "TestHeatingRates.heating_ratesREFERENCE.xml")
ws.Compare(ws.heating_rates, ws.heating_rates0, 1e-09)
