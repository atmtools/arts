# DEFINITIONS:  -*-sh-*-
# --------------------------------------------------------------------
# !!! WARNING !!! Untested !!! USE WITH CARE !!!
# Test batch calculations for SEVIRI.
# Based on HIRS test by Viju Oommen John and Stefan Buehler, August to
# October 2008.
# --------------------------------------------------------------------

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Copy( ppath_step_agenda,ppath_step_agenda__RefractedPath )
# Copy( refr_index_air_agenda, refr_index_air_agenda__NoRefrac)
# blackbody surface with skin temperature interpolated from t_surface field
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.StringCreate("satellite")
ws.ArrayOfIndexCreate("channels")
ws.ArrayOfIndexCreate("views")
ws.StringCreate("hitran_file")
ws.NumericCreate("f_grid_spacing")
# Select here which satellite you want
# ---
ws.StringSet(ws.satellite, "MET9")
# Select here which channels you want
# ---
#
# MVIRI Channel-1 is WV2
#              -2 is IR1
#              -3 is VIS
#
# WATCH OUT! We use the usual zero-based ARTS indexing here, so index
# 1 is actually channel 2, etc.
ws.ArrayOfIndexSet(ws.channels, [3, 4, 5, 6, 7, 8, 9, 10, 11])
# Select here which views you want
# ---
ws.ArrayOfIndexSet(
    ws.views,
    [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ],
)
# You have to give here the location of the HITRAN catalogue
# ---
# FIXME: This has to be replaced by a little piece of the catalog
# explicitly included here, or?
# We have to discuss with Oliver (and perhaps Patrick) how to handle
# this.
# We can also include a pre-calculated absorption table for HIRS with
# ARTS.
# StringSet(hitran_file,"/storage3/data/catalogue/hitran/hitran2004/HITRAN04.par")
ws.StringSet(
    ws.hitran_file,
    "/scratch/uni/u237/data/catalogue/hitran/hitran2012_140407/HITRAN2012.par",
)
# StringSet(hitran_file,"/home/patrick/Data/HITRAN_2004/HITRAN04.par")
# Set frequency grid spacing
# (See comments in hirs_reference.arts concerning useful values here)
# ---
ws.StringCreate("f_grid_spacing_str")
ws.NumericSet(ws.f_grid_spacing, 6000000000.0)
# default 5e8
# StringSet(f_grid_spacing_str,"6e8")
ws.StringSet(ws.f_grid_spacing_str, "6e9_fast")
# Basic settings (already needed in sensor part)
# ---
# This example assumes 1D
ws.AtmosphereSet1D()
# scalar RT
ws.IndexSet(ws.stokes_dim, 1)
ws.execute_controlfile("seviri_fast.arts")
# INCLUDE "seviri_reference.arts"
# Set up absorption
# =================
# Atmospheric profiles
ws.ReadXML(ws.batch_atm_fields_compact, "testdata/garand_profiles.xml.gz")
# add constant profiles for O2 and N2
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2", value=0.2095)
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2", value=0.7808)
# Set parameters for lookup table
# ---
# Arguments omitted for better maintainability of this test file.
ws.abs_lookupSetupBatch()
# Optional manual override of T and VMR perturbations
# ---
# If your input data contains extreme outliers, the lookup table will
# get unreasonably large. It is suggested that you instead set them
# manually here. The Chevallier 91L data (humidity set) contains
# temperature variations from -70 to +55 (K) and humidity variations from
# 0 to 6 (fractional units). This should be the reasonable range of
# atmospheric variability. You will get failures from individual jobs in
# the batch, that are outside this range.
# VectorLinSpace( abs_t_pert, -70, 55, 5 )
# VectorLinSpace( abs_nls_pert, 0, 6, 0.5 )
# Create the lookup table
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.abs_lookupCalc()
# WriteXML("binary", abs_lookup)
# Test (and print) lookup table accuracy
# ---
# This method is not necessary, it just tests and prints the lookup
# table accuracy. Comment it out if you do not want this
# information. The test should take approximately as much time as the
# lookup table generation, perhaps even more. So, it is not cheap!
# abs_lookupTestAccuracy
# Set propmat_clearsky_agenda to use lookup table
# ---
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Set up RT calculation
# =====================
# Definition of sensor position and LOS
# ---
# Optionally set sensor_pos
# ---
# The sensor altitude is predefined in amsu.arts to 850 km above the geoid.
# Comment out the next line if you want to set it to something else.
ws.MatrixSetConstant(ws.sensor_pos, 26, 1, 36000000.0)
# Set the agenda for batch calculations:
# ---
#
@arts_agenda
def ybatch_calc_agenda(ws):
    # Extract the atmospheric profiles for this case:
    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    # Split up *atm_fields_compact* to
    # generate p_grid, t_field, z_field, vmr_field:
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
    # Optionally set Jacobian parameters.
    # uncomment this for NO jacobian calculations
    ws.jacobianOff()
    # Uncomment this block if you want Jacobians. Attention, it slows down the
    # computation a lot.
    # Also, you can add other Jacobians here, for example for temperature.
    #
    #   jacobianInit
    #   jacobianAddAbsSpecies( jacobian_quantities, jacobian_agenda,
    #                          atmosphere_dim,
    #                          p_grid, lat_grid, lon_grid,
    #                          p_grid, lat_grid, lon_grid,
    #                          "H2O, H2O-SelfContCKDMT100, H2O-ForeignContCKDMT100",
    #                          "rel")
    #   jacobianClose
    # No scattering
    ws.cloudboxOff()
    # get some surface properties from corresponding atmospheric fields
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)
    # sensorOff
    # Perform RT calculations
    # ---
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.yCalc()
    # Optionally write out jacobian:
    #  WriteXMLIndexed( output_file_format, ybatch_index, jacobian, "" )
    ws.StringSet(ws.iy_unit, "PlanckBT")
    ws.yApplyUnit()


ws.ybatch_calc_agenda = ybatch_calc_agenda

# Set number of batch cases:
ws.nelemGet(ws.ybatch_n, ws.batch_atm_fields_compact)
# IndexSet(ybatch_n, 1)
# Execute the batch calculations:
# ---
ws.propmat_clearsky_agenda_checkedCalc()
ws.StringCreate("out_file_name_sensor_response")
ws.StringCreate("out_file_name_ybatch")
ws.StringCreate("out_file_name_f_grid")
ws.StringSet(ws.out_file_name_sensor_response, "TestSEVIRI.sensor_response_")
ws.StringSet(ws.out_file_name_ybatch, "TestSEVIRI.ybatch_")
ws.StringSet(ws.out_file_name_f_grid, "TestSEVIRI.f_grid_")
ws.Append(ws.out_file_name_sensor_response, ws.satellite)
ws.Append(ws.out_file_name_ybatch, ws.satellite)
ws.Append(ws.out_file_name_f_grid, ws.satellite)
ws.StringCreate("extent")
ws.StringSet(ws.extent, "_")
ws.StringSet(ws.dummy, ".xml")
ws.Append(ws.extent, ws.f_grid_spacing_str)
ws.Append(ws.extent, ws.dummy)
ws.Append(ws.out_file_name_sensor_response, ws.extent)
ws.Append(ws.out_file_name_ybatch, ws.extent)
ws.Append(ws.out_file_name_f_grid, ws.extent)
ws.WriteXML("ascii", ws.sensor_response, ws.out_file_name_sensor_response)
ws.ybatchCalc()
# Store result matrix:
# ---
ws.WriteXML("ascii", ws.ybatch, ws.out_file_name_ybatch)
ws.WriteXML("ascii", ws.f_grid, ws.out_file_name_f_grid)
# Verify results
ws.ArrayOfVectorCreate("ybatch_ref")
ws.ReadXML(ws.ybatch_ref, "TestSEVIRI.ybatch_MET9_6e9_fastREFERENCE.xml")
ws.Compare(
    ws.ybatch, ws.ybatch_ref, 0.01, "Total BT should be close to the reference values"
)
