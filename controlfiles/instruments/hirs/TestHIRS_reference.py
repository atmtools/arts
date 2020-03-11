# DEFINITIONS:  -*-sh-*-
# --------------------------------------------------------------------
# Test batch calculations for HIRS.
# Viju Oommen John and Stefan Buehler, August to October 2008.
# --------------------------------------------------------------------

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

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
# blackbody surface with skin temperature interpolated from t_surface field
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.jacobianOff()
ws.StringCreate("satellite")
ws.ArrayOfIndexCreate("channels")
ws.ArrayOfIndexCreate("views")
ws.NumericCreate("f_grid_spacing")
# Select here which satellite you want
# ---
ws.StringSet(ws.satellite, "NOAA14")
# Select here which channels you want
# ---
#
# HIRS Channels 13-19 are shortwave channels. Simulating them with
# ARTS for thermal radiation only is pointless.
#
# WATCH OUT! We use the usual zero-based ARTS indexing here, so index
# 11 is actually channel 12, etc.
# ArrayOfIndexSet( channels, [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] )
ws.ArrayOfIndexSet(ws.channels, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
# Select here which views you want
# ---
ws.ArrayOfIndexSet(ws.views, [0, 7, 14, 21, 27])
# Set frequency grid spacing
# ---
# Concerning the frequency grid spacing:
# I tested the spacing on HIRS channel 12. With 5e8 Hz spacing I got
# an RMS error of 0.01 K and a maximum error of 0.02 K. This should be
# good enough. We will have 8402 frequencies for Channel 12.
# A coarser spacing of 5e9 Hz gave an RMS error of 0.03 and a maximum
# error 0.05.
# All these number are against a reference calculation with 5e7 Hz
# spacing, which was shown to give nearly identical result to a somewhat
# coarser calculation.
# Recommendation:
# reference calculation: f_grid_spacing 5e8
# faster calculations:   f_grid_spacing 5e9
ws.NumericSet(ws.f_grid_spacing, 500000000.0)
# Basic settings
# ---
# This example assumes 1D
ws.AtmosphereSet1D()
# Scalar RT
ws.IndexSet(ws.stokes_dim, 1)
# HIRS spectral response functions assume radiance.
ws.StringSet(ws.iy_unit, "1")
# Set up absorption
# ---
ws.abs_speciesSet(
    species=[
        "H2O, H2O-SelfContCKDMT100, H2O-ForeignContCKDMT100",
        "O3",
        "CO2, CO2-CKDMT100",
        "N2O",
        "CO",
        "CH4",
        "O2, O2-CIAfunCKDMT100",
        "N2, N2-CIAfunCKDMT100, N2-CIArotCKDMT100",
    ]
)
ws.abs_lineshapeDefine(ws.abs_lineshape, "Voigt_Kuntz6", "VVH", 750000000000.0)
# deriving abs_lines separate, because one might want to use a lookup table
#
# no use of HITRAN in test cases
# (uncomment when you want to use HITRAN)
ws.abs_linesReadFromHitran(
    ws.abs_lines, "HITRAN2012.par", 2002223400000.0, 798557330000000.0
)
# instead we use a HITRAN-extract converted to ARTSCAT
# (comment out when you want to use HITRAN)
# ReadXML( abs_lines, "testdata/abs_lines_IR.xml.gz" )
ws.abs_lines_per_speciesCreateFromLines()
# Sensor setup
# ---
# Definition of sensor position and LOS
ws.ReadXML(ws.sensor_los, "hirs.sensor_los.xml")
# Select those views that are requested by the user
ws.Select(ws.sensor_los, ws.sensor_los, ws.views)
ws.nrowsGet(ws.nrows, ws.sensor_los)
ws.ncolsGet(ws.ncols, ws.sensor_los)
ws.MatrixSetConstant(ws.sensor_pos, ws.nrows, ws.ncols, 850000.0)
# Start sensor response setup
# Normalise the sensor response
ws.IndexSet(ws.sensor_norm, 1)
ws.AntennaOff()
# Construct names of sensor description files
# Nominal channel frequencies
ws.StringCreate("f_backend_file")
ws.StringJoin(ws.f_backend_file, ws.satellite, "_HIRS.f_backend.xml")
# Weights associated with each frequency
ws.StringCreate("backend_channel_response_file")
ws.StringJoin(
    ws.backend_channel_response_file, ws.satellite, "_HIRS.backend_channel_response.xml"
)
# Spectrometer
ws.ReadXML(ws.f_backend, ws.f_backend_file)
ws.ReadXML(ws.backend_channel_response, ws.backend_channel_response_file)
# Select the desired channels
ws.Select(ws.f_backend, ws.f_backend, ws.channels)
ws.Select(ws.backend_channel_response, ws.backend_channel_response, ws.channels)
# Frequency grid
ws.f_gridFromSensorHIRS(
    ws.f_grid, ws.f_backend, ws.backend_channel_response, ws.f_grid_spacing
)
# Initialize sensor variables
ws.sensor_responseInit()
# Set up backend sensor response
ws.sensor_responseBackend()
# End of sensor response setup
# Compact line list, to kick out lines that are outside there
# cutoff range for all frequencies.
ws.abs_lines_per_speciesCompact()
# Atmospheric profiles
# ---
# Atmospheric profiles are stored in an ArrayOfGriddedField4.
# It contains one GriddedField4 for each atmospheric state.
#
ws.ReadXML(ws.batch_atm_fields_compact, "testdata/garand_profiles.xml.gz")
# add constant profiles for O2 and N2
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2", value=0.2095)
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2", value=0.7808)
# Set up lookup table
# ---
# Set parameters for lookup table
# Arguments omitted for better maintainability of this test file.
ws.abs_lookupSetupBatch()
# Optional manual override of T and VMR perturbations
#
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
ws.abs_xsec_agenda_checkedCalc()
ws.abs_lookupCalc()
# WriteXML( "binary", abs_lookup )
# Test (and print) lookup table accuracy
#
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
#
# The sensor altitude is predefined in amsu.arts to 850 km above the geoid.
# Comment out the next line if you want to set it to something else.
# MatrixSetConstant( sensor_pos, 1, 1, 850e3 )
# Set the agenda for batch calculations:
# ---
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
    # Perform RT calculations
    # ---
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.yCalc()
    # Optionally write out jacobian:
    # WriteXMLIndexed( output_file_format, ybatch_index, jacobian, "" )
    # Convert the measurement from radiance units to Planck Tb:
    ws.StringSet(ws.iy_unit, "PlanckBT")
    ws.yApplyUnit()


ws.ybatch_calc_agenda = ybatch_calc_agenda

# Set number of batch cases:
ws.nelemGet(ws.ybatch_n, ws.batch_atm_fields_compact)
# IndexSet( ybatch_n, 1 )
# Execute the batch calculations:
# ---
ws.propmat_clearsky_agenda_checkedCalc()
ws.ybatchCalc()
# Store result matrix:
# ---
ws.WriteXML("ascii", ws.ybatch)
# Compare with reference:
ws.ArrayOfVectorCreate("ybatch_ref")
ws.ReadXML(ws.ybatch_ref, "TestHIRS.NOAA14.ybatch.ref.xml")
ws.Compare(ws.ybatch_ref, ws.ybatch, 0.01)
