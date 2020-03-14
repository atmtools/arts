import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Basic settings of simulation
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
ws.StringSet(ws.iy_unit, "PlanckBT")
# Selecting the channels and viewing angles of sensor
ws.ArrayOfIndexCreate("channels")
ws.ArrayOfIndexCreate("viewing_angles")
ws.ArrayOfIndexSet(ws.channels, [-1])
# ArrayOfIndexSet(viewing_angles, [-1]) # Pick all viewing angles
# We use the usual zero-based ARTS indexing here, so index
# 11 is actually channel 12, etc.
# For testing, only calculate one viewing angle
ws.ArrayOfIndexSet(ws.viewing_angles, [0])
# Definition of sensor position
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 850000.0)
# Definition of sensor LOS (180 is Nadir looking direction of the sensor)
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 180.0)
# This index (met_mm_accuracy) selects the number of frequencies to be used
# from the available accuracies. This is set by the user in the main controlfile
# Must be documented here:
# N  Accuracy   Speed           Number of frequencies/channel
# 0: fast       very-very fast  1
# 1: normal     fast            varying
# 2: high       slow            varying
# 3: reference  extremely slow  varying
#
# More information on met_mm_accuracy can be found in the ARTS-Documentation
# web-page in the section Technical Reports.
# http://arts.mi.uni-hamburg.de/docs/met_mm_setup_documentation.pdf
#
ws.IndexCreate("met_mm_accuracy")
ws.IndexSet(ws.met_mm_accuracy, 1)
# ====================================================================
# Sensor characteristics
ws.execute_controlfile("instruments/metmm/sensor_descriptions/prepare_metmm.arts")
# Change the filename to switch sensors.
ws.execute_controlfile("instruments/metmm/sensor_descriptions/sensor_amsub.arts")
ws.execute_controlfile("instruments/metmm/sensor_descriptions/apply_metmm.arts")
# Common microwave sensor settings
ws.execute_controlfile("instruments/metmm/common_metmm.arts")
# ====================================================================
# Spectroscopy
ws.abs_speciesSet(
    species=[
        "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
        "O2-66, O2-CIAfunCKDMT100",
        "N2,  N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
        "O3",
    ]
)
# Read HITRAN catalog:
#    abs_linesReadFromHitran( abs_lines,
#                             "/scratch/uni/u237/data/catalogue/hitran/hitran2012_140407/HITRAN2012.par",
#                             0,
#                             1000e9 )
#    abs_linesArtscat5FromArtscat34( abs_lines )
ws.ReadARTSCAT(ws.abs_lines, "instruments/metmm/abs_lines_metmm.xml.gz")
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_lines_per_speciesSetCutoffForSpecies(
    option="ByLine", value=5000000000.0, species_tag="O3"
)
#    WriteXML( "zascii", abs_lines, "instruments/metmm/abs_lines_metmm.xml.gz" )
# ====================================================================
# Set surface reflectivity
# Reflectivity = 0.4; emissivity = 0.6
ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.4)
# Atmospheric profiles
ws.ReadXML(ws.batch_atm_fields_compact, "testdata/garand_profiles.xml.gz")
# add constant profiles for O2 and N2
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2", value=0.2095)
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2", value=0.7808)
# ====================================================================
# Absorption lookup table
ws.abs_lookupSetupBatch()
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.abs_lookupCalc()
# Setting the agenda for batch calculation
# Garand profiles have 42 different. We will make RT calculations for all of them.
@arts_agenda
def ybatch_calc_agenda(ws):
    # Extract the atmospheric profiles for this case:
    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    # Split up *atm_fields_compact* to
    # generate p_grid, t_field, z_field, vmr_field:
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
    # Optionally set Jacobian parameters.
    ws.jacobianOff()
    # No scattering
    ws.cloudboxOff()
    # get some surface properties from corresponding atmospheric fields
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)
    # Checks
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    # Perform RT calculations
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

# How many simulations do we want to perform?
# All atmospheres, or manually set the number (uncomment IndexSet line)
ws.nelemGet(ws.ybatch_n, ws.batch_atm_fields_compact)
# IndexSet(ybatch_n, 1)
# ====================================================================
# Execute the batch calculations:
# First check, then execute the batch RT calculations
ws.propmat_clearsky_agenda_checkedCalc()
ws.ybatchCalc()
# ====================================================================
# Store results
ws.WriteXML("ascii", ws.ybatch)
ws.WriteXML("ascii", ws.f_grid)
# Verify results
ws.ArrayOfVectorCreate("ybatch_ref")
ws.ReadXML(ws.ybatch_ref, "ybatchREFERENCE.xml")
ws.Compare(
    ws.ybatch, ws.ybatch_ref, 0.01, "Total BT should be close to the reference values"
)
