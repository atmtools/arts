# DEFINITIONS:  -*-sh-*-
# This control file handles a clearsky batch calculation including calculation of
# an absorption lookup table. Six selected atmospheres from Chevallier91L data
# are used for atmospheric input. No sensor characteristics are applied. Only two
# frequencies and two line of sights are used, for calculation time reasons.
#
# Author: Jana Mendrok
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/agendasDOIT.arts")
ws.execute_controlfile("general/planet_earth.arts")
# 1.General Settings:---------------------------------------------
# -----------------------------------------------------------------
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
# Set out put file format
# ------------------------
ws.output_file_formatSetAscii()
# Define f_grid
# --------------
ws.VectorSet(ws.f_grid, np.array([9.0e10, 1.9e11]))
# ReadXML(f_grid, "f_grid.xml")
# VectorNLinSpace (f_grid, 10, 120e9, 300e9)
# WriteXML (output_file_format, f_grid, "f_grid.xml")
# Set stokes dim
# --------------
ws.IndexSet(ws.stokes_dim, 1)
# def of atmosphere
# -----------------
ws.IndexSet(ws.atmosphere_dim, 1)
# Choose *y* output unit
# ----------------------
ws.StringSet(ws.iy_unit, "PlanckBT")
# No jacobian calculations
# -----------------
ws.jacobianOff()
# Modifiy the maximum propagation step, from the default(10.e3)
# to 250 m:---------------------------------------------------
ws.NumericSet(ws.ppath_lmax, 250.0)
# Surface properties
# -------------------
# Set surface reflectivity (=1-emissivity)
# corresponds to emissivity=0.75
ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.25)
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
)
# 2. Sensor:---------------------------------------------------------
# --------------------------------------------------------------------
# Definition of sensor position and LOS
# ------------------------------------
# Line of sight
ws.MatrixSet(ws.sensor_los, np.array([[131.0], [179.0]]))
# Sensor position
ws.nrowsGet(ws.nrows, ws.sensor_los)
ws.ncolsGet(ws.ncols, ws.sensor_los)
ws.MatrixSetConstant(ws.sensor_pos, ws.nrows, ws.ncols, 850000.0)
# No sensor characteristics are specified
ws.sensorOff()
# 3. Read chevallier atmospheric profiles for batch calc--------------
# ---------------------------------------------------------------------
ws.ReadXML(ws.batch_atm_fields_compact, "testdata/chevallierl91_all_extract.xml")
# add constant profiles for O2 and N2
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2", value=0.2095)
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2", value=0.7808)
# 4. Absorption-------------------------------------------------
# ---------------------------------------------------------------
ws.abs_speciesSet(species=["H2O-PWR98", "O3", "O2-PWR93", "N2-SelfContStandardType"])
# Creation of abs_lookup table
# -----------------------------
ws.ReadSplitARTSCAT(basename="spectroscopy/Perrin/", fmin=0.0, fmax=3000000000000.0)
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_lookupSetupBatch()
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.abs_lookupCalc()
ws.abs_lookupAdapt()
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# no scattering
ws.cloudboxOff()
# Batch Agenda-----------------------------------------------------------
# ------------------------------------------------------------------------
@arts_agenda
def ybatch_calc_agenda(ws):
    # Extract the atmospheric profiles for current atmosphere:
    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    # Split up *atm_fields_compact* to
    # generate p_grid, t_field, z_field, massdensity_field, vmr_field:
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
    # Get some surface properties from corresponding atmospheric fields
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)
    # Consistency checks
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    # Calculate complete measurement vector
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

# Set number of batch cases:
ws.nelemGet(ws.ybatch_n, ws.batch_atm_fields_compact)
# IndexSet(ybatch_start, 0)
# IndexSet(ybatch_n, 4)
# ===========start batch calc=====================
# Execute the batch calculations:
# This test control file can be run multi-threaded, since it was approved
# that none of the jobs fails.
# If settings are changed, especially if the input atmospheres are altered
# or exchanged, the robust option in *ybatchCalc* should be used.
ws.propmat_clearsky_agenda_checkedCalc()
ws.ybatchCalc()
# Call *ybatchCalc* robust:
# Set robust flag to 1. If one individual job fails, ARTS will continue
# with the next batch job.
# ybatchCalc( robust=1 )
# WriteXML( in=ybatch, filename="TestBatch.ybatch.ref.xml" )
# Verify results
ws.ArrayOfVectorCreate("ybatch_ref")
ws.ReadXML(out=ws.ybatch_ref, filename="TestBatch.ybatch.ref.xml")
ws.Compare(ws.ybatch, ws.ybatch_ref, 1e-06)
# ==================stop==========================
# End of Main
