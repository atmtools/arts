# DEFINITIONS:  -*-sh-*-
# Simple simulations of ground-based measurements of ozone at 110.8 GHz,
# mainly to test different sensor response methods.
# User defined variables are also used heavily, to set up consistent
# frequency grids inside the control file (instead of loading files).
# This includes the generation of a compact, non-uniform, f_grid.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Number of Stokes components to be computed
ws.IndexSet(ws.stokes_dim, 1)
# ---- f_grid ----------------------------------------------------------------
ws.NumericCreate("v0")
ws.NumericCreate("fw")
ws.IndexCreate("nlogpart")
ws.NumericCreate("fw_fine")
ws.NumericCreate("df_fine")
# Centre frequency
ws.NumericSet(ws.v0, 110836040000.0)
# One sided width of f_grid
ws.NumericSet(ws.fw, 330000000.0)
# Numer of points (on each side) of logarithmic part
ws.IndexSet(ws.nlogpart, 35)
# One sided width of fine grid at centre of f_grid
ws.NumericSet(ws.fw_fine, 240000.0)
# Spacing of this fine grid
ws.NumericSet(ws.df_fine, 40000.0)
# A logarithmically spaced grid between [fw_fine,fw]
ws.NumericCreate("f1")
ws.NumericCreate("f2")
ws.VectorCreate("flog")
ws.Copy(ws.f1, ws.fw_fine)
ws.Copy(ws.f2, ws.fw)
ws.VectorNLogSpace(ws.flog, ws.nlogpart, ws.f1, ws.f2)
# First part of f_grid is flog "mirrored"
ws.VectorFlip(ws.f_grid, ws.flog)
ws.VectorScale(ws.f_grid, ws.f_grid, -1.0)
# Append an equidistant grid between [-fw_fine+df_fine,fw_fine-df_fine]
ws.VectorCreate("feqd")
ws.Copy(ws.f1, ws.fw_fine)
ws.NumericScale(ws.f1, ws.f1, -1.0)
ws.NumericAdd(ws.f1, ws.f1, ws.df_fine)
ws.NumericScale(ws.f2, ws.f1, -1.0)
ws.VectorLinSpace(ws.feqd, ws.f1, ws.f2, ws.df_fine)
ws.Append(ws.f_grid, ws.feqd)
# Append flog
ws.Append(ws.f_grid, ws.flog)
# Add v0
ws.VectorAddScalar(ws.f_grid, ws.f_grid, ws.v0)
# ---- Species ---------------------------------------------------------------
ws.abs_speciesSet(species=["O3", "H2O"])
# ---- Atmospheric scenario --------------------------------------------------
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# A pressure grid rougly matching 0 to 80 km in 500 m steps.
ws.VectorNLogSpace(ws.p_grid, 160, 101300.0, 1.0)
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
# ---- Absorption ------------------------------------------------------------
ws.ReadARTSCAT(ws.abs_lines, "testdata/ozone_line.xml")
ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", 750000000000.0)
ws.abs_linesSetNormalization(ws.abs_lines, "VVH")
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_cont_descriptionInit()
ws.AbsInputFromAtmFields()
ws.abs_speciesSet(abs_species=ws.abs_nls, species=[])
ws.VectorSet(ws.abs_nls_pert, array([], dtype=float64))
ws.VectorSet(ws.abs_t_pert, array([], dtype=float64))
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.jacobianOff()
ws.abs_lookupCalc()
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# ---- Observation geometry --------------------------------------------------
ws.NumericCreate("z_platform")
ws.NumericCreate("za")
# Platform altitude
ws.NumericSet(ws.z_platform, 50.0)
# Zenith angle
ws.NumericSet(ws.za, 60.0)
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, ws.z_platform)
ws.Copy(ws.z_surface, ws.sensor_pos)
ws.MatrixSetConstant(ws.sensor_los, 1, 1, ws.za)
# ---- Finalise atmosphere ---------------------------------------------------
ws.cloudboxOff()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
# Important to run HSE before yCalc if temperature jacobians with HSE
# will be used. A latitude and longitude must here be specified.
ws.Extract(ws.p_hse, ws.p_grid, 0)
ws.NumericSet(ws.z_hse_accuracy, 0.1)
ws.VectorSet(ws.lat_true, array([58.0]))
ws.VectorSet(ws.lon_true, array([-12.0]))
#
ws.z_fieldFromHSE()
# ---- Turn off cosmic background radiation  ---------------------------------
# This to faciliate comparison of spectra from the different observation modes
@arts_agenda
def iy_space_agenda(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.nelemGet(ws.nelem, ws.f_grid)
    ws.MatrixSetConstant(ws.iy, ws.nelem, ws.stokes_dim, 0.0)


ws.iy_space_agenda = iy_space_agenda

# --- Common sensor settings -------------------------------------------------
ws.FlagOn(ws.sensor_norm)
ws.StringSet(ws.iy_unit, "RJBT")
ws.NumericCreate("f_resolution")
ws.NumericCreate("f_switch")
# Resolution (and also channel spacing) of spectrometer
ws.NumericSet(ws.f_resolution, 500000.0)
# Size of frequency throw
ws.NumericSet(ws.f_switch, 10000000.0)
ws.VectorCreate("f_resolution_v")
ws.VectorSetConstant(ws.f_resolution_v, 1, ws.f_resolution)
ws.backend_channel_responseGaussian(
    ws.backend_channel_response, ws.f_resolution_v, array([2.0])
)
# Calculate where first channel can start (considering f_switch and
# channel widths)
ws.Copy(ws.f1, ws.fw)
ws.NumericScale(ws.f1, ws.f1, -1.0)
ws.NumericAdd(ws.f1, ws.f1, ws.f_switch)
ws.NumericAdd(ws.f1, ws.f1, ws.f_resolution)
ws.Copy(ws.f2, ws.f1)
ws.NumericScale(ws.f2, ws.f2, -1.0)
ws.VectorLinSpace(ws.f_backend, ws.f1, ws.f2, ws.f_resolution)
ws.VectorAddScalar(ws.f_backend, ws.f_backend, ws.v0)
ws.VectorCreate("yREFERENCE")
# --- Spectrum for "direct" observation (load switching) ---------------------
ws.AntennaOff()
ws.sensor_responseInit()
ws.sensor_responseBackend()
#
ws.jacobianOff()
#
ws.propmat_clearsky_agenda_checkedCalc()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
ws.WriteXML(ws.output_file_format, ws.y, "TestGbased.y1.xml")
ws.WriteXML(ws.output_file_format, ws.y_f, "TestGbased.f.xml")
ws.ReadXML(ws.yREFERENCE, "TestGbased.y1REFERENCE.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.0001)
# --- Beam switching ---------------------------------------------------------
ws.NumericCreate("za_negative")
ws.Copy(ws.za_negative, ws.za)
ws.NumericScale(ws.za_negative, ws.za_negative, -1.0)
ws.VectorCreate("za_vector")
ws.VectorNLinSpace(ws.za_vector, 2, ws.za_negative, 0.0)
ws.Matrix1ColFromVector(ws.mblock_dlos_grid, ws.za_vector)
ws.sensor_responseInit()
ws.sensor_responseBeamSwitching()
ws.sensor_responseBackend()
ws.sensor_checkedCalc()
ws.yCalc()
#
ws.WriteXML(ws.output_file_format, ws.y, "TestGbased.y2.xml")
ws.ReadXML(ws.yREFERENCE, "TestGbased.y2REFERENCE.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.0001)
# --- Frequency switching -----------------------------------------------------
ws.Copy(ws.f1, ws.f_switch)
ws.NumericScale(ws.f1, ws.f1, -1.0)
ws.AntennaOff()
ws.sensor_responseInit()
ws.sensor_responseBackendFrequencySwitching(
    ws.sensor_response,
    ws.sensor_response_f,
    ws.sensor_response_pol,
    ws.sensor_response_dlos,
    ws.sensor_response_f_grid,
    ws.sensor_response_pol_grid,
    ws.sensor_response_dlos_grid,
    ws.f_backend,
    ws.backend_channel_response,
    ws.sensor_norm,
    ws.f1,
    ws.f_switch,
)
ws.sensor_checkedCalc()
ws.yCalc()
#
ws.WriteXML(ws.output_file_format, ws.y, "TestGbased.y3.xml")
ws.ReadXML(ws.yREFERENCE, "TestGbased.y3REFERENCE.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.0001)
