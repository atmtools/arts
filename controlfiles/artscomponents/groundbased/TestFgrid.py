# A control file to test polynomial "filling" of spectra

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
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
# ---- Species ---------------------------------------------------------------
ws.abs_speciesSet(species=["O3", "H2O"])
# ---- Atmospheric scenario --------------------------------------------------
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# A pressure grid rougly matching 0 to 80 km in 1 km steps.
ws.VectorNLogSpace(ws.p_grid, 81, 101300.0, 1.0)
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
# ---- Absorption ------------------------------------------------------------
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
ws.ReadARTSCAT(ws.abs_lines, "testdata/ozone_line.xml")
ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", 750000000000.0)
ws.abs_linesSetNormalization(ws.abs_lines, "VVH")
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_cont_descriptionInit()
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
# ---- Create fast f_grid ----------------------------------------------------
ws.NumericCreate("v0")
ws.NumericCreate("fw")
ws.IndexCreate("nlogpart")
ws.NumericCreate("fw_fine")
ws.NumericCreate("df_fine")
# Centre frequency
ws.NumericSet(ws.v0, 110836040000.0)
# One sided width of f_grid
ws.NumericSet(ws.fw, 500000000.0)
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
# --- Checks & Misc settings --------------------------------------------------
ws.IndexSet(ws.stokes_dim, 1)
ws.jacobianOff()
ws.cloudboxOff()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.lbl_checkedCalc()
ws.StringSet(ws.iy_unit, "RJBT")
# -- Calculate spectrum for monochromatic grid and store as y1(x1)
ws.sensorOff()
ws.sensor_checkedCalc()
ws.yCalc()
ws.VectorCreate("x1")
ws.VectorCreate("y1")
ws.Copy(ws.x1, ws.f_grid)
ws.Copy(ws.y1, ws.y)
#
ws.WriteXML(ws.output_file_format, ws.x1)
ws.WriteXML(ws.output_file_format, ws.y1)
# ---- Create reference f_grid -----------------------------------------------
ws.Copy(ws.f1, ws.fw)
ws.NumericScale(ws.f1, ws.f1, -1.0)
ws.VectorLinSpace(ws.f_grid, ws.f1, ws.fw, 30000.0)
ws.VectorAddScalar(ws.f_grid, ws.f_grid, ws.v0)
# -- Calculate spectrum for reference grid and store as ye(x2)
ws.sensorOff()
ws.sensor_checkedCalc()
ws.yCalc()
ws.VectorCreate("x2")
ws.VectorCreate("ye")
ws.Copy(ws.x2, ws.f_grid)
ws.Copy(ws.ye, ws.y)
#
ws.WriteXML(ws.output_file_format, ws.x2)
ws.WriteXML(ws.output_file_format, ws.ye)
# ---- Create spectra with polynomial filling ---------------------------------
ws.Copy(ws.f_grid, ws.x1)
ws.FlagOn(ws.sensor_norm)
ws.AntennaOff()
ws.sensor_responseInit()
ws.sensor_responseFillFgrid(
    ws.sensor_response,
    ws.sensor_response_f,
    ws.sensor_response_pol,
    ws.sensor_response_dlos,
    ws.sensor_response_f_grid,
    ws.sensor_response_pol_grid,
    ws.sensor_response_dlos_grid,
    3,
    2,
)
ws.sensor_checkedCalc()
ws.yCalc()
ws.WriteXML(ws.output_file_format, ws.sensor_response_f_grid, "TestFgrid.x3.xml")
ws.WriteXML(ws.output_file_format, ws.y, "TestFgrid.y3.xml")
ws.sensor_responseInit()
ws.sensor_responseFillFgrid(
    ws.sensor_response,
    ws.sensor_response_f,
    ws.sensor_response_pol,
    ws.sensor_response_dlos,
    ws.sensor_response_f_grid,
    ws.sensor_response_pol_grid,
    ws.sensor_response_dlos_grid,
    5,
    4,
)
ws.sensor_checkedCalc()
ws.yCalc()
ws.WriteXML(ws.output_file_format, ws.sensor_response_f_grid, "TestFgrid.x5.xml")
ws.WriteXML(ws.output_file_format, ws.y, "TestFgrid.y5.xml")
