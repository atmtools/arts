# DEFINITIONS:  -*-sh-*-
#
# Simple simulations of ground-based measurements of ozone at 110.8 GHz,
# to test impact of winds.
#
# Author: Patrick Eriksson

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
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
# Number of Stokes components to be computed
ws.IndexSet(ws.stokes_dim, 1)
# ---- f_grid ----------------------------------------------------------------
ws.NumericCreate("v0")
ws.NumericSet(ws.v0, 110836040000.0)
ws.VectorLinSpace(ws.f_grid, -5000000.0, 5000000.0, 50000.0)
ws.VectorAddScalar(ws.f_grid, ws.f_grid, ws.v0)
# ---- Species ---------------------------------------------------------------
ws.abs_speciesSet(species=["O3", "H2O"])
# ---- Atmospheric scenario --------------------------------------------------
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# A pressure grid rougly matching 0 to 88 km in about 500 m steps.
ws.IndexCreate("np")
ws.IndexSet(ws.np, 180)
ws.VectorNLogSpace(ws.p_grid, ws.np, 101300.0, 0.5)
ws.AtmRawRead(basename="testdata/tropical")
# All settings here:
#
# ---- Select atmosphere_dim, LOS angles  and winds --------------------------
ws.VectorSet(ws.lat_grid, np.array([-10.0, 10.0]))
ws.VectorSet(ws.lon_grid, np.array([-10.0, 10.0]))
# AtmosphereSet1D
# AtmosphereSet2D
ws.AtmosphereSet3D()
ws.MatrixSet(ws.sensor_los, np.array([[30.0, 50.0]]))
# AtmFieldsCalc
ws.AtmFieldsCalcExpand1D()
ws.nrowsGet(ws.nrows, ws.t_field)
ws.ncolsGet(ws.ncols, ws.t_field)
ws.Tensor3SetConstant(ws.wind_u_field, ws.np, ws.nrows, ws.ncols, 50.0)
ws.Tensor3SetConstant(ws.wind_v_field, ws.np, ws.nrows, ws.ncols, 100.0)
ws.Tensor3SetConstant(ws.wind_w_field, ws.np, ws.nrows, ws.ncols, 2.0)
# ---- Absorption ------------------------------------------------------------
ws.ReadARTSCAT(ws.abs_lines, "testdata/ozone_line.xml")
ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", 750000000000.0)
ws.abs_linesSetNormalization(ws.abs_lines, "VVH")
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_cont_descriptionInit()
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# We have to also set the absorption lookup table interpolation order, even though
# we are not using the lookup table here. atmfields_checkedCalc will otherwise throw an error.
ws.IndexSet(ws.abs_f_interp_order, 1)
# ---- The surface -----------------------------------------------------
ws.MatrixSetConstant(ws.z_surface, ws.nrows, ws.ncols, 0.0)
# ---- Observation position ---------------------------------------------------
ws.MatrixSetConstant(ws.sensor_pos, 1, ws.atmosphere_dim, 0.0)
# ---- Final stuff -----------------------------------------------------------
ws.sensorOff()
ws.jacobianOff()
ws.cloudboxOff()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
ws.StringSet(ws.iy_unit, "RJBT")
# ---- Calculate -------------------------------------------------------------
ws.yCalc()
# WriteXML( "ascii", f_grid, "f.xml" )
# WriteXML( "ascii", y, "yREFERENCE.xml" )
# Expected results
#
ws.VectorCreate("yREFERENCE")
#
ws.ReadXML(ws.yREFERENCE, "yREFERENCE.xml")
#
ws.Compare(ws.y, ws.yREFERENCE, 0.0001)
# ---- Without winds ----------------------------------------------------
# Tensor3SetConstant( wind_u_field,np, nrows, ncols, 0 )
# Tensor3SetConstant( wind_v_field,np, nrows, ncols, 0 )
# Tensor3SetConstant( wind_w_field,np, nrows, ncols, 0 )
# yCalc
# WriteXML( output_file_format, y, "TestWinds.y0.xml" )
