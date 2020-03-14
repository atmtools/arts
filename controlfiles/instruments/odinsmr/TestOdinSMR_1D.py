# DEFINITIONS:  -*-sh-*-
# This is a test doing Odin-SMR simulations. Same test as in
# testOdinSMR, but the calculations are here performed inside a single
# measurement block. This approach can only be used for 1D.
#
# Author: Patrick Eriksson

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Basic settings (already needed in sensor part)
# ---
# This example assumes 1D
ws.AtmosphereSet1D()
# scalar RT
ws.IndexSet(ws.stokes_dim, 1)
# Select frequency band here:
#
ws.execute_controlfile("odinsmr_501_1D.arts")
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
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Possibility to change considered species
# For example
# abs_speciesSet( species=[
#  "H2O,H2O-ForeignContStandardType,H2O-SelfContStandardType",
#  "O3"
# ] )
# ---- Atmospheric scenario --------------------------------------------------
# A pressure grid rougly matching 0 to 80 km in 250 m steps.
# The pressure grid is for the SMR processing not uniform. It is there
# created to be most dense over the actual range of tangent altitudes.
#
ws.VectorNLogSpace(ws.p_grid, 321, 100000.0, 1.0)
# Atmospheric profiles here taken from Fascod
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalc()
# Get ground altitude (z_surface) from z_field
ws.Extract(ws.z_surface, ws.z_field, 0)
# No jacobian calculations
ws.jacobianOff()
# No scattering
ws.cloudboxOff()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# ---- Create absorption table -----------------------------------------------
ws.abs_lines_per_speciesCreateFromLines()
ws.AbsInputFromAtmFields()
ws.abs_speciesSet(abs_species=ws.abs_nls, species=[])
ws.VectorSet(ws.abs_nls_pert, [])
ws.VectorSet(ws.abs_t_pert, [])
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.abs_lookupCalc()
# ---- Sensor position and LOS -----------------------------------------------
# Tangent altitudes are here determined when creating sensor response matrix.
# The tangent altitudes are stored in *vector_2*.
# ---- Calculate and save
ws.propmat_clearsky_agenda_checkedCalc()
ws.yCalc()
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "yREFERENCE_1D.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.01)
ws.WriteXML(ws.output_file_format, ws.y)
ws.WriteXML(ws.output_file_format, ws.z_tan)
ws.WriteXML(ws.output_file_format, ws.sensor_response_f_grid)
ws.WriteXML(ws.output_file_format, ws.y_f)
ws.WriteXML(ws.output_file_format, ws.y_pol)
ws.WriteXML(ws.output_file_format, ws.y_pos)
ws.WriteXML(ws.output_file_format, ws.y_los)
