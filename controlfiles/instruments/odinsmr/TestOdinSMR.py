# DEFINITIONS:  -*-sh-*-
# This is a test doing Odin-SMR simulations.
#
# The calculations are using multiple measuremeent blocks and the control file
# can easily be modified to work for 2D or 3D atmospheres.
# For 1D limb sounding it can be more effecient to do the simulations inside
# a single measurement block. Especially if the number tangent altitudes are
# high and the spacing between tangent altitudes small. See TestOdinSMR_1D.arts
# for an example on such calculations.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Basic settings (already needed in sensor part)
# ---
# This example assumes 1D
ws.AtmosphereSet1D()
# scalar RT
ws.IndexSet(ws.stokes_dim, 1)
# Select frequency band here:
#
ws.execute_controlfile("odinsmr_501.arts")
# INCLUDE "odinsmr_544.arts"
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
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
# Number of tangent altitudes
ws.IndexCreate("n_tan")
ws.IndexSet(ws.n_tan, 4)
# Sensor position, platform altitude set to 600 km
ws.MatrixSetConstant(ws.sensor_pos, ws.n_tan, 1, 600000.0)
# LOS, specified by the corresponding geometrical tangent altitudes
# Tangent altitudes will be equally spaced between 50 and 20 km
ws.VectorCreate("z_tan")
ws.VectorNLinSpace(ws.z_tan, ws.n_tan, 50000.0, 20000.0)
ws.VectorCreate("za")
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.z_tan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
ws.sensor_checkedCalc()
# ---- Calculate spectra and save
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.propmat_clearsky_agenda_checkedCalc()
ws.yCalc()
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "yREFERENCE.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.01)
ws.WriteXML(ws.output_file_format, ws.y)
ws.WriteXML(ws.output_file_format, ws.z_tan)
ws.WriteXML(ws.output_file_format, ws.sensor_response_f_grid)
ws.WriteXML(ws.output_file_format, ws.y_f)
ws.WriteXML(ws.output_file_format, ws.y_pol)
ws.WriteXML(ws.output_file_format, ws.y_pos)
ws.WriteXML(ws.output_file_format, ws.y_los)
