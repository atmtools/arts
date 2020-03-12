# DEFINITIONS:  -*-sh-*-
#
# This file tests temperature Jacobian calculations. With some focus on
# inclusion of HSE, and mainly for 1D atmosphere and Stokes dim 1. At the end it
# is also tested that a 3D case with latitude and longitude retrieval grids of
# length 1 gives the same Jacobian as 1D.
#
# Three frequencies are done, one with high suface sensitivity and two around
# 118 GHz (with Jacobian peaking around 30 and 70 km.
#
# 2018-11-18, Patrick Eriksson

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
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# Geometrical path calculation (i.e., refraction neglected)
#
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Standard RT agendas
#
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# Definition of species
#
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Atmosphere
#
ws.AtmosphereSet1D()
ws.VectorNLogSpace(ws.p_grid, 161, 101300.0, 1.0)
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalc()
# Surface
#
# Don't use interpolation of t_field to set surface temperature. That will
# cause a difference between analytical and perturbation Jacobian
#
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
ws.VectorSet(ws.surface_scalar_reflectivity, np.array([0.4]))
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
)
# Frequencies and Stokes dim.
#
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, np.array([3.5000e10, 1.1875e11, 1.1880e11]))
# Sensor pos and los
#
ws.MatrixSet(ws.sensor_pos, np.array([[820000.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[140.0]]))
# Define analytical Jacobian
#
ws.jacobianInit()
ws.jacobianAddTemperature(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, hse="off")
ws.jacobianClose()
# Deactive parts not used
#
ws.cloudboxOff()
ws.sensorOff()
# Checks
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
# HSE
#
ws.VectorSet(ws.lat_true, np.array([0.0]))
ws.VectorSet(ws.lon_true, np.array([0.0]))
#
ws.Extract(ws.p_hse, ws.p_grid, 0)
ws.NumericSet(ws.z_hse_accuracy, 0.5)
ws.z_fieldFromHSE()
# Run RT calcs
#
ws.StringSet(ws.iy_unit, "RJBT")
#
ws.yCalc()
# Check y against reference
#
ws.VectorCreate("yref")
ws.ReadXML(ws.yref, "yREF1.xml")
ws.Compare(
    ws.y, ws.yref, 0.0001, "Calculated *y* does not agree with saved reference values."
)
# Copy Jacobian
#
ws.MatrixCreate("jcopy")
ws.Copy(ws.jcopy, ws.jacobian)
# Save
#
# output_file_formatSetAscii
# WriteXML( output_file_format, y, "yREF1.xml" )
# WriteXML( output_file_format, f_grid, "f.xml" )
# WriteXML( output_file_format, z_field, "z.xml" )
# WriteXML( output_file_format, jacobian, "Ja_off.xml" )
# Re-do, HSE off, perturbation
#
ws.NumericCreate("dt")
ws.NumericSet(ws.dt, 0.1)
ws.IndexNumberOfAtmosphericPoints(n=ws.ybatch_n)
#
ws.jacobianOff()
#
@arts_agenda
def ybatch_calc_agenda(ws):
    ws.AtmFieldPerturbAtmGrids(
        perturbed_field=ws.t_field,
        original_field=ws.t_field,
        pert_index=ws.ybatch_index,
        pert_size=ws.dt,
    )
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

#
ws.ybatchCalc(ybatch_start=0)
ws.jacobianFromYbatch(pert_size=ws.dt)
#
# WriteXML( output_file_format, jacobian, "Jp_off.xml" )
#
ws.Compare(
    ws.jacobian,
    ws.jcopy,
    0.0001,
    "Analytical and perturbation Jacobian disagree with HSE=off",
)
# Re-do, HSE on, analytical
#
ws.jacobianInit()
ws.jacobianAddTemperature(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, hse="on")
ws.jacobianClose()
#
ws.yCalc()
#
# WriteXML( output_file_format, jacobian, "Ja_on.xml" )
#
ws.Copy(ws.jcopy, ws.jacobian)
# Re-do, HSE on, perturbation
#
ws.jacobianOff()
#
@arts_agenda
def ybatch_calc_agenda(ws):
    ws.AtmFieldPerturbAtmGrids(
        perturbed_field=ws.t_field,
        original_field=ws.t_field,
        pert_index=ws.ybatch_index,
        pert_size=ws.dt,
    )
    ws.z_fieldFromHSE()
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

#
ws.ybatchCalc(ybatch_start=0)
ws.jacobianFromYbatch(pert_size=ws.dt)
#
# WriteXML( output_file_format, jacobian, "Jp_on.xml" )
#
ws.Compare(
    ws.jcopy,
    ws.jacobian,
    0.0001,
    "Analytical and perturbation Jacobian disagree with HSE=on",
)
# Move to a 3D view and redo analytical with HSE=on
#
ws.AtmosphereSet3D()
ws.VectorSet(ws.lat_grid, np.array([-10.0, 10.0]))
ws.Copy(ws.lon_grid, ws.lat_grid)
ws.AtmFieldsCalcExpand1D()
#
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
#
ws.MatrixSet(ws.sensor_pos, np.array([[820000.0, 0.0, 0.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[140.0, 20.0]]))
#
ws.VectorCreate("lat0")
ws.VectorCreate("lon0")
ws.VectorSet(ws.lat0, np.array([0.0]))
ws.VectorSet(ws.lon0, np.array([0.0]))
ws.jacobianInit()
ws.jacobianAddTemperature(g1=ws.p_grid, g2=ws.lat0, g3=ws.lon0, hse="on")
ws.jacobianClose()
#
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
# WriteXML( output_file_format, jacobian, "J3d.xml" )
#
# There is a bit of "noise" in 118.75 GHz Jacobian (for unknown reason)
# and we must allow a bit higher deviation
ws.Compare(
    ws.jacobian, ws.jcopy, 0.002, "Jacobians for 1D and matching 3D view do not agree."
)
