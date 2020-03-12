# DEFINITIONS:  -*-sh-*-
#
# This file tests temperature Jacobian calculations for fully polarised
# simulations. 3D calculations applied, to allow Zeeman.
#
# Two frequencies are done, one with high suface sensitivity and one at 118 GHz
# with a significant Zeeman signature for all Stokes elements. Surface is
# modelled by FASTEM.
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
# on-the-fly absorption, with Zeeman
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly_ZeemanPreCalc)
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
# Frequencies and Stokes dim.
#
ws.IndexSet(ws.stokes_dim, 4)
ws.VectorSet(ws.f_grid, np.array([3.50000e10, 1.18751e11]))
# Definition of species
#
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-Z-66", "H2O-PWR98"])
# Line data
#
ws.ReadARTSCAT(filename="line_118ghz.xml", localquantumnumbers="J")
ws.abs_lines_per_speciesCreateFromLines()
ws.Wigner6Init(ws.wigner_initialized, 40000, 100)
# Atmosphere
#
ws.AtmosphereSet3D()
ws.VectorNLogSpace(ws.p_grid, 201, 101300.0, 0.05)
ws.VectorSet(ws.lat_grid, np.array([-10.0, 10.0]))
ws.Copy(ws.lon_grid, ws.lat_grid)
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalcExpand1D()
# Magnetic field
#
# Craete a synthetic field, where all componenets have the same value
#
ws.Copy(ws.mag_u_field, ws.t_field)
ws.Tensor3Scale(ws.mag_u_field, ws.mag_u_field, 0.0)
ws.Tensor3AddScalar(ws.mag_u_field, ws.mag_u_field, 2.5e-05)
ws.Copy(ws.mag_v_field, ws.mag_u_field)
ws.Copy(ws.mag_w_field, ws.mag_u_field)
# Surface
#
# Don't use interpolation of t_field to set surface temperature. That will
# cause a difference between analytical and perturbation Jacobian
#
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
ws.VectorCreate("trv")
ws.nelemGet(v=ws.f_grid)
ws.VectorSetConstant(ws.trv, ws.nelem, 0.9)


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceFastem(
        wind_speed=5.0, wind_direction=45.0, fastem_version=5, transmittance=ws.trv
    )


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Sensor pos and los
#
ws.MatrixSet(ws.sensor_pos, np.array([[820000.0, 0.0, 0.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[140.0, 45.0]]))
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
ws.ReadXML(ws.yref, "yREF4.xml")
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
# WriteXML( output_file_format, f_grid, "f.xml" )
# WriteXML( output_file_format, z_field, "z.xml" )
# WriteXML( output_file_format, y, "yREF4.xml" )
# WriteXML( output_file_format, jacobian, "Ja.xml" )
# Re-do by external perturbations
#
ws.NumericCreate("dt")
ws.NumericSet(ws.dt, 0.1)
#
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
# WriteXML( output_file_format, jacobian, "Jp.xml" )
# To make test more stringent, apply different tolerances for different Stokes
#
ws.VectorCreate("v1")
ws.VectorCreate("v2")
# Compare frequency 1, Stokes element 1
#
ws.VectorExtractFromMatrix(ws.v1, ws.jcopy, 0, "row")
ws.VectorExtractFromMatrix(ws.v2, ws.jacobian, 0, "row")
ws.Compare(
    ws.v1, ws.v2, 0.0001, "Disagreement for first frequency and first Stokes element"
)
# Compare frequency 1, Stokes element 2
#
ws.VectorExtractFromMatrix(ws.v1, ws.jcopy, 1, "row")
ws.VectorExtractFromMatrix(ws.v2, ws.jacobian, 1, "row")
ws.Compare(
    ws.v1, ws.v2, 1e-05, "Disagreement for first frequency and second Stokes element"
)
# Compare frequency 1, Stokes element 3
#
ws.VectorExtractFromMatrix(ws.v1, ws.jcopy, 2, "row")
ws.VectorExtractFromMatrix(ws.v2, ws.jacobian, 2, "row")
ws.Compare(
    ws.v1, ws.v2, 1e-05, "Disagreement for first frequency and third Stokes element"
)
# Compare frequency 1, Stokes element 4
#
ws.VectorExtractFromMatrix(ws.v1, ws.jcopy, 3, "row")
ws.VectorExtractFromMatrix(ws.v2, ws.jacobian, 3, "row")
ws.Compare(
    ws.v1, ws.v2, 1e-05, "Disagreement for first frequency and fourth Stokes element"
)
