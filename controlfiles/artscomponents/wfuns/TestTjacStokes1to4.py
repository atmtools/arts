# DEFINITIONS:  -*-sh-*-
#
# This file tests consistency between Stokes 1, 2, 3 and 4, involving a
# polaristaion signature, across the Stokes elements, created by the surface
#
# One frequency is done. The surface is modelled by FASTEM.
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
# Frequencies and Stokes dim.
#
ws.IndexSet(ws.stokes_dim, 4)
ws.VectorSet(ws.f_grid, np.array([3.5e10]))
# Definition of species
#
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Atmosphere
#
ws.AtmosphereSet1D()
ws.VectorNLogSpace(ws.p_grid, 101, 101300.0, 10000.0)
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
ws.MatrixSet(ws.sensor_pos, np.array([[820000.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[140.0]]))
# Define analytical Jacobian
#
ws.jacobianInit()
ws.jacobianAddTemperature(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, hse="on")
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
# Copy results
#
ws.NumericCreate("i4")
ws.NumericCreate("q4")
ws.NumericCreate("u4")
ws.Extract(ws.i4, ws.y, 0)
ws.Extract(ws.q4, ws.y, 1)
ws.Extract(ws.u4, ws.y, 2)
#
ws.VectorCreate("ji4")
ws.VectorCreate("jq4")
ws.VectorCreate("ju4")
ws.VectorExtractFromMatrix(ws.ji4, ws.jacobian, 0, "row")
ws.VectorExtractFromMatrix(ws.jq4, ws.jacobian, 1, "row")
ws.VectorExtractFromMatrix(ws.ju4, ws.jacobian, 2, "row")
# Required agreement
#
ws.NumericCreate("delta")
ws.NumericSet(ws.delta, 1e-09)
# Stokes = 3
#
ws.IndexSet(ws.stokes_dim, 3)
ws.sensorOff()
ws.yCalc()
#
ws.NumericCreate("i")
ws.NumericCreate("q")
ws.NumericCreate("u")
ws.Extract(ws.i, ws.y, 0)
ws.Extract(ws.q, ws.y, 1)
ws.Extract(ws.u, ws.y, 2)
ws.Compare(ws.i, ws.i4, ws.delta, "Disagreement for I beteen Stokes 4 and 3.")
ws.Compare(ws.q, ws.q4, ws.delta, "Disagreement for Q beteen Stokes 4 and 3.")
ws.Compare(ws.u, ws.u4, ws.delta, "Disagreement for U beteen Stokes 4 and 3.")
#
ws.VectorCreate("ji")
ws.VectorCreate("jq")
ws.VectorCreate("ju")
ws.VectorExtractFromMatrix(ws.ji, ws.jacobian, 0, "row")
ws.VectorExtractFromMatrix(ws.jq, ws.jacobian, 1, "row")
ws.VectorExtractFromMatrix(ws.ju, ws.jacobian, 2, "row")
ws.Compare(
    ws.ji,
    ws.ji4,
    ws.delta,
    "Disagreement for I-part of Jacobian beteen Stokes 4 and 3.",
)
ws.Compare(
    ws.jq,
    ws.jq4,
    ws.delta,
    "Disagreement for Q-part of Jacobian beteen Stokes 4 and 3.",
)
ws.Compare(
    ws.ju,
    ws.ju4,
    ws.delta,
    "Disagreement for U-part of Jacobian beteen Stokes 4 and 3.",
)
# Stokes = 2
#
ws.IndexSet(ws.stokes_dim, 2)
ws.sensorOff()
ws.yCalc()
#
ws.Extract(ws.i, ws.y, 0)
ws.Extract(ws.q, ws.y, 1)
ws.Compare(ws.i, ws.i4, ws.delta, "Disagreement for I beteen Stokes 4 and 2.")
ws.Compare(ws.q, ws.q4, ws.delta, "Disagreement for Q beteen Stokes 4 and 2.")
#
ws.VectorExtractFromMatrix(ws.ji, ws.jacobian, 0, "row")
ws.VectorExtractFromMatrix(ws.jq, ws.jacobian, 1, "row")
ws.Compare(
    ws.ji,
    ws.ji4,
    ws.delta,
    "Disagreement for I-part of Jacobian beteen Stokes 4 and 2.",
)
ws.Compare(
    ws.jq,
    ws.jq4,
    ws.delta,
    "Disagreement for Q-part of Jacobian beteen Stokes 4 and 2.",
)
# Stokes = 1
#
ws.IndexSet(ws.stokes_dim, 1)
ws.sensorOff()
ws.yCalc()
#
ws.Extract(ws.i, ws.y, 0)
ws.Compare(ws.i, ws.i4, ws.delta, "Disagreement for I beteen Stokes 4 and 2.")
#
ws.VectorExtractFromMatrix(ws.ji, ws.jacobian, 0, "row")
ws.Compare(
    ws.ji,
    ws.ji4,
    ws.delta,
    "Disagreement for I-part of Jacobian beteen Stokes 4 and 1.",
)
