# DEFINITIONS:  -*-sh-*-
#
# A simple demonstration and test of iyIndependentBeamApproximation
#
# 2019-09-17, Patrick Eriksson

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
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Basic settings
#
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, array([3.10e10, 8.90e10, 1.60e11, 1.84e11]))
ws.StringSet(ws.iy_unit, "PlanckBT")
# no jacobian calculation
ws.jacobianOff()
# no scattering
ws.cloudboxOff()
# no sensor
ws.sensorOff()
# Definition of species
#
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"])
# No transitions needed
#
ws.abs_lines_per_speciesSetEmpty()
# some checks
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.lbl_checkedCalc()
# Surface = Ocean by FASTEM
#
ws.VectorCreate("trv")
ws.nelemGet(v=ws.f_grid)
ws.VectorSetConstant(ws.trv, ws.nelem, 0.9)


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.surfaceFastem(
        wind_speed=5.0, wind_direction=45.0, fastem_version=5, transmittance=ws.trv
    )


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Read a 1D atmospheric case
#
ws.AtmRawRead(basename="testdata/tropical")
# 1D
#
ws.AtmosphereSet1D()
ws.VectorNLogSpace(ws.p_grid, 81, 105000.0, 10000.0)
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)
#
ws.AtmFieldsCalc()
#
ws.MatrixSet(ws.sensor_pos, array([[800000.0]]))
ws.MatrixSet(ws.sensor_los, array([[130.0]]))
#
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
ws.VectorCreate("y1d")
ws.Copy(ws.y1d, ws.y)
# Expand to 2D and run IBA
#
ws.AtmosphereSet2D()
ws.IndexCreate("nlat")
ws.IndexSet(ws.nlat, 3)
ws.VectorNLinSpace(ws.lat_grid, ws.nlat, -5.0, 45.0)
ws.MatrixSetConstant(ws.z_surface, ws.nlat, 1, 0.0)
#
ws.Copy(ws.lat_true, ws.lat_grid)
ws.VectorSetConstant(ws.lon_true, ws.nlat, 0.0)
#
ws.AtmFieldsCalcExpand1D()
#
ws.MatrixSet(ws.sensor_pos, array([[800000.0, 0.0]]))
ws.MatrixSet(ws.sensor_los, array([[130.0]]))
#
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
#
@arts_agenda
def iy_main_agenda(ws):
    ws.iyIndependentBeamApproximation()


ws.iy_main_agenda = iy_main_agenda


@arts_agenda
def iy_independent_beam_approx_agenda(ws):
    ws.Ignore(ws.lat_grid)
    ws.Ignore(ws.lon_grid)
    ws.Ignore(ws.lat_true)
    ws.Ignore(ws.lon_true)
    ws.Ignore(ws.z_surface)
    ws.Ignore(ws.z_field)
    ws.Ignore(ws.cloudbox_limits)
    ws.Ignore(ws.pnd_field)
    #
    ws.ppathCalc()
    ws.iyEmissionStandard()


ws.iy_independent_beam_approx_agenda = iy_independent_beam_approx_agenda

#
ws.yCalc()
# Comapre
ws.Compare(ws.y, ws.y1d, 0.001)
