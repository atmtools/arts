# DEFINITIONS:  -*-sh-*-
#
# A simple, basic test of 1D and 2D antenna patterns. The setting of *y_geo* is
# also tested inside this file.
#
# 2018-12-12, Patrick Eriksson

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

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
# geo-positioning
@arts_agenda
def geo_pos_agenda(ws):
    ws.geo_posEndOfPpath()


ws.geo_pos_agenda = geo_pos_agenda

# Basic settings
#
ws.AtmosphereSet3D()
ws.IndexSet(ws.stokes_dim, 2)
ws.VectorSet(ws.f_grid, np.array([1.8e10, 3.1e10]))
ws.StringSet(ws.iy_unit, "PlanckBT")
# no jacobian calculation
ws.jacobianOff()
# no scattering
ws.cloudboxOff()
# lat and lon true can be left empty for 3D
ws.VectorSet(ws.lat_true, [])
ws.VectorSet(ws.lon_true, [])
# Definition of species
#
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"])
# No transitions needed
#
ws.abs_lines_per_speciesSetEmpty()
# Atmospheric grids
#
ws.VectorNLogSpace(ws.p_grid, 41, 105000.0, 10000.0)
ws.VectorNLinSpace(ws.lat_grid, 3, -30.0, 30.0)
ws.VectorNLinSpace(ws.lon_grid, 3, -30.0, 30.0)
# Read a 1D atmospheric case and expand to *atmosphere_dim*
#
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalcExpand1D()
# Surface = Ocean by FASTEM
#
ws.IndexCreate("nlat")
ws.IndexCreate("nlon")
ws.nelemGet(ws.nlat, ws.lat_grid)
ws.nelemGet(ws.nlon, ws.lon_grid)
#
ws.MatrixSetConstant(ws.z_surface, ws.nlat, ws.nlon, 0.0)
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

# Check data and settings (beside sensor)
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.lbl_checkedCalc()
# Sensor pos/los
#
ws.MatrixSet(ws.sensor_pos, np.array([[800000.0, 0.0, 0.0], [800000.0, 0.0, 0.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[130.0, 20.0], [130.0, 15.0]]))
# First do without antenna to get reference for geo_pos
#
ws.sensorOff()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
ws.MatrixCreate("geo_ref")
ws.Copy(ws.geo_ref, ws.y_geo)
# Define 1D antenna
#
ws.NumericCreate("xwidth")
ws.NumericSet(ws.xwidth, 3.0)
ws.VectorCreate("dza")
ws.VectorCreate("zeros")
ws.IndexCreate("ndlos")
#
ws.antenna_responseVaryingGaussian(
    ws.antenna_response, 1.5, ws.xwidth, 0.1, 5, 10000000000.0, 40000000000.0, 0
)
ws.IndexSet(ws.sensor_norm, 1)
ws.IndexSet(ws.antenna_dim, 1)
ws.MatrixSet(ws.antenna_dlos, np.array([[0.0]]))
#
ws.IndexSet(ws.ndlos, 21)
ws.VectorNLinSpace(ws.dza, ws.ndlos, -2.0, 2.0)
ws.VectorSetConstant(ws.zeros, ws.ndlos, 0.0)
ws.Matrix1ColFromVector(ws.mblock_dlos_grid, ws.dza)
# Calculate
#
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
ws.VectorCreate("y_ref")
ws.Copy(ws.y_ref, ws.y)
#
ws.Compare(ws.geo_ref, ws.y_geo, 1e-05)
# Re-do with 2D antenna (but still "1D" mblock_dlos_grid)
#
ws.antenna_responseVaryingGaussian(
    ws.antenna_response, 1.5, ws.xwidth, 0.1, 5, 10000000000.0, 40000000000.0, 1
)
ws.WriteXML("ascii", ws.antenna_response, "R.xml")
ws.IndexSet(ws.antenna_dim, 2)
# Calculate
#
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
ws.Compare(ws.geo_ref, ws.y_geo, 1e-05)
# Results for *y* should not be identical, but fairly close
ws.Compare(ws.y, ws.y_ref, 0.001)
# Repeat with circular mblock_dlos_grid
#
ws.mblock_dlos_gridUniformCircular(spacing=0.2, width=2.0, centre=1)
ws.WriteXML("ascii", ws.mblock_dlos_grid, "mdg.xml")
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_checkedCalc()
ws.WriteXML("ascii", ws.sensor_response, "H.xml")
#
ws.yCalc()
#
ws.Compare(ws.geo_ref, ws.y_geo, 1e-05)
# Some deviation expected here for *y*
ws.Compare(ws.y, ws.y_ref, 0.01)
# reset y_ref
ws.Copy(ws.y_ref, ws.y)
# Repeat with rectangular mblock_dlos_grid
#
ws.mblock_dlos_gridUniformRectangular(spacing=0.2, za_width=2.0, aa_width=2.0, centre=1)
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_checkedCalc()
#
ws.yCalc()
#
ws.Compare(ws.geo_ref, ws.y_geo, 1e-05)
# *y*should basically be identical
ws.Compare(ws.y, ws.y_ref, 1e-06)
