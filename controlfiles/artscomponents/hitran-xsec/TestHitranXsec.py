# DEFINITIONS:  -*-sh-*-
#
# Simple test for HITRAN cross-section species
#
# Author: Oliver Lemke

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.ReadXML(ws.hitran_xsec_data, "CFC11.xml.gz")
# Agenda for scalar gas absorption calculation
@arts_agenda
def abs_xsec_agenda(ws):
    ws.Ignore(ws.abs_nlte)
    ws.Ignore(ws.abs_vmrs)
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddHitranXsec()


ws.abs_xsec_agenda = abs_xsec_agenda

# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# No jacobian calculation
#
ws.jacobianOff()
# Clearsky = No scattering
#
ws.cloudboxOff()
ws.VectorNLinSpace(ws.f_grid, 1000, 24200000000000.0, 33500000000000.0)
ws.WriteXML("ascii", ws.f_grid)
ws.VectorNLogSpace(ws.p_grid, 20, 100000.0, 1.0)
# Definition of species
# ---
ws.abs_speciesSet(species=["CFC11-HXSEC"])
# Atmospheric scenario
# ---
ws.AtmRawRead(basename="testdata/tropical")
# Weakly reflecting surface
# ---
ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.0)
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
)
# No sensor properties
# ---
ws.sensorOff()
ws.StringSet(ws.iy_unit, "1")
# StringSet( iy_unit, "PlanckBT" )
#########################################################################
# 1D
#########################################################################
# Atmosphere and surface
# ---
ws.AtmosphereSet1D()
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
# Definition of sensor position and LOS
# ---
# MatrixSetConstant( sensor_pos, 3, 1, 600e3 )
# MatrixSet( sensor_los, [ 95; 113; 135] )
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 600000.0)
ws.MatrixSet(ws.sensor_los, np.array([[135.0]]))
# Perform RT calculations
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.timerStart()
ws.yCalc()
ws.timerStop()
ws.Print(ws.timer, 1)
# OK?
# ---
ws.WriteXML(ws.output_file_format, ws.y, "", 0)
ws.WriteXML(ws.output_file_format, ws.f_grid, "", 0)
# WriteXML( "ascii", y, "yREFERENCE.xml" )
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "yREFERENCE.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.01)
# End of Main
