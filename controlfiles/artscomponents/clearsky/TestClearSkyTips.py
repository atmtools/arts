# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of simple ARTS 1D, 2D and 3D clear sky calculations.
#
# Observations from a satellite is treated, with three viewing directions:
#   1: Cold space (ie. above the model atmosphere)
#   2: Limb sounding
#   3: Downward observation.
#
# For the test sequence (ie. "make check") this is a first test on that the
# full chain around yCalc is working. Calculation of optical depth as an
# auxilary variable is included and tested.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Read TIPS partition functions from arts-xml-data
ws.ReadXML(ws.partition_functions, "spectroscopy/PartitionSums/TIPS/tips.xml")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
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
# Read a line file and a matching small frequency grid
# ---
ws.ReadARTSCAT(ws.abs_lines, "abs_lines.xml")
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.VectorNLinSpace(ws.f_grid, 5, 320000000000.0, 322000000000.0)
# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
# ---
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
# Definition of species
# ---
ws.abs_speciesSet(
    species=[
        "H2O-SelfContStandardType, H2O-ForeignContStandardType, H2O",
        "N2-SelfContStandardType",
        "O3",
    ]
)
# Sort the line file according to species
# ---
ws.abs_lines_per_speciesCreateFromLines()
# Atmospheric scenario
# ---
ws.AtmRawRead(basename="testdata/tropical")
# Weakly reflecting surface
# ---
ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.8)
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
)
# No sensor properties
# ---
ws.sensorOff()
# We select here to use Rayleigh-Jean brightness temperatures
# ---
ws.StringSet(ws.iy_unit, "RJBT")
# Extract radiative background and optical depth as auxiliary variables
# ---
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Optical depth", "Radiative background"])
# Create vector container for the optical depth
ws.VectorCreate("odepth")
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
ws.MatrixSetConstant(ws.sensor_pos, 3, 1, 600000.0)
ws.MatrixSet(ws.sensor_los, array([[95.0], [113.0], [135.0]]))
# Perform RT calculations
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
ws.yCalc()
# OK?
# ---
ws.Extract(ws.odepth, ws.y_aux, 0)
# WriteXML( "ascii", y, "yREFERENCE_1D.xml" )
# WriteXML( "ascii", odepth, "y_auxREFERENCE_1D.xml" )
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "yREFERENCE_1D.xml")
# Increased by 10% compared to test using
# builtin partition function coeffs
ws.Compare(ws.y, ws.yREFERENCE, 0.011)
ws.ReadXML(ws.yREFERENCE, "y_auxREFERENCE_1D.xml")
ws.Compare(ws.odepth, ws.yREFERENCE, 0.001)
#########################################################################
# 2D
#########################################################################
# Atmosphere and surface
# ---
ws.AtmosphereSet2D()
ws.VectorLinSpace(ws.lat_grid, -45.0, 45.0, 1.0)
ws.AtmFieldsCalcExpand1D()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
ws.refellipsoidEarth(ws.refellipsoid, "WGS84")
# Definition of sensor position and LOS
# ---
ws.MatrixCreate("zeros")
ws.MatrixSetConstant(ws.zeros, 3, 1, 0.0)
ws.Append(ws.sensor_pos, ws.zeros, "trailing")
# Perform RT calculations
# ---
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.yCalc()
# OK?
# ---
ws.Extract(ws.odepth, ws.y_aux, 0)
# WriteXML( "ascii", y, "yREFERENCE_2D.xml" )
# WriteXML( "ascii", odepth, "y_auxREFERENCE_2D.xml" )
ws.ReadXML(ws.yREFERENCE, "yREFERENCE_2D.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.01)
ws.ReadXML(ws.yREFERENCE, "y_auxREFERENCE_2D.xml")
ws.Compare(ws.odepth, ws.yREFERENCE, 0.001)
#########################################################################
# 3D
#########################################################################
# Atmosphere and surface
# ---
ws.AtmosphereSet3D()
ws.VectorLinSpace(ws.lon_grid, -45.0, 45.0, 1.0)
ws.AtmFieldsCalcExpand1D()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
# Definition of sensor position and LOS
# ---
ws.Append(ws.sensor_pos, ws.zeros, "trailing")
ws.Append(ws.sensor_los, ws.zeros, "trailing")
# Perform RT calculations
# ---
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.yCalc()
# OK?
# ---
ws.Extract(ws.odepth, ws.y_aux, 0)
# WriteXML( "ascii", y, "yREFERENCE_3D.xml" )
# WriteXML( "ascii", odepth, "y_auxREFERENCE_3D.xml" )
ws.ReadXML(ws.yREFERENCE, "yREFERENCE_3D.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.01)
ws.ReadXML(ws.yREFERENCE, "y_auxREFERENCE_3D.xml")
ws.Compare(ws.odepth, ws.yREFERENCE, 0.001)
# End of Main
