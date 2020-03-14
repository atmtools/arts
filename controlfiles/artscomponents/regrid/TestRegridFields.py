# DEFINITIONS:  -*-sh-*-
#
# Testing AtmFieldRefinePgrid.
#
# 1D clear sky with an initial and a refined p_grid. Based on
# artscomponents/clearky/TestClearSky.arts.
#
# Author: Jana Mendrok

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
ws.ReadARTSCAT(ws.abs_lines, "artscomponents/clearsky/abs_lines.xml")
ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", 750000000000.0)
ws.abs_linesSetNormalization(ws.abs_lines, "VVH")
ws.VectorNLinSpace(ws.f_grid, 5, 320000000000.0, 322000000000.0)
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
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Optical depth"])
# Definition of sensor position and LOS
# ---
ws.MatrixSetConstant(ws.sensor_pos, 3, 1, 600000.0)
ws.MatrixSet(ws.sensor_los, np.array([[95.0], [113.0], [135.0]]))
#########################################################################
# initial vertical spacing
#########################################################################
# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
# ---
ws.VectorCreate("p_init")
ws.VectorNLogSpace(ws.p_init, 41, 100000.0, 1.0)
ws.Copy(ws.p_grid, ws.p_init)
# Atmosphere and surface
# ---
ws.AtmosphereSet1D()
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
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
# Create containers for comparison data
ws.VectorCreate("odepth_init")
ws.VectorCreate("y_init")
ws.Extract(ws.odepth_init, ws.y_aux, 0)
ws.Copy(ws.y_init, ws.y)
# WriteXML( in=y_init )
# WriteXML( in=odepth_init )
#########################################################################
# refined vertical spacing I: all clearsky fields in one WSM
#########################################################################
ws.AtmFieldsRefinePgrid(p_step=0.01)
ws.VectorCreate("p_ref")
ws.Tensor3Create("z_ref")
ws.Tensor3Create("t_ref")
ws.Tensor4Create("vmr_ref")
ws.Copy(ws.p_ref, ws.p_grid)
# WriteXML( in=p_ref )
ws.Copy(ws.z_ref, ws.z_field)
# WriteXML( in=z_ref )
ws.Copy(ws.t_ref, ws.t_field)
ws.Copy(ws.vmr_ref, ws.vmr_field)
# Perform RT calculations
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.yCalc()
# Create containers for comparison data
ws.VectorCreate("odepth_refine")
ws.VectorCreate("y_refine")
ws.Extract(ws.odepth_refine, ws.y_aux, 0)
ws.Copy(ws.y_refine, ws.y)
# WriteXML( in=y_refine )
# WriteXML( in=odepth_refine )
# OK?
# ---
ws.VectorCreate("odepth_refineREFERENCE")
ws.VectorCreate("y_refineREFERENCE")
ws.ReadXML(out=ws.y_refineREFERENCE)
ws.ReadXML(out=ws.odepth_refineREFERENCE)
ws.Compare(ws.y_refine, ws.y_refineREFERENCE, 0.01)
ws.Compare(ws.odepth_refine, ws.odepth_refineREFERENCE, 0.01)
ws.Compare(ws.y_refine, ws.y_init, 0.5)
ws.Compare(ws.odepth_refine, ws.odepth_init, 2.0)
#########################################################################
# refined vertical spacing II: all fields separately
#########################################################################
# need to recreate the original fields
ws.Copy(ws.p_grid, ws.p_init)
ws.AtmFieldsCalc()
ws.VectorCreate("p_grid_old")
ws.Copy(ws.p_grid_old, ws.p_grid)
ws.p_gridRefine(p_grid_old=ws.p_grid_old, p_step=0.01)
# WriteXML( in=p_grid )
ws.Compare(ws.p_grid, ws.p_ref, 1e-10)
ws.AtmFieldPRegrid(ws.z_field, ws.z_field, ws.p_grid, ws.p_grid_old, 1)
# WriteXML( in=z_field )
ws.Compare(ws.z_field, ws.z_ref, 1e-06)
ws.AtmFieldPRegrid(ws.t_field, ws.t_field, ws.p_grid, ws.p_grid_old, 1)
ws.Compare(ws.t_field, ws.t_ref, 1e-06)
ws.AtmFieldPRegrid(ws.vmr_field, ws.vmr_field, ws.p_grid, ws.p_grid_old, 1)
ws.Compare(ws.vmr_field, ws.vmr_ref, 1e-12)
# Perform RT calculations
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.yCalc()
# Create containers for comparison data
ws.VectorCreate("odepth_refineII")
ws.VectorCreate("y_refineII")
ws.Extract(ws.odepth_refineII, ws.y_aux, 0)
ws.Copy(ws.y_refineII, ws.y)
# WriteXML( in=y_refineII )
# WriteXML( in=odepth_refineII )
ws.Compare(ws.y_refineII, ws.y_refine, 1e-06)
ws.Compare(ws.odepth_refineII, ws.odepth_refine, 1e-06)
# End of Main
