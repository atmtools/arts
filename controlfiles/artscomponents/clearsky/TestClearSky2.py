# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of a simple ARTS 1D clear sky calculations.
#
# Observations from a satellite is treated, with a single observation angle.
# RT is done by iyCalc and ySimpleSpectrometer is used.
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
ws.ReadARTSCAT(abs_lines=ws.abs_lines, filename="abs_lines.xml")
ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", 750000000000.0)
ws.abs_linesSetNormalization(ws.abs_lines, "VVH")
ws.VectorNLinSpace(ws.f_grid, 4001, 321000000000.0, 329000000000.0)
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
# We select here to use Rayleigh-Jean brightness temperatures
# ---
ws.StringSet(ws.iy_unit, "RJBT")
# Atmosphere and surface
# ---
ws.AtmosphereSet1D()
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
# Definition of sensor position and LOS
# ---
ws.VectorSet(ws.rte_pos, array([600000.0]))
ws.VectorSet(ws.rte_los, array([113.4]))
ws.VectorSet(ws.rte_pos2, array([], dtype=float64))
# A dummy value
# Perform RT calculations
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.lbl_checkedCalc()
#
ws.iyCalc()
ws.WriteXML("ascii", ws.f_grid, "f_grid.xml")
ws.WriteXML("ascii", ws.iy, "iy.xml")
ws.ySimpleSpectrometer(ws.y, ws.y_f, ws.iy, ws.stokes_dim, ws.f_grid, 200000000.0)
ws.WriteXML("ascii", ws.y_f, "y_f.xml")
ws.WriteXML("ascii", ws.y, "y.xml")
