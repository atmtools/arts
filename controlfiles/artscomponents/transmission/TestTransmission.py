# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of a transmission calculation in a
# refractive 3D atmosphere. The case treats ground-based sensor, observing
# at a high/low zenith/elevation angle.
#
# The control file performs also a test of iyLoopFrequencies
#
# 2012-04-02, Patrick Eriksson

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
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 4)
# Reference ellipsoid
#
ws.refellipsoidEarth(ws.refellipsoid, "WGS84")
# Frequency grid
#
ws.VectorSet(ws.f_grid, np.array([1.0e10, 2.0e10]))
# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
#
ws.VectorNLogSpace(ws.p_grid, 41, 101300.0, 1.0)
# Atmospheric dimensionality and lat/lon grids
#
ws.VectorNLinSpace(ws.lat_grid, 11, 5.0, 13.0)
ws.VectorNLinSpace(ws.lon_grid, 11, -14.0, -10.0)
ws.AtmosphereSet3D()
# Definition of species
#
ws.abs_speciesSet(species=["H2O-PWR98", "N2-SelfContStandardType", "O2-PWR93"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Atmospheric profiles
#
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalcExpand1D()
# Get ground altitude (z_surface) from z_field
ws.MatrixSetConstant(ws.z_surface, 11, 11, 0.0)
# No jacobian calculations
#
ws.jacobianOff()
# No scattering
#
ws.cloudboxOff()
# Check model atmosphere
#
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
# Propagation path agendas and variables
#
# refracted path
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesEarth)
ws.NumericSet(ws.ppath_lmax, 2000.0)
ws.NumericSet(ws.ppath_lraytrace, 500.0)
# Postion and line-of-sight of sensor
#
ws.VectorSet(ws.rte_pos, np.array([0.0, 5.1, -13.82]))
ws.VectorSet(ws.rte_los, np.array([80.0, 24.0]))
ws.VectorSet(ws.rte_pos2, [])
# No transmitter position defined
# Radiative transfer agendas
#
ws.Copy(ws.iy_transmitter_agenda, ws.iy_transmitter_agenda__UnitUnpolIntensity)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Transmission)
# Auxiliary variables
#
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Radiative background", "Optical depth"])
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.lbl_checkedCalc()
# Calculate
#
ws.iyCalc()
# To save calculated spectrum and transmission
#
# output_file_formatSetAscii
# WriteXML( output_file_format, iy, "iyREFERENCE.xml" )
# WriteXML( output_file_format, iy_aux, "iy_auxREFERENCE.xml" )
# Check that results are OK with respect to an older reference calculation
#
ws.MatrixCreate("iy0")
#
ws.ReadXML(ws.iy0, "iyREFERENCE.xml")
ws.Compare(ws.iy, ws.iy0, 0.0001)
#
# Repeat, but running in "dispersion mode".
#
# Result shall here be identical.
#
ws.Copy(ws.iy0, ws.iy)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Freqloop)
ws.Copy(ws.iy_loop_freqs_agenda, ws.iy_loop_freqs_agenda__Transmission)
# No along-the-path-variables can here be included
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Optical depth"])
ws.iyCalc()
ws.Compare(ws.iy, ws.iy0, 1e-06)
