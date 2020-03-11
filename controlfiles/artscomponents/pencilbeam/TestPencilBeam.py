# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of performing a single monochromatic pencil beam
# calculation, without involving any sensor characteristics. That is, how to
# calculate monochromatic pencil beam spectra without using yCalc.
#
# The case treats 1D limb sounding. The result is converted to a measurement
# vector, y, holding brightness temperatures.
#
# 2012-02-17, Patrick Eriksson

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
# Frequency grid
#
ws.VectorNLinSpace(ws.f_grid, 201, 325000000000.0, 327000000000.0)
# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
#
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
# Definition of species
#
ws.abs_speciesSet(species=["H2O-PWR98", "N2-SelfContStandardType", "O2-PWR93"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Dimensionality of the atmosphere
#
ws.AtmosphereSet1D()
# Atmospheric profiles
#
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalc()
# Get ground altitude (z_surface) from z_field
ws.Extract(ws.z_surface, ws.z_field, 0)
# Definition of position and LOS (simulating limb sounding from 600 km)
#
ws.VectorSet(ws.rte_pos, array([600000.0]))
ws.VectorSet(ws.rte_los, array([113.3]))
ws.VectorSet(ws.rte_pos2, array([], dtype=float64))
# No transmitter involved
# Define auxiliary data (the order between the quantities is free)
#
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Radiative background", "Optical depth"])
# No jacobian calculation
#
ws.jacobianOff()
# No scattering
#
ws.cloudboxOff()
# Perform RT calculations
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.lbl_checkedCalc()
#
ws.iyCalc()
# Convert to Tb
#
ws.StringSet(ws.iy_unit, "RJBT")
#
ws.iyApplyUnit()
# To save calculated spectrum and transmission
#
# output_file_formatSetAscii
# WriteXML( output_file_format, f_grid, "f_grid.xml" )
# WriteXML( output_file_format, iy, "iyREFERENCE.xml" )
# WriteXML( output_file_format, iy_aux, "iy_auxREFERENCE.xml" )
# WriteXML( output_file_format, ppath, "ppath.xml" )
# Check that results are OK with respect to an older reference calculation
#
ws.MatrixCreate("iy0")
#
ws.ReadXML(ws.iy0, "iyREFERENCE.xml")
ws.Compare(ws.iy, ws.iy0, 0.01)
