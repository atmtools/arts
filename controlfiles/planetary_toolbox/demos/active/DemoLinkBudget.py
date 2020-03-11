# DEFINITIONS:  -*-sh-*-
#
# The script demonstrates a 3D radio link calculation. The settings mimics a
# satellite-to-ground link, but by changing the receiver position a radio
# occultation geometry can be ac hived. Calculations are done for three
# frequencies (f_centre+-f_bwidth/2), where dispersion is considered.
#
# The basic atmosphere is Fascode tropical. An ionosphere is added. By
# changing an include file, a non-ionosphere case can be obtained (resulting in
# faster calculations).
#
# Run with report level 0:
#   arts -r0 link_budget.arts
# to obtain formatted screen output with results.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# ============================================================================
# User settings
# ============================================================================
#
# Position of receiver
#
ws.VectorCreate("reciever")
# z [m], lat [deg], lon  [deg]
#
ws.VectorSet(ws.reciever, array([0.0, 15.0, 78.0]))
#
# Position of transmitter
#
# (If propagation close to the poles, modifications below can be needed.)
#
ws.VectorCreate("transmitter")
# z [m], lat [deg], lon  [deg]
#
ws.VectorSet(ws.transmitter, array([8.0e05, 1.7e01, 8.8e01]))
#
# Frequency
#
ws.NumericCreate("f_centre")
# Centre frequency [Hz]
ws.NumericCreate("f_bwidth")
# Bandwidth (full) [Hz]
#
ws.NumericSet(ws.f_centre, 2000000000.0)
ws.NumericSet(ws.f_bwidth, 100000000.0)
#
# Transmitted polarisation
#
# Do 'arts -d instrument_pol' for information on the coding of polarisation states.
# For example, V and H are coded as 5 and 5, respectively.
#
ws.ArrayOfIndexSet(ws.instrument_pol, [5])
#
# Partition function handling
#
ws.IndexCreate("bad_partition_functions_ok")
# ok to use extrapolation of partition functions (i.e. fitting range of
#  partition function paameterization does not cover the actual atmospheric
#  temperatures)?
# ( 0: not ok, 1: ok )
ws.IndexSet(ws.bad_partition_functions_ok, 0)
# ============================================================================
# Basic and RT settings
# ============================================================================
#
# General initialisation
#
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
#
# Propgation path agendas
#
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__TransmitterReceiverPath)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesEarth)
#
# Radiative transfer step lengths
#
ws.NumericSet(ws.ppath_lmax, 10000.0)
ws.NumericSet(ws.ppath_lraytrace, 500.0)
#
# Here we consider dispersion (overkill if not ionopshere included)
#
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Freqloop)
ws.Copy(ws.iy_sub_agenda, ws.iy_sub_agenda__Radiolink)
#
# Transmitted signal
#
ws.Copy(ws.iy_transmitter_agenda, ws.iy_transmitter_agenda__UnitPolIntensity)
#
# f_grid
#
ws.NumericCreate("grid_hwidth")
ws.NumericCreate("grid_start")
ws.NumericCreate("grid_stop")
#
ws.NumericScale(ws.grid_hwidth, ws.f_bwidth, 0.5)
ws.NumericAdd(ws.grid_stop, ws.f_centre, ws.grid_hwidth)
ws.NumericScale(ws.grid_hwidth, ws.grid_hwidth, -1.0)
ws.NumericAdd(ws.grid_start, ws.f_centre, ws.grid_hwidth)
#
ws.VectorNLinSpace(ws.f_grid, 3, ws.grid_start, ws.grid_stop)
# Postion and line-of-sight of sensor
#
ws.VectorSet(ws.rte_los, array([], dtype=float64))
# Dummy value
#
ws.Copy(ws.rte_pos, ws.reciever)
ws.Copy(ws.rte_pos2, ws.transmitter)
#
# Basic stuff
#
ws.IndexSet(ws.stokes_dim, 4)
ws.AtmosphereSet3D()
#
ws.jacobianOff()
ws.cloudboxOff()
ws.sensorOff()
#
ws.VectorSet(ws.lat_true, array([], dtype=float64))
ws.VectorSet(ws.lon_true, array([], dtype=float64))
# ============================================================================
# Generate lat_grid and lon_grid
# ============================================================================
#
# Both grids: 21 points covering = +-10 degs around receiver
#
ws.NumericCreate("grid0")
ws.IndexCreate("grid_length")
#
ws.NumericSet(ws.grid_hwidth, 15.0)
# Half-width [deg]
ws.IndexSet(ws.grid_length, 21)
# Total number of grid points
#
ws.Extract(ws.grid0, ws.reciever, 1)
ws.NumericAdd(ws.grid_stop, ws.grid0, ws.grid_hwidth)
ws.NumericScale(ws.grid_hwidth, ws.grid_hwidth, -1.0)
ws.NumericAdd(ws.grid_start, ws.grid0, ws.grid_hwidth)
#
ws.VectorNLinSpace(ws.lat_grid, ws.grid_length, ws.grid_start, ws.grid_stop)
#
ws.Extract(ws.grid0, ws.reciever, 2)
ws.NumericAdd(ws.grid_start, ws.grid0, ws.grid_hwidth)
ws.NumericScale(ws.grid_hwidth, ws.grid_hwidth, -1.0)
ws.NumericAdd(ws.grid_stop, ws.grid0, ws.grid_hwidth)
#
ws.VectorNLinSpace(ws.lon_grid, ws.grid_length, ws.grid_start, ws.grid_stop)
# ============================================================================
# Load atmosphere and absorption lines
# ============================================================================
# INCLUDE "earth_tropical_false3d_0to64km_lowfreq.arts"
ws.execute_controlfile("earth_tropical_false3d_0to1000km_lowfreq.arts")
# ============================================================================
# Calculate and output
# ============================================================================
#
# Run data checks
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=ws.bad_partition_functions_ok)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
#
# Auxilary variables
#
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Faraday rotation", "Extra path delay"])
#
# Determine transmitted signal (just for report purpose)
#
ws.MatrixCreate("iy_trans")
#
ws.VectorSet(ws.rtp_pos, array([], dtype=float64))
# Dummy values
ws.VectorSet(ws.rtp_los, array([], dtype=float64))
ws.AgendaExecute(ws.iy_transmitter_agenda)
ws.Copy(ws.iy_trans, ws.iy)
#
# Perform RT
#
ws.iyCalc()
#
# Extract total Faraday rotation
#
ws.Tensor4Create("farrot_total")
ws.Extract(ws.farrot_total, ws.iy_aux, 0)
#
# Report
#
ws.StringCreate("message")
ws.StringCreate("seperator")
ws.StringCreate("empty_line")
ws.StringSet(ws.seperator, "-")
ws.StringSet(ws.empty_line, "")
#
ws.StringSet(ws.message, "---- Result of radio link calculation ---")
ws.Print(ws.empty_line, 0)
ws.Print(ws.message, 0)
#
ws.StringSet(ws.message, "Frequencies [GHz]:")
ws.VectorCreate("f_ghz")
ws.VectorScale(ws.f_ghz, ws.f_grid, 1e-09)
ws.Print(ws.empty_line, 0)
ws.Print(ws.message, 0)
ws.Print(ws.seperator, 0)
ws.Print(ws.f_ghz, 0)
#
ws.StringSet(ws.message, "Transmitted Stokes vector (each row a frequency):")
ws.Print(ws.empty_line, 0)
ws.Print(ws.message, 0)
ws.Print(ws.seperator, 0)
ws.Print(ws.iy_trans, 0)
#
ws.StringSet(ws.message, "Received Stokes vector:")
ws.Print(ws.empty_line, 0)
ws.Print(ws.message, 0)
ws.Print(ws.seperator, 0)
ws.Print(ws.iy, 0)
#
ws.StringSet(ws.message, "Total Faraday rotation [deg]:")
ws.Print(ws.empty_line, 0)
ws.Print(ws.message, 0)
ws.Print(ws.seperator, 0)
ws.Print(ws.farrot_total, 0)
#
ws.Print(ws.empty_line, 0)
