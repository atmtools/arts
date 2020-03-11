# DEFINITIONS:  -*-sh-*-
# An example ARTS controlfile that calculates absorption
# coefficients with Doppler shift.
# 2011-05-13 Stefan Buehler

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Turn on frequency interpolation
# (we could also try higher order interpolation here)
ws.IndexSet(ws.abs_f_interp_order, 1)
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# On-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Read the spectroscopic line data from the ARTS catalogue and
# create the workspace variable `lines'.
ws.ReadARTSCAT(
    abs_lines=ws.abs_lines, filename="lines.xml", fmin=1000000000.0, fmax=200000000000.0
)
# This test catalogue was generated from the HITRAN catalogue in the
# following way:
# abs_linesReadFromHitran( abs_lines,
#        "PATH_TO_ARTS-DATA/spectroscopy/hitran96/hitran96_lowfreq.par",
#        1e9,
#        200e9 )
# Define the list of absorption species:
ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR93", "N2-SelfContStandardType"])
# Alternatively select all species that we can find in the scenario:
# abs_speciesDefineAllInScenario( basename="testdata/tropical" )
# This separates the lines into the different tag groups and creates
# the workspace variable `abs_lines_per_species':
ws.abs_lines_per_speciesCreateFromLines()
# Dimensionality of the atmosphere
#
ws.AtmosphereSet1D()
ws.VectorNLogSpace(ws.p_grid, 10, 100000.0, 10.0)
# Atmospheric profiles
ws.AtmRawRead(basename="testdata/tropical")
# Now interpolate all the raw atmospheric input onto the pressure
# grid and create the atmospheric variables `t_field', `z_field', `vmr_field'
ws.AtmFieldsCalc()
# Initialize the input variables of propmat_clearsky_fieldCalc from the Atm fields:
ws.AbsInputFromAtmFields()
# Create the frequency grid `f_grid':
ws.VectorNLinSpace(ws.f_grid, 500, 50000000000.0, 150000000000.0)
# Calculate field of absorption coefficients:
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.lbl_checkedCalc()
ws.jacobianOff()
ws.propmat_clearsky_fieldCalc()
# Write out:
ws.Tensor7Create("abs_field_nodoppler")
ws.Copy(ws.abs_field_nodoppler, ws.propmat_clearsky_field)
ws.WriteXML(ws.output_file_format, ws.abs_field_nodoppler)
# Now with Doppler (still LBL)!
# Create vector of Doppler shift values
ws.VectorCreate("doppler")
# Make ramp with same dimension as p_grid
ws.nelemGet(ws.nelem, ws.p_grid)
ws.VectorNLinSpace(ws.doppler, ws.nelem, 0.0, 1000000000.0)
# Calculate field of absorption coefficients:
ws.propmat_clearsky_fieldCalc(doppler=ws.doppler)
# Write out:
ws.Tensor7Create("abs_field_doppler")
ws.Copy(ws.abs_field_doppler, ws.propmat_clearsky_field)
ws.WriteXML(ws.output_file_format, ws.abs_field_doppler)
# Now with Doppler and lookup table!
# Make the absorption lookupt table frequency grid cover a larger range
# (we will move outside the original f_grid due to Doppler).
# At the same time, make it also denser for better accuracy.
ws.VectorCreate("f_grid_backup")
ws.Copy(ws.f_grid_backup, ws.f_grid)
ws.VectorNLinSpace(ws.f_grid, 1000, 49000000000.0, 151000000000.0)
ws.abs_lookupSetup()
ws.abs_lookupCalc()
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Restore original frequency grid
ws.Copy(ws.f_grid, ws.f_grid_backup)
# Calculate field of absorption coefficients:
ws.propmat_clearsky_fieldCalc(doppler=ws.doppler)
# Write out:
ws.Tensor7Create("abs_field_doppler_lookup")
ws.Copy(ws.abs_field_doppler_lookup, ws.propmat_clearsky_field)
ws.WriteXML(ws.output_file_format, ws.abs_field_doppler_lookup)
# Write out aux variables for plotting:
ws.WriteXML(ws.output_file_format, ws.p_grid)
ws.WriteXML(ws.output_file_format, ws.f_grid)
