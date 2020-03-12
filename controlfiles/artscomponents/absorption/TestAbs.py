# DEFINITIONS:  -*-sh-*-
#
# An example ARTS controlfile that calculates absorption
# coefficients.
# SAB 16.06.2000

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
# --------------------< A specific method >--------------------
#                      -------------------
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
# This defines the list of absorption species.
# Spectral lines will be assigned to the species in the order as the species
# are specified here. That means if you do ["H2O-181","H2O"], the last
# group H2O gets assigned all the H2O lines that do not fit in the
# first group.
#
# The continuum tags are special, since continua are not added by
# default. Thus, just selecting "H2O" will give you no continuum.
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
# Initialize the input variables of abs_coefCalcFromXsec from the Atm fields:
ws.AbsInputFromAtmFields()
# Create the frequency grid `f_grid':
ws.VectorNLinSpace(ws.f_grid, 100, 50000000000.0, 150000000000.0)
# Select species for non-linear treatment in lookup table:
ws.abs_speciesSet(abs_species=ws.abs_nls, species=[])
# abs_speciesSet( abs_species=abs_nls, species=["H2O-PWR98", "O2-PWR93"] )
# Set tempertature perturbation vector for lookup table:
ws.VectorSet(ws.abs_t_pert, [])
# VectorLinSpace( abs_t_pert, -10, 10, 1 )
# Set non-linear species VMR perturbation vector for lookup table:
ws.VectorSet(ws.abs_nls_pert, [])
# VectorNLogSpace( abs_nls_pert, 5, 0.01, 10 )
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.jacobianOff()
ws.abs_lookupCalc()
# Optionally write these to files:
# WriteXML( output_file_format, abs_t )
# WriteXML( output_file_format, vmrs )
# Write absorption coefficients to files:
# WriteXML( output_file_format, abs_coef )
# WriteXML( output_file_format, abs_coef_per_species )
# Write absorption lookup table to file:
ws.WriteXML(ws.output_file_format, ws.abs_lookup)
