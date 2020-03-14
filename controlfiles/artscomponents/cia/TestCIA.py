# DEFINITIONS:  -*-sh-*-
#
# Calculate absorption cross sections for all CIA continua
# 2013-03-08 SAB

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.ReadXML(ws.abs_cia_data, "spectroscopy/cia/hitran2011/hitran_cia2012_adapted.xml.gz")
ws.Print(ws.abs_cia_data)
# Output:
# - Print
#   CIA tag; Spectral range [cm-1]; Temp range [K]; # of sets
#   N2-CIA-N2-0; 0.02 - 554.00; 40.00 - 400.00; 10
#   N2-CIA-N2-1; 1850.00 - 3000.09; 228.20 - 362.50; 10
#   N2-CIA-H2-0; 0.02 - 1886.00; 40.00 - 400.00; 10
#   N2-CIA-CH4-0; 0.02 - 1379.00; 40.00 - 400.00; 10
#   H2-CIA-H2-0; 20.00 - 10000.00; 200.00 - 3000.00; 113
#   H2-CIA-He-0; 20.00 - 20000.00; 200.00 - 9900.00; 334
#   H2-CIA-CH4-0; 0.02 - 1946.00; 40.00 - 400.00; 10
#   H2-CIA-H-0; 100.00 - 10000.00; 1000.00 - 2500.00; 4
#   He-CIA-H-0; 50.00 - 11000.00; 1500.00 - 10000.00; 10
#   O2-CIA-O2-0; 1150.00 - 1950.00; 193.40 - 353.40; 15
#   CO2-CIA-CO2-0; 1.00 - 250.00; 200.00 - 800.00; 7
#   CH4-CIA-CH4-0; 0.02 - 990.00; 40.00 - 400.00; 10
#   CH4-CIA-Ar-0; 1.00 - 697.00; 70.00 - 296.00; 5
# Define all available CIA tags:
ws.abs_speciesSet(
    species=[
        "N2-CIA-N2-0, N2-CIA-N2-1",
        "N2-CIA-H2-0",
        "N2-CIA-CH4-0",
        "H2-CIA-H2-0",
        "H2-CIA-He-0",
        "H2-CIA-CH4-0",
        "O2-CIA-O2-0",
        "CO2-CIA-CO2-0",
        "CH4-CIA-CH4-0",
        "CH4-CIA-Ar-0",
        "H",
        "Ar",
        "He",
    ]
)
ws.VectorNLinSpace(ws.f_grid, 100, 0.0, 75000000000000.0)
ws.VectorSet(ws.abs_p, np.array([100000.0]))
ws.VectorSet(ws.abs_t, np.array([310.0]))
# Set all VMR values to the same value:
ws.nelemGet(ws.nrows, ws.abs_species)
ws.nelemGet(ws.ncols, ws.abs_t)
# Set all VMRs to 1. This is unphysical, since actually the total
# should add up to 1. But does not play a role here. (As long as I
# don't forget about it.)
ws.MatrixSetConstant(ws.abs_vmrs, ws.nrows, ws.ncols, 1.0)
ws.ArrayOfIndexSet(ws.abs_species_active, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
# We don't have an abs_xsec_agenda here. This is specific
# to this test case that checks the consistency of AddCIA directly.
# In normal life, AddCIA is always called from inside the abs_xsec_agenda.
# We can therefore just set the abs_xsec_agenda_checked manually to true.
ws.FlagOn(ws.abs_xsec_agenda_checked)
ws.jacobianOff()
ws.abs_xsec_per_speciesInit()
ws.abs_xsec_per_speciesAddCIA()
# WriteXML ( "ascii", abs_xsec_per_species )
# WriteXML ( "ascii", f_grid )
# WriteXML ( "ascii", abs_t )
# WriteXML ( "ascii", abs_p )
# WriteXML ( "ascii", abs_species )
# WriteXML ( "zascii", abs_cia_data )
# Compare results to reference calculation:
ws.ArrayOfMatrixCreate("abs_xsec_per_species_reference")
ws.ReadXML(
    ws.abs_xsec_per_species_reference, "TestCIA.abs_xsec_per_species_reference.xml"
)
ws.Compare(
    ws.abs_xsec_per_species,
    ws.abs_xsec_per_species_reference,
    1e-40,
    "CIA Continua differ from reference calculation",
)
