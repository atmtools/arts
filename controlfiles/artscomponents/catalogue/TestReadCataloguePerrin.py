#
# Simple test to read all lines for all species from the
# Perrin catalogue in arts-xml-data.
#
# Author: Oliver Lemke

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# set case (for location where to look for data)
ws.StringCreate("case")
ws.StringSet(ws.case, "spectroscopy/Perrin/")
# derive all species we have line data for
ws.abs_speciesDefineAllInScenario(basename=ws.case)
# alternatively, specify the species by hand
# abs_speciesSet(species=
#  ["C2H2", "CH4", "H2CO", "HCN", "HO2", "NO+", "OCS", "C2H4", "CO", "H2O",
#   "HCOOH", "HOBr", "NO", "OH", "C2H6", "CO2", "H2O2", "HCl", "HOCl", "NO2",
#   "PH3", "C3H8", "COF2", "H2S", "HF", "N2", "O", "SF6", "CH3CN", "ClO",
#   "H2SO4", "HI", "N2O", "O2", "SO", "CH3Cl", "ClONO2", "HBr", "HNO3",
#   "NH3", "O3", "SO2"])
# Read the line files and prepare line data
ws.abs_linesReadFromSplitArtscat(
    ws.abs_lines, ws.abs_species, "spectroscopy/Perrin/", 0.0, 10000000000000.0
)
ws.abs_lines_per_speciesCreateFromLines()
