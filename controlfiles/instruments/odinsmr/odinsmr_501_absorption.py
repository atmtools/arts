import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("odinsmr.arts")
# A number of species with very weak emission are here neglected
ws.abs_speciesSet(
    species=[
        "H2O,H2O-ForeignContStandardType,H2O-SelfContStandardType",
        "N2-SelfContMPM93",
        "O2,O2-PWR98",
        "O3",
        "ClO",
        "N2O",
    ]
)
# Line file
ws.ReadARTSCAT(filename="linefile.SM_AC2ab.xml", fmin=0.0, fmax=2000000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
