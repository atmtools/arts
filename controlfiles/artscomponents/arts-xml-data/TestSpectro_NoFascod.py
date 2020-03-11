#
# Check functionality of Toolbox spectroscopic line data and consistency of
#  these with HITRAN spectroscopic catalogues.
#
# Version for non-Fascod species.
#
# CAUTION:
#   - HITRAN data is NOT included in the toolbox! Hence, for running the part
#      that uses HITRAN it is required that the user gets his/her own copy of
#      the HITRAN catalogue (and adapt path to the catalogue accordingly).
#   - This is a time consuming test!
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.StringCreate("ext")
ws.StringSet(ws.ext, "NoFascode")
# set atmospheric scenario
ws.StringCreate("atmcase")
ws.StringSet(ws.atmcase, "planetary_toolbox/validation/common/spectroscopy/data/")
# set absspecies to be considered (O2, N2 needed as broadening species with ARTS4)
ws.abs_speciesSet(
    species=[
        "C2H4",
        "C3H8",
        "CF4",
        "CH3Br",
        "CH3CN",
        "CH3OH",
        "ClONO2",
        "H2S",
        "H2SO4",
        "HCOOH",
        "HO2",
        "HOBr",
        "NO+",
        "O",
        "OCS",
        "SO",
        "O2",
        "N2",
        "H2O",
    ]
)
ws.execute_controlfile("TestSpectro_core.arts")
