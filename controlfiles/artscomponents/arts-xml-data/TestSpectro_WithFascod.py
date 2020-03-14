#
# Check functionality of Toolbox spectroscopic line data and consistency of
#  these with HITRAN spectroscopic catalogues.
#
# Version for Fascod species.
#
# CAUTION:
#   - HITRAN data is NOT included in the toolbox! Hence, for running the part
#      that uses HITRAN it is required that the user gets his/her own copy of
#      the HITRAN catalogue (and adapt path to the catalogue accordingly).
#   - This is a time consuming test!
#
# Jana Mendrok 2013-02-26

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.StringCreate("ext")
ws.StringSet(ws.ext, "WithFascode")
# set atmospheric scenario
ws.StringCreate("atmcase")
ws.StringSet(ws.atmcase, "planets/Earth/Fascod/tropical/tropical")
# derive abs species from scenario data
ws.abs_speciesDefineAllInScenario(basename=ws.atmcase)
ws.execute_controlfile("TestSpectro_core.arts")
