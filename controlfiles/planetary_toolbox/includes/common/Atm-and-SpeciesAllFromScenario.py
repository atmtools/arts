#
# include file for deriving abs_species and atm fields from available scenario
# data
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# derive abs species from scenario data
ws.abs_speciesDefineAllInScenario(basename=ws.atmcase)
# get corresponding atmospheric data
ws.AtmRawRead(basename=ws.atmcase)
