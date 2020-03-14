#
# include file for deriving manually setting abs_species and deriving respective
# atmospheric field data
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.abs_speciesSet(species=ws.aostrtmp)
# get atmospheric data
ws.AtmRawRead(basename=ws.atmcase)
