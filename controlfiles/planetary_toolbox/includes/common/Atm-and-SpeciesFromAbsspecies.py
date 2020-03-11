#
# include file for deriving manually setting abs_species and deriving respective
# atmospheric field data
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.abs_speciesSet(species=ws.aostrtmp)
# get atmospheric data
ws.AtmRawRead(basename=ws.atmcase)
