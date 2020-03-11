#
# include file for adding a single species to abs_species and add a
# corresponding vmr field
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.abs_speciesAdd(species=ws.aostrtmp)
ws.Copy(ws.casefull, ws.atmcase)
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
