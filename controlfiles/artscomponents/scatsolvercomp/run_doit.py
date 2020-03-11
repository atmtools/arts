import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# DOIT
ws.DoitInit()
ws.DoitGetIncoming(rigorous=0)
ws.cloudbox_fieldSetClearsky()
ws.DoitCalc()
ws.yCalc(y=ws.y_doit)
