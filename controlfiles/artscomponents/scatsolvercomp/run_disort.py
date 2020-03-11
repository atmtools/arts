import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# DISORT
ws.DisortCalc(pfct_method="interpolate")
ws.yCalc(y=ws.y_disort)
