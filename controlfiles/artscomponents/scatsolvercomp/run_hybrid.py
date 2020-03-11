import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Hybrid (here we must switch to another iy_main_agenda)
# Note that i_field comes fro last other called method
ws.Copy(ws.iy_main_agenda, ws.iy_hybrid_agenda)
ws.yCalc(y=ws.y_hybrid)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
