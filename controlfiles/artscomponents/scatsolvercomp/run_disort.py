import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# DISORT
ws.DisortCalc(pfct_method="interpolate")
ws.yCalc(y=ws.y_disort)
