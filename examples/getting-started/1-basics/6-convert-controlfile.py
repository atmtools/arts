
import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda
ws = Workspace(verbosity=0)
ws.VectorNLinSpace(ws.f_grid, 5, 320000000000.0, 322000000000.0)
