import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# RT4
ws.Copy(ws.za_grid_copy, ws.za_grid)
ws.Copy(ws.aa_grid_copy, ws.aa_grid)
ws.RT4Calc(
    nstreams=16,
    auto_inc_nstreams=32,
    robust=1,
    quad_type="l",
    pfct_aa_grid_size=37,
    pfct_method="interpolate",
)
ws.yCalc(y=ws.y_rt4)
ws.Copy(ws.za_grid, ws.za_grid_copy)
ws.Copy(ws.aa_grid, ws.aa_grid_copy)
