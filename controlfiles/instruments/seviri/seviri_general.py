# Included by seviri_fast.arts and seviri_reference.arts

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.output_file_formatSetZippedAscii()
# ???????? FIXME ???????
# SEVIRI spectral response functions assume radiance.
# ---
ws.StringSet(ws.iy_unit, "1")
