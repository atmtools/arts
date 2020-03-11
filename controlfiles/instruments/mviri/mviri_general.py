# Included by mviri_fast.arts and mviri_reference.arts

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.output_file_formatSetZippedAscii()
# ???????? FIXME ???????
# HIRS spectral response functions assume radiance.
# ---
ws.StringSet(ws.iy_unit, "1")
