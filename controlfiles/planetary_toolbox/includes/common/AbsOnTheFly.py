#
# include file for getting absorption on-the-fly
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Read the line files and prepare line data
ws.abs_linesReadFromSplitArtscat(
    ws.abs_lines, ws.abs_species, "spectroscopy/Perrin/", 0.0, 3000000000000.0
)
ws.abs_lines_per_speciesCreateFromLines()
# do on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
