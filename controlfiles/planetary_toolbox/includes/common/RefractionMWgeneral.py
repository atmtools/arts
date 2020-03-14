#
# include file for doing path calculations without considering refraction
# (purely geometric path)
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
# use microwave refractive index method for general atmospheres (from Newell&Baird)
ws.Copy(ws.refr_index_agenda, ws.refr_index_agenda__GasMicrowavesGeneral)
