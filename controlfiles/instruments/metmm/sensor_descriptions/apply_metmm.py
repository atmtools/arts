# DEFINITIONS:  -*-sh-*-
# Creates sensor response following the met-mm system
#
# This file expects that one met-mm definition file has been called, e.g.
# sensor_mhs.arts (which in turn requires prepare_metmm.arts been called
# before).
# In addition, the following workspace variables must exist and be set:
#
# channels    (ArrayOfIndex) Which channels you want.
#                            Note that this array uses zero-based ARTS
#                            indexing. It can be set to [-1] to select all
#                            channels.
# viewing_angles (ArrayOfIndex) Which views you want.

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.Select(ws.antenna_dlos, ws.antenna_dlos, ws.viewing_angles)
ws.Select(ws.met_mm_backend, ws.met_mm_backend, ws.channels)
ws.Select(ws.met_mm_polarisation, ws.met_mm_polarisation, ws.channels)
ws.Select(ws.met_mm_freq_number, ws.met_mm_freq_number, ws.channels)
ws.Select(ws.met_mm_freq_spacing, ws.met_mm_freq_spacing, ws.channels)
ws.f_gridMetMM(freq_spacing=ws.met_mm_freq_spacing, freq_number=ws.met_mm_freq_number)
ws.sensor_responseMetMM()
