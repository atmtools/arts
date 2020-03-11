# DEFINITIONS:  -*-sh-*-
# ARTS setup file for fast AVHRR simulations.
#
# This expects a number of workspace variables to exist and to be set:
#
# satellite   (String) 	     The name of the satellite. Is
#                      	     used internally to construct
#                      	     file names of the instrument
#                      	     description files
#                      	     f_backend_file and
#                      	     backend_channel_response_file.
# channels    (ArrayOfIndex) Which channels you want.
#                            Only AVHRR channels 3B, 4 and 5 are
#                            implemented. As ARTS is zero-based so to
#                            get a specific channel, consider:
#                            ARTS 0 -> AVHRR 3B
#                            ARTS 1 -> AVHRR 4
#                            ARTS 2 -> AVHRR 5
# views       (ArrayOfIndex) Which views you want. Implemented for AVHRR
#                            FRAC/LAC, e.g. 1024 views on each side of the
#                            scan.
# hitran_file (String)       Name of HITRAN catalogue file.

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# 1. General
# -----------
# this part is the same as for HIRS
ws.execute_controlfile("hirs/hirs_general.arts")
ws.execute_controlfile("hirs/hirs_spectroscopy.arts")
ws.execute_controlfile("hirs/hirs_hitran.arts")
# 3. Sensor:
# -----------
ws.execute_controlfile("avhrr_sensor_common.arts")
ws.execute_controlfile("avhrr_sensor_fast.arts")
# See comment in hirs/hirs_reference.arts
ws.abs_lines_per_speciesCompact()
