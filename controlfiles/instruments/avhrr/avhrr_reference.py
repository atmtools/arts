# DEFINITIONS:  -*-sh-*-
# ARTS setup file for AVHRR reference simulations.
#
# Based on hirs_reference.arts.
#
# So far this is the only include file for AVHRR, but in the future a fast
# setup may be added, in analogy to hirs_fast.arts.
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
# f_grid_spacing (Numeric)   Frequency grid spacing.

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# 1. General
# -----------
# this part is the same as for HIRS
ws.execute_controlfile("instruments/hirs/hirs_general.arts")
ws.execute_controlfile("instruments/hirs/hirs_spectroscopy.arts")
# no use of HITRAN in test cases
# (uncomment when you want to use HITRAN)
# INCLUDE "instruments/hirs/hirs_hitran.arts"
# instead we use a HITRAN-extract converted to ARTSCAT
# (outcomment when you want to use HITRAN)
ws.ReadXML(ws.abs_lines, "testdata/abs_lines_IR.xml.gz")
ws.abs_lines_per_speciesCreateFromLines()
# 3. Sensor:
# -----------
ws.execute_controlfile("avhrr_sensor_common.arts")
ws.execute_controlfile("avhrr_sensor_reference.arts")
# See comment in hirs/hirs_reference.arts
ws.abs_lines_per_speciesCompact()
