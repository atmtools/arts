# DEFINITIONS:  -*-sh-*-
# ARTS setup file for fast SEVIRI simulations.
#
# Use this if you want to do fast SEVIRI simulations.  It should be a
# good approximation of the exact simulations.
#
# This expects a number of workspace variables to exist and to be set:
#
# satellite   (String) 	     The name of the satellite. Is
#                      	     used internally to construct
#                      	     file names of the instrument
#                      	     description files.
# channels    (ArrayOfIndex) Which channels you want.
#                            HIRS Channels 13-19 are shortwave
#                            channels. Simulating them with
#                            ARTS for thermal radiation only is
#                            pointless. So you probably want 0 to 11
#                            (zero based ARTS indexing)
# views       (ArrayOfIndex) Which views you want.
# hitran_file (String)       Name of HITRAN catalogue file.

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# 1. General
# -----------
ws.execute_controlfile("seviri_general.arts")
# 2. Spectroscopy
# ----------------
ws.execute_controlfile("seviri_spectroscopy.arts")
# deriving abs_lines separate, because one might want to use a lookup table
#
# no use of HITRAN in test cases
# (uncomment when you want to use HITRAN)
# INCLUDE "instruments/hirs/hirs_hitran.arts"
# instead we use a HITRAN-extract converted to ARTSCAT
# (outcomment when you want to use HITRAN)
ws.ReadXML(ws.abs_lines, "testdata/abs_lines_IR.xml.gz")
ws.abs_lines_per_speciesCreateFromLines()
# 3. Sensor:
# -----------
ws.execute_controlfile("seviri_sensor_common.arts")
ws.execute_controlfile("seviri_sensor_fast.arts")
# see comment in hirs_reference.arts
ws.abs_lines_per_speciesCompact()
