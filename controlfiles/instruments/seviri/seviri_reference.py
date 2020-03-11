# DEFINITIONS:  -*-sh-*-
# ARTS setup file for SEVIRI reference simulations.
#
# Use this if you want to do exact calculations, for example to
# validate fast calculations.
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
#                            HIRS Channels 13-19 are shortwave
#                            channels. Simulating them with
#                            ARTS for thermal radiation only is
#                            pointless. So you probably want 0 to 11
#                            (zero based ARTS indexing)
# views       (ArrayOfIndex) Which views you want.
# hitran_file (String)       Name of HITRAN catalogue file.
# f_grid_spacing (Numeric)   Frequency grid spacing.
# ????????????   TODO   ???????????????
# Concerning the frequency grid spacing:
# I tested the spacing on HIRS channel 12. With 5e8 Hz spacing I got
# an RMS error of 0.01 K and a maximum error of 0.02 K. This should be
# good enough. We will have 8402 frequencies for Channel 12.
# A coarser spacing of 5e9 Hz gave an RMS error of 0.03 and a maximum
# error 0.05.
# All these number are against a reference calculation with 5e7 Hz
# spacing, which was shown to give nearly identical result to a somewhat
# coarser calculation.
# Recommendation:
# reference calculation: f_grid_spacing 5e8
# faster calculations:   f_grid_spacing 5e9

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# 1. General
# -----------
ws.execute_controlfile("seviri_general.arts")
# 2. Spectroscopy
# ----------------
ws.execute_controlfile("seviri_spectroscopy.arts")
# hitran separate because one might want to use a lookup table
ws.execute_controlfile("seviri_hitran.arts")
# 3. Sensor:
# -----------
ws.execute_controlfile("seviri_sensor_common.arts")
ws.execute_controlfile("seviri_sensor_reference.arts")
# Compact line list, to kick out lines that are outside there
# cutoff range for all frequencies.
ws.abs_lines_per_speciesCompact()
# ???????? TODO ????????????
# If there was a method to convert abs_lines_per_species back to
# abs_lines, then we could do that, and save that file. The only
# problem is, though, that the lines to include are not exactly the same
# for all the different versions of the HIRS sensor. Thus, one would
# also have to perform a merge of all the different resulting line
# lists. :-(
