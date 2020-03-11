# Included by seviri_fast.arts

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Fast frequency grid
ws.StringCreate("f_grid_file")
ws.StringSet(ws.f_grid_file, "")
ws.Append(ws.f_grid_file, ws.satellite)
ws.StringSet(ws.dummy, "_SEVIRI.f_grid_fast.xml")
ws.Append(ws.f_grid_file, ws.dummy)
# Weights associated with each frequency
ws.StringCreate("weights_file")
ws.StringSet(ws.weights_file, "")
ws.Append(ws.weights_file, ws.satellite)
ws.StringSet(ws.dummy, "_SEVIRI.W_fast.xml")
ws.Append(ws.weights_file, ws.dummy)
# Spectrometer:
#
ws.ReadXML(ws.f_backend, ws.f_backend_file)
# Read in optimized frequency grid
ws.ReadXML(ws.f_grid, ws.f_grid_file)
# Read in and WMRF weights:
ws.ReadXML(ws.wmrf_weights, ws.weights_file)
# Select only the active channels:
ws.Copy(ws.wmrf_channels, ws.channels)
# The method acts on several variables:
# f_grid, wmrf_weights, and f_backend
ws.WMRFSelectChannels()
# Initialize sensor variables.
ws.sensor_responseInit()
# Add WMRF weigths to sensor response:
ws.sensor_responseWMRF()
# End of sensor response setup
# Compact line list, to kick out lines that are outside there
# cutoff range for all frequencies.
# abs_lines_per_speciesCompact
# ???????? TODO ????????????
# If there was a method to convert abs_lines_per_species back to
# abs_lines, then we could do that, and save that file. The only
# problem is, though, that the lines to include are not exactly the same
# for all the different versions of the HIRS sensor. Thus, one would
# also have to perform a merge of all the different resulting line
# lists. :-(
