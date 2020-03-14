# Included by mviri_reference.arts

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Channel response functions:
ws.StringCreate("backend_channel_response_file")
ws.StringSet(ws.backend_channel_response_file, "")
ws.Append(ws.backend_channel_response_file, ws.satellite)
# FIXME: Remove this when using constant in Append works:
ws.StringSet(ws.dummy, "_MVIRI.backend_channel_response.xml")
#
ws.Append(ws.backend_channel_response_file, ws.dummy)
# Spectrometer:
#
ws.ReadXML(ws.f_backend, ws.f_backend_file)
ws.ReadXML(ws.backend_channel_response, ws.backend_channel_response_file)
# Select the desired channels:
#
ws.Select(ws.f_backend, ws.f_backend, ws.channels)
ws.Select(ws.backend_channel_response, ws.backend_channel_response, ws.channels)
# Frequency grid:
#
ws.f_gridFromSensorHIRS(
    ws.f_grid, ws.f_backend, ws.backend_channel_response, ws.f_grid_spacing
)
# Construct sensor response from all these inputs:
#
ws.sensor_responseInit()
ws.sensor_responseBackend()
# End of sensor response setup
