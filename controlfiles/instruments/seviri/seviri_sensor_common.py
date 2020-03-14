# Included by seviri_fast.arts and seviri_reference.arts

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Definition of sensor position and LOS
# ---
ws.ReadXML(ws.sensor_los, "seviri.sensor_los.xml")
# Select those views that are requested by the user
ws.Select(ws.sensor_los, ws.sensor_los, ws.views)
ws.nrowsGet(ws.nrows, ws.sensor_los)
ws.ncolsGet(ws.ncols, ws.sensor_los)
ws.MatrixSetConstant(ws.sensor_pos, ws.nrows, ws.ncols, 850000.0)
# Start sensor response setup
# ---
# Normalise the sensor response
# ---
ws.IndexSet(ws.sensor_norm, 1)
# Antenna
# ---
ws.AntennaOff()
# See setup_input.m for details around other sensor variables
# Construct names of sensor description files:
# Nominal channel frequencies:
ws.StringCreate("f_backend_file")
ws.StringSet(ws.f_backend_file, "")
ws.Append(ws.f_backend_file, ws.satellite)
ws.StringCreate("dummy")
ws.StringSet(ws.dummy, "_SEVIRI.f_backend.xml")
ws.Append(ws.f_backend_file, ws.dummy)
