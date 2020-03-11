# DEFINITIONS:  -*-sh-*-

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Sensor response setup (reset for CloudIce)
# ---
ws.MatrixSet(
    ws.sensor_description_amsu,
    array(
        [
            [1.8331e11, 2.0000e08, 2.0000e08],
            [1.8331e11, 1.0000e09, 5.0000e08],
            [1.8331e11, 3.0000e09, 1.0000e09],
            [1.8331e11, 5.0000e09, 1.5000e09],
            [1.8331e11, 7.0000e09, 2.0000e09],
            [1.8331e11, 1.1000e10, 3.0000e09],
            [2.4320e11, 2.5000e09, 3.0000e09],
            [3.2515e11, 1.5000e09, 1.6000e09],
            [3.2515e11, 3.5000e09, 2.4000e09],
            [3.2515e11, 9.5000e09, 3.0000e09],
            [4.4800e11, 1.4000e09, 1.2000e09],
            [4.4800e11, 3.0000e09, 2.0000e09],
            [4.4800e11, 7.2000e09, 3.0000e09],
            [6.6400e11, 4.2000e09, 5.0000e09],
        ]
    ),
)
ws.sensor_responseSimpleAMSU()
# Replace f_grid and sensor_response by optimized ones.
ws.ReadXML(ws.f_grid, "instruments/smr/cloudice.f_grid_fast.xml")
ws.ReadXML(ws.sensor_response, "instruments/smr/cloudice.sensor_response_fast.xml")
ws.ReadXML(ws.sensor_response_f, "instruments/smr/cloudice.sensor_response_f_fast.xml")
