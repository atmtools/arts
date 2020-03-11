# DEFINITIONS:  -*-sh-*-

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Observation geometry
# ---
# EPS-SG-B[123] shall be 817 km above the surface.  ICI shall be a conical
# scanner with a 53.1° incidence angle (at surface) = 45° nadir angle =
# 135° zenith angle
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 817000.0)
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 135.0)
# Sensor response setup
# ---
# For EPS-SG ICI.
#
# Source: EUMETSAT EPS-SG End User Requirements Document,
# EUM/PEPS/REQ/09/0151, v3B Draft, 8 April 2013.  Table 18.
#
# Note that currently (2013-12-11), channel 1 is documented differently
# at WMO-OSCAR, where channel 1 is 183.31±8.4 GHz with a 2x3000 MHz
# bandwidth (see http://www.wmo-sat.info/oscar/instruments/view/342 ).
ws.MatrixSet(
    ws.sensor_description_amsu,
    array(
        [
            [1.8331e11, 7.0000e09, 2.0000e09],
            [1.8331e11, 3.4000e09, 1.5000e09],
            [1.8331e11, 2.0000e09, 1.5000e09],
            [2.4320e11, 2.5000e09, 3.0000e09],
            [3.2515e11, 9.5000e09, 3.0000e09],
            [3.2515e11, 3.5000e09, 2.4000e09],
            [3.2515e11, 1.5000e09, 1.6000e09],
            [4.4800e11, 7.2000e09, 3.0000e09],
            [4.4800e11, 3.0000e09, 2.0000e09],
            [4.4800e11, 1.4000e09, 1.2000e09],
            [6.6400e11, 4.2000e09, 5.0000e09],
        ]
    ),
)
ws.sensor_responseSimpleAMSU()
