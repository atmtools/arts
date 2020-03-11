# ARTS setup file for simulations of Odin-SMR measurements around 501.8 GHz.
# This band is part of the stratospheric mode, denoted as SM_AC2ab.
#
# The simulations are intended to match the operational settings, but there
# are some differences. See notes below.
#
# Input files are found in the folder tests/OdinSMR.
# See setup_input.m for source of input files.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.execute_controlfile("odinsmr_544_absorption.arts")
# Frequency grid
ws.ReadXML(ws.f_grid, "f_mono.SM_AC1e.xml")
#
# Sensor characteristics
#
ws.IndexSet(ws.sensor_norm, 1)
# Antenna:
#
ws.IndexSet(ws.antenna_dim, 1)
#
ws.VectorSet(
    ws.za_grid,
    array(
        [
            -0.2,
            -0.15,
            -0.1,
            -0.05,
            -0.04,
            -0.03,
            -0.02,
            -0.01,
            0.0,
            0.01,
            0.02,
            0.03,
            0.04,
            0.05,
            0.1,
            0.15,
            0.2,
        ]
    ),
)
ws.Matrix1ColFromVector(ws.mblock_dlos_grid, ws.za_grid)
#
ws.MatrixSetConstant(ws.antenna_dlos, 1, 1, 0.0)
# The folder contains antenna pattern for different integration times.
# The pattern for smallest integration time is selected here.
ws.ReadXML(ws.antenna_response, "antenna.SM_AC1e.875ms.xml")
# Mixer:
#
ws.ReadXML(ws.lo, "lo.SM_AC1e.xml")
ws.ReadXML(ws.sideband_response, "sideband.SM_AC1e.xml")
ws.StringSet(ws.sideband_mode, "lower")
# Spectrometer:
#
ws.ReadXML(ws.f_backend, "f_backend.SM_AC1e.xml")
ws.ReadXML(ws.backend_channel_response, "backend_channel_response.xml")
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_responseMixer()
ws.sensor_responseIF2RF()
ws.sensor_responseBackend()
# End Arts
