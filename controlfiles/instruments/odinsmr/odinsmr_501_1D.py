# As odinsmr_501.arts but the calculations are here done using a single
# measurement block.
#
# The creation of the sensor response matrix takes here a longer time
# and this option is most beneficial if the number of spectra is high
# and there is a considerable overlap of the antenna pattern for
# different tangent altitudes.
#
# A small part of an Odin-SMR scan is here considered. The tangent altitude
# spacing is 1.5 km for the shortest integration time.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.execute_controlfile("odinsmr_501_absorption.arts")
# Frequency grid
ws.ReadXML(ws.f_grid, "f_mono.SM_AC2ab.xml")
#
# Sensor characteristics
#
ws.IndexSet(ws.sensor_norm, 1)
# --- Antenna: ----------------------------------------------------------------
#
ws.IndexSet(ws.antenna_dim, 1)
# Number of pencil beam directions and tangent altitudes
ws.IndexCreate("n_pbeams")
ws.IndexSet(ws.n_pbeams, 55)
ws.IndexCreate("n_tan")
ws.IndexSet(ws.n_tan, 5)
# Final sensor_pos and los will here have length 1. The selection of
# sensor_los is arbitrary as long as the sum of sensor_los and
# antenna_dlos gives correct observation directions. It is here selected to
# put sensor_los to 0.
# The folder contains antenna pattern for different integration times.
# The pattern for smallest integration time is selected here.
ws.ReadXML(ws.antenna_response, "antenna.SM_AC2ab.875ms.xml")
# We need sensor_pos as input to *VectorZtanToZa1D*
ws.MatrixSetConstant(ws.sensor_pos, ws.n_pbeams, 1, 600000.0)
# Create vector of pencil beam directions
ws.VectorCreate("z_pbeams")
ws.VectorNLinSpace(ws.z_pbeams, ws.n_pbeams, 37000.0, 9000.0)
ws.VectorZtanToZa1D(
    ws.za_grid, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.z_pbeams
)
ws.Matrix1ColFromVector(ws.mblock_dlos_grid, ws.za_grid)
# We need sensor_pos as input to *VectorZtanToZa1D*
ws.MatrixSetConstant(ws.sensor_pos, ws.n_tan, 1, 600000.0)
# Create vector of zenith angles for selected tangent altitudes
ws.VectorCreate("z_tan")
ws.VectorNLinSpace(ws.z_tan, ws.n_tan, 26000.0, 20000.0)
ws.VectorCreate("za")
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.z_tan)
ws.Matrix1ColFromVector(ws.antenna_dlos, ws.za)
# Set *sensor_pos*
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 600000.0)
# Set *sensor_los*
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 0.0)
# ----- End of  antenna part --------------------------------------------------
# Mixer:
#
ws.ReadXML(ws.lo, "lo.SM_AC2ab.xml")
ws.ReadXML(ws.sideband_response, "sideband.SM_AC2ab.xml")
ws.StringSet(ws.sideband_mode, "upper")
# Spectrometer:
#
ws.ReadXML(ws.f_backend, "f_backend.SM_AC2ab.xml")
ws.ReadXML(ws.backend_channel_response, "backend_channel_response.xml")
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_responseMixer()
ws.sensor_responseIF2RF()
ws.sensor_responseBackend()
# End Arts
