# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for MHS simulations
#
# This requires to run prepare_metmm.arts beforehand.
#
# This expects the following workspace variables to exist and to be set:
#    met_mm_accuracy (Index)    Selection of accuracy level.
#
# The following variables are set:
#    antenna_dlos
#    met_mm_backend
#    met_mm_polarisation
#    met_mm_freq_number
#    met_mm_freq_spacing
#    met_mm_antenna

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Sensor characteristics based on KLM User's Guide at
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c3/sec3-9.htm
# Viewing angles
# There are 45 different angles, corresponding to one side of the MHS scan.
ws.MatrixSet(
    ws.antenna_dlos,
    array(
        [
            [-49.444444],
            [-48.333333],
            [-47.222222],
            [-46.111111],
            [-45.0],
            [-43.888888],
            [-42.777777],
            [-41.666666],
            [-40.555555],
            [-39.444444],
            [-38.333333],
            [-37.222222],
            [-36.111111],
            [-35.0],
            [-33.888889],
            [-32.777777],
            [-31.666666],
            [-30.555555],
            [-29.444444],
            [-28.333333],
            [-27.222222],
            [-26.111111],
            [-25.0],
            [-23.888889],
            [-22.777778],
            [-21.666666],
            [-20.555555],
            [-19.444444],
            [-18.333333],
            [-17.222222],
            [-16.111111],
            [-15.0],
            [-13.888889],
            [-12.777778],
            [-11.666667],
            [-10.555555],
            [-9.444444],
            [-8.333333],
            [-7.222222],
            [-6.111111],
            [-5.0],
            [-3.888889],
            [-2.777778],
            [-1.666667],
            [-0.555556],
        ]
    ),
)
# Sensor response setup
# ---
ws.MatrixSet(
    ws.met_mm_backend,
    array(
        [
            [8.90000e10, 0.00000e00, 0.00000e00, 2.80000e09],
            [1.57000e11, 0.00000e00, 0.00000e00, 2.80000e09],
            [1.83311e11, 1.00000e09, 0.00000e00, 5.00000e08],
            [1.83311e11, 3.00000e09, 0.00000e00, 1.00000e09],
            [1.90311e11, 0.00000e00, 0.00000e00, 2.20000e09],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation, ["AMSU-V", "AMSU-V", "AMSU-H", "AMSU-H", "AMSU-V"]
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, array([], dtype=float64))
ws.ArrayOfIndexSet(ws.met_mm_freq_number, [12, 12, 12, 12, 12])
ws.VectorSet(ws.freq_spacing_tmp, array([1.0e10, 1.0e09, 1.0e09, 1.0e09]))
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
