# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for ATMS simulations
#
# The following variables are set:
#
#    antenna_dlos
#    met_mm_backend
#    met_mm_polarisation
#    met_mm_freq_number
#    met_mm_antenna

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Sensor characteristics based on the WMO-OSCAR at
# https://www.wmo-sat.info/oscar/instruments/view/53
# OSCAR stands for Observing Systems Capability Analysis and Review Tool
# Central frequency of Channel 16 (88.2GHz) is unclear (different sources
# give different numbers: either 88.2 or 89.5.
# Viewing angles
# There are 48 different angles, corresponding to one side of the ATMS scan.
ws.MatrixSet(
    ws.antenna_dlos,
    array(
        [
            [-52.777777],
            [-51.666666],
            [-50.555555],
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
            [2.3800000e10, 0.0000000e00, 0.0000000e00, 2.7000000e08],
            [3.1400000e10, 0.0000000e00, 0.0000000e00, 1.8000000e08],
            [5.0300000e10, 0.0000000e00, 0.0000000e00, 1.8000000e08],
            [5.1760000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.2800000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.3596000e10, 1.1500000e08, 0.0000000e00, 1.7000000e08],
            [5.4400000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.4940000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.5500000e10, 0.0000000e00, 0.0000000e00, 3.3000000e08],
            [5.7290344e10, 0.0000000e00, 0.0000000e00, 3.3000000e08],
            [5.7290344e10, 2.1700000e08, 0.0000000e00, 7.8000000e07],
            [5.7290344e10, 3.2220000e08, 4.8000000e07, 3.6000000e07],
            [5.7290344e10, 3.2220000e08, 2.2000000e07, 1.6000000e07],
            [5.7290344e10, 3.2220000e08, 1.0000000e07, 8.0000000e06],
            [5.7290344e10, 3.2220000e08, 4.5000000e06, 3.0000000e06],
            [8.8200000e10, 0.0000000e00, 0.0000000e00, 5.0000000e09],
            [1.6550000e11, 0.0000000e00, 0.0000000e00, 3.0000000e09],
            [1.8331000e11, 7.0000000e09, 0.0000000e00, 2.0000000e09],
            [1.8331000e11, 4.5000000e09, 0.0000000e00, 2.0000000e09],
            [1.8331000e11, 3.0000000e09, 0.0000000e00, 1.0000000e09],
            [1.8331000e11, 1.8000000e09, 0.0000000e00, 1.0000000e09],
            [1.8331000e11, 1.0000000e09, 0.0000000e00, 5.0000000e08],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation,
    [
        "AMSU-V",
        "AMSU-V",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-V",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
    ],
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, array([], dtype=float64))
ws.ArrayOfIndexSet(
    ws.met_mm_freq_number,
    [
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
        12,
    ],
)
ws.VectorSet(
    ws.freq_spacing_tmp,
    array(
        [
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e09,
            1.0e10,
        ]
    ),
)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
