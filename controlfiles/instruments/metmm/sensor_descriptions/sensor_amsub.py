# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for AMSU-B simulations
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
# Sensor characteristics based on KLM User's Guide at
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c3/sec3-4.htm
# Viewing angles
# There are 45 different angles, corresponding to one side of the AMSU-B Scan.

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.MatrixSet(
    ws.antenna_dlos,
    np.array(
        [
            [-48.95],
            [-47.85],
            [-46.75],
            [-45.65],
            [-44.55],
            [-43.45],
            [-42.35],
            [-41.25],
            [-40.15],
            [-39.05],
            [-37.95],
            [-36.85],
            [-35.75],
            [-34.65],
            [-33.55],
            [-32.45],
            [-31.35],
            [-30.25],
            [-29.15],
            [-28.05],
            [-26.95],
            [-25.85],
            [-24.75],
            [-23.65],
            [-22.55],
            [-21.45],
            [-20.35],
            [-19.25],
            [-18.15],
            [-17.05],
            [-15.95],
            [-14.85],
            [-13.75],
            [-12.65],
            [-11.55],
            [-10.45],
            [-9.35],
            [-8.25],
            [-7.15],
            [-6.05],
            [-4.95],
            [-3.85],
            [-2.75],
            [-1.65],
            [-0.55],
        ]
    ),
)
# Sensor response setup
# ---
ws.MatrixSet(
    ws.met_mm_backend,
    np.array(
        [
            [8.9000e10, 9.0000e08, 0.0000e00, 1.0000e09],
            [1.5000e11, 9.0000e08, 0.0000e00, 1.0000e09],
            [1.8331e11, 1.0000e09, 0.0000e00, 5.0000e08],
            [1.8331e11, 3.0000e09, 0.0000e00, 1.0000e09],
            [1.8331e11, 7.0000e09, 0.0000e00, 2.0000e09],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation, ["AMSU-V", "AMSU-V", "AMSU-V", "AMSU-V", "AMSU-V"]
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, [])
# How many monochromatic frequencies to simulate the channel
ws.Touch(ws.met_mm_available_accuracies)
ws.Delete(ws.met_mm_available_accuracies)
# Number of frequencies for first accuracy (fast)
ws.ArrayOfIndexSet(ws.freq_number_tmp, [1, 1, 1, 1, 1])
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for second accuracy (normal)
ws.ArrayOfIndexSet(ws.freq_number_tmp, [1, 2, 2, 2, 3])
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for third accuracy (high)
ws.ArrayOfIndexSet(ws.freq_number_tmp, [1, 18, 20, 7, 10])
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for fourth accuracy (reference)
ws.ArrayOfIndexSet(ws.freq_number_tmp, [2, 23, 67, 19, 25])
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
ws.VectorSet(ws.freq_spacing_tmp, np.array([1.0e10, 1.0e09, 1.0e09, 1.0e09]))
ws.Delete(ws.freq_number_tmp)
ws.Extract(ws.met_mm_freq_number, ws.met_mm_available_accuracies, ws.met_mm_accuracy)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
