# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for AMSU-A simulations
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
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c3/sec3-3.htm
# Viewing angles
# There are 15 different angles, corresponding to one side of the AMSU-A scan.

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.MatrixSet(
    ws.antenna_dlos,
    array(
        [
            [-48.33],
            [-44.996897],
            [-41.663793],
            [-38.33069],
            [-34.997586],
            [-31.664483],
            [-28.331379],
            [-24.998276],
            [-21.665172],
            [-18.332069],
            [-14.998966],
            [-11.665862],
            [-8.332759],
            [-4.999655],
            [-1.666552],
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
            [5.2800000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.3596115e10, 1.1500000e08, 0.0000000e00, 1.7000000e08],
            [5.4400000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.4940000e10, 0.0000000e00, 0.0000000e00, 4.0000000e08],
            [5.5500000e10, 0.0000000e00, 0.0000000e00, 3.3000000e08],
            [5.7290344e10, 0.0000000e00, 0.0000000e00, 3.3000000e08],
            [5.7290344e10, 2.1700000e08, 0.0000000e00, 7.8000000e07],
            [5.7290344e10, 3.2220000e08, 4.8000000e07, 3.6000000e07],
            [5.7290344e10, 3.2220000e08, 2.2000000e07, 1.6000000e07],
            [5.7290344e10, 3.2220000e08, 1.0000000e07, 8.0000000e06],
            [5.7290344e10, 3.2220000e08, 4.5000000e06, 3.0000000e06],
            [8.9000000e10, 0.0000000e00, 0.0000000e00, 2.0000000e09],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation,
    [
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-H",
        "AMSU-H",
        "AMSU-V",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-H",
        "AMSU-V",
    ],
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, array([], dtype=float64))
# How many monochromatic frequencies to simulate the channel
ws.Touch(ws.met_mm_available_accuracies)
ws.Delete(ws.met_mm_available_accuracies)
# Number of frequencies for first accuracy (fast)
ws.ArrayOfIndexSet(ws.freq_number_tmp, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for second accuracy (normal)
ws.ArrayOfIndexSet(ws.freq_number_tmp, [1, 1, 1, 3, 3, 5, 5, 4, 4, 3, 3, 3, 4, 2, 1])
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for third accuracy (high)
ws.ArrayOfIndexSet(
    ws.freq_number_tmp, [1, 1, 1, 8, 8, 16, 15, 11, 13, 9, 9, 9, 11, 6, 1]
)
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for fourth accuracy (reference)
ws.ArrayOfIndexSet(
    ws.freq_number_tmp, [6, 1, 3, 23, 24, 44, 43, 34, 38, 26, 26, 27, 31, 17, 4]
)
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
ws.VectorSet(ws.freq_spacing_tmp, array([1.0e10, 1.0e09, 1.0e09, 1.0e09]))
ws.Delete(ws.freq_number_tmp)
ws.Extract(ws.met_mm_freq_number, ws.met_mm_available_accuracies, ws.met_mm_accuracy)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
