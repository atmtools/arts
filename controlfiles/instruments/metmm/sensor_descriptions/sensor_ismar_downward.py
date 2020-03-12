# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for ISMAR simulations
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
# Viewing angles
ws.MatrixSet(
    ws.antenna_dlos,
    np.array(
        [
            [-180.0],
            [-170.0],
            [-110.0],
            [-50.0],
            [-40.0],
            [-30.0],
            [-20.0],
            [-10.0],
            [0.0],
            [10.0],
            [70.0],
            [140.0],
            [150.0],
            [160.0],
            [170.0],
        ]
    ),
)
# Sensor response setup
# ---
ws.MatrixSet(
    ws.met_mm_backend,
    np.array(
        [
            [1.187503e11, 1.100000e09, 0.000000e00, 4.000000e08],
            [1.187503e11, 1.500000e09, 0.000000e00, 4.000000e08],
            [1.187503e11, 2.100000e09, 0.000000e00, 8.000000e08],
            [1.187503e11, 3.000000e09, 0.000000e00, 1.000000e09],
            [1.187503e11, 5.000000e09, 0.000000e00, 2.000000e09],
            [2.432000e11, 2.500000e09, 0.000000e00, 3.000000e09],
            [2.432000e11, 2.500000e09, 0.000000e00, 3.000000e09],
            [3.251500e11, 1.500000e09, 0.000000e00, 1.600000e09],
            [3.251500e11, 3.500000e09, 0.000000e00, 2.400000e09],
            [3.251500e11, 9.500000e09, 0.000000e00, 3.000000e09],
            [4.240000e11, -1.000000e00, -1.000000e00, -1.000000e00],
            [4.240000e11, -1.000000e00, -1.000000e00, -1.000000e00],
            [4.240000e11, -1.000000e00, -1.000000e00, -1.000000e00],
            [4.240000e11, -1.000000e00, -1.000000e00, -1.000000e00],
            [4.480000e11, 1.400000e09, 0.000000e00, 1.200000e09],
            [4.480000e11, 3.000000e09, 0.000000e00, 2.000000e09],
            [4.480000e11, 7.200000e09, 0.000000e00, 3.000000e09],
            [6.640000e11, 4.200000e09, 0.000000e00, 5.000000e09],
            [6.640000e11, 4.200000e09, 0.000000e00, 5.000000e09],
            [8.740000e11, -1.000000e00, -1.000000e00, -1.000000e00],
            [8.740000e11, -1.000000e00, -1.000000e00, -1.000000e00],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation,
    [
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-H",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-V",
        "?",
        "?",
        "?",
        "?",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-V",
        "ISMAR-H",
        "ISMAR-V",
        "?",
        "?",
    ],
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, [])
# How many monochromatic frequencies to simulate the channel
ws.Touch(ws.met_mm_available_accuracies)
ws.Delete(ws.met_mm_available_accuracies)
# Number of frequencies for first accuracy (fast)
ws.ArrayOfIndexSet(
    ws.freq_number_tmp,
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1.0, -1.0, -1.0, -1.0, 1, 1, 1, 1, 1, -1.0, -1.0],
)
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for second accuracy (normal)
ws.ArrayOfIndexSet(
    ws.freq_number_tmp,
    [2, 2, 2, 2, 4, 2, 2, 4, 3, 5, -1.0, -1.0, -1.0, -1.0, 2, 3, 5, 10, 10, -1.0, -1.0],
)
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for third accuracy (high)
ws.ArrayOfIndexSet(
    ws.freq_number_tmp,
    [
        5,
        4,
        6,
        5,
        5,
        11,
        11,
        21,
        8,
        8,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        10,
        17,
        30,
        24,
        24,
        -1.0,
        -1.0,
    ],
)
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
# Number of frequencies for fourth accuracy (reference)
ws.ArrayOfIndexSet(
    ws.freq_number_tmp,
    [
        14,
        13,
        19,
        14,
        25,
        -1.0,
        -1.0,
        30,
        38,
        84,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        66,
        95,
        -1.0,
        92,
        92,
        -1.0,
        -1.0,
    ],
)
ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
ws.Delete(ws.freq_number_tmp)
ws.VectorSet(ws.freq_spacing_tmp, np.array([1.0e10, 1.0e09, 1.0e09, 1.0e09]))
ws.Extract(ws.met_mm_freq_number, ws.met_mm_available_accuracies, ws.met_mm_accuracy)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
