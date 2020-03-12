# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for MARSS simulations
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
            [-160.0],
            [-150.0],
            [-140.0],
            [-90.0],
            [-40.0],
            [-30.0],
            [-20.0],
            [-10.0],
            [0.0],
            [10.0],
            [20.0],
            [30.0],
            [40.0],
            [90.0],
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
            [8.89920e10, 1.07500e09, 0.00000e00, 6.50000e08],
            [1.57075e11, 2.60000e09, 0.00000e00, 2.60000e09],
            [1.83248e11, 9.75000e08, 0.00000e00, 4.50000e08],
            [1.83248e11, 3.00000e09, 0.00000e00, 1.00000e09],
            [1.83248e11, 7.00000e09, 0.00000e00, 2.00000e09],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation, ["MARSS-V", "MARSS-H", "MARSS-H", "MARSS-H", "MARSS-H"]
)
ws.VectorSet(ws.met_mm_antenna, [])
ws.ArrayOfIndexSet(ws.met_mm_freq_number, [12, 12, 12, 12, 12])
ws.VectorSet(ws.freq_spacing_tmp, np.array([1.0e10, 1.0e09, 1.0e09, 1.0e09]))
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
