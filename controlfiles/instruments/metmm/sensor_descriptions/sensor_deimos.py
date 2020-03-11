# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for Deimos simulations
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
ws.Error("The MetMM file for DEIMOS is not yet finalised!!!")
# Viewing angles
# FIXME!!!
ws.MatrixSet(ws.antenna_dlos, array([], shape=(1, 0), dtype=float64))
# Sensor response setup
# ---
ws.MatrixSet(
    ws.met_mm_backend,
    array(
        [
            [2.38e10, 7.00e07, 0.00e00, 1.27e08],
            [2.38e10, 7.00e07, 0.00e00, 1.27e08],
            [5.01e10, 8.00e07, 0.00e00, 8.20e07],
            [5.01e10, 8.00e07, 0.00e00, 8.20e07],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation, ["ISMAR-V", "ISMAR-H", "ISMAR-V", "ISMAR-H"]
)
ws.VectorSet(ws.met_mm_antenna, array([], dtype=float64))
ws.ArrayOfIndexSet(ws.met_mm_freq_number, [12, 12, 12, 12])
ws.VectorSet(ws.freq_spacing_tmp, array([1.0e10, 1.0e09, 1.0e09, 1.0e09]))
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
