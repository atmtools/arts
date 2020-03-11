# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for SAPHIR L1A2 simulations
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
# Sensor characteristics based
# Viewing angles
# There are 65 different angles, corresponding to one side of the SAPHIR Scan.
# Viewing angles definition from Table 3.2-2 Scan angle coverage. Distance
# between pixels is always 0.66:
# https://cnes.fr/fr/media/20130117level-1productdefed3rev4pdf

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.MatrixSet(
    ws.antenna_dlos,
    array(
        [
            [-42.96],
            [-42.3],
            [-41.64],
            [-40.98],
            [-40.31],
            [-39.65],
            [-38.99],
            [-38.33],
            [-37.67],
            [-37.01],
            [-36.34],
            [-35.68],
            [-35.02],
            [-34.36],
            [-33.7],
            [-33.04],
            [-32.37],
            [-31.71],
            [-31.05],
            [-30.39],
            [-29.73],
            [-29.07],
            [-28.41],
            [-27.74],
            [-27.08],
            [-26.42],
            [-25.76],
            [-25.1],
            [-24.44],
            [-23.77],
            [-23.11],
            [-22.45],
            [-27.79],
            [-21.13],
            [-20.47],
            [-19.8],
            [-17.14],
            [-18.48],
            [-17.82],
            [-17.16],
            [-16.5],
            [-15.83],
            [-15.17],
            [-14.51],
            [-13.85],
            [-13.19],
            [-12.53],
            [-11.87],
            [-11.2],
            [-10.54],
            [-9.88],
            [-9.22],
            [-8.56],
            [-7.9],
            [-7.23],
            [-6.57],
            [-5.91],
            [-5.25],
            [-4.59],
            [-3.93],
            [-3.26],
            [-2.6],
            [-1.94],
            [-1.27],
            [-0.61],
        ]
    ),
)
ws.MatrixSet(
    ws.met_mm_backend,
    array(
        [
            [1.8331e11, 2.0000e08, 0.0000e00, 2.0000e08],
            [1.8331e11, 1.1000e09, 0.0000e00, 3.5000e08],
            [1.8331e11, 2.8000e09, 0.0000e00, 5.0000e08],
            [1.8331e11, 4.2000e09, 0.0000e00, 7.0000e08],
            [1.8331e11, 6.8000e09, 0.0000e00, 1.2000e09],
            [1.8331e11, 1.1000e10, 0.0000e00, 2.0000e09],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation, ["AMSU-V", "AMSU-V", "AMSU-V", "AMSU-V", "AMSU-V", "AMSU-V"]
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, array([], dtype=float64))
ws.ArrayOfIndexSet(ws.met_mm_freq_number, [12, 12, 12, 12, 12, 12])
ws.VectorSet(
    ws.freq_spacing_tmp, array([1.0e09, 1.0e09, 1.0e09, 1.0e09, 1.0e09, 1.0e09])
)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
