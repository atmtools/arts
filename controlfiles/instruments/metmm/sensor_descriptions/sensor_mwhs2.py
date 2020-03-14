# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for MWS-2 simulations onboard FY-3C satellite.
# Lanuched in Sept. 2013
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
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Viewing angles
# Viewing angle information is taken from Table 15:
# https://directory.eoportal.org/web/eoportal/satellite-missions/f/fy-3#3Y4Y413eKram
# There are 48 different angles, corresponding to one side of the MWHS-2 scan.
ws.MatrixSet(
    ws.antenna_dlos,
    np.array(
        [
            [-53.35],
            [-52.25],
            [-51.15],
            [-50.05],
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
            [8.90000e10, 0.00000e00, 0.00000e00, 1.50000e09],
            [1.18750e11, 5.00000e09, 0.00000e00, 2.00000e09],
            [1.18750e11, 3.00000e09, 0.00000e00, 1.00000e09],
            [1.18750e11, 2.50000e09, 0.00000e00, 2.00000e08],
            [1.18750e11, 1.10000e09, 0.00000e00, 2.00000e08],
            [1.18750e11, 8.00000e08, 0.00000e00, 2.00000e08],
            [1.18750e11, 3.00000e08, 0.00000e00, 1.65000e08],
            [1.18750e11, 2.00000e08, 0.00000e00, 1.00000e08],
            [1.18750e11, 8.00000e07, 0.00000e00, 2.00000e07],
            [1.50000e11, 0.00000e00, 0.00000e00, 1.50000e09],
            [1.83311e11, 7.00000e09, 0.00000e00, 2.00000e09],
            [1.83311e11, 4.50000e09, 0.00000e00, 2.00000e09],
            [1.83311e11, 3.00000e09, 0.00000e00, 1.00000e09],
            [1.83311e11, 1.80000e09, 0.00000e00, 7.00000e08],
            [1.83311e11, 1.00000e09, 0.00000e00, 5.00000e08],
        ]
    ),
)
ws.ArrayOfStringSet(
    ws.met_mm_polarisation,
    [
        "AMSU-V",
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
    ],
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, [])
ws.ArrayOfIndexSet(
    ws.met_mm_freq_number, [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
)
ws.VectorSet(
    ws.freq_spacing_tmp,
    np.array(
        [
            1.0e10,
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
        ]
    ),
)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
