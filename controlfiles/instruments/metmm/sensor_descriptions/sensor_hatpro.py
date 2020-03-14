# DEFINITIONS:  -*-sh-*-
# ARTS sensor description for RPG HATPRO and HF simulations
# These are groundbased radiometers deployed at Summit Greenland
# ICECAPS campaign.
# Contact persons: Dave Turner (NOAA), Ralf Bennartz (Uni Wisconsin)
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
# Sensor characteristics based on the pdf-presentation
# ftp://ftp.etl.noaa.gov/psd3/arctic/summit/mwr/0_docs/Summit_Datagrams_MicroWaveRadiometer.pdf
# Viewing angles( !Caution ground-based instrument! )
# This instrument has only 1 viewing angle
ws.MatrixSet(ws.antenna_dlos, np.array([[180.0]]))
# Sensor response setup
# ---
ws.MatrixSet(
    ws.met_mm_backend,
    np.array(
        [
            [2.224e10, 0.000e00, 0.000e00, 2.300e08],
            [2.304e10, 0.000e00, 0.000e00, 2.300e08],
            [2.384e10, 0.000e00, 0.000e00, 2.300e08],
            [2.544e10, 0.000e00, 0.000e00, 2.300e08],
            [2.624e10, 0.000e00, 0.000e00, 2.300e08],
            [2.784e10, 0.000e00, 0.000e00, 2.300e08],
            [3.140e10, 0.000e00, 0.000e00, 2.300e08],
            [5.126e10, 0.000e00, 0.000e00, 1.820e08],
            [5.228e10, 0.000e00, 0.000e00, 1.790e08],
            [5.386e10, 0.000e00, 0.000e00, 1.880e08],
            [5.494e10, 0.000e00, 0.000e00, 1.700e08],
            [5.666e10, 0.000e00, 0.000e00, 7.040e08],
            [5.730e10, 0.000e00, 0.000e00, 9.270e08],
            [5.800e10, 0.000e00, 0.000e00, 1.854e09],
            [9.000e10, 0.000e00, 0.000e00, 2.000e09],
            [1.500e11, 0.000e00, 0.000e00, 2.000e09],
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
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
        "AMSU-V",
    ],
)
# Antenna is not supported for now
ws.VectorSet(ws.met_mm_antenna, [])
ws.ArrayOfIndexSet(
    ws.met_mm_freq_number,
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
)
ws.VectorSet(
    ws.freq_spacing_tmp,
    np.array(
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
        ]
    ),
)
ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
ws.Delete(ws.current_spacing)
