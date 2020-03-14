# DEFINITIONS:  -*-sh-*-
# Prepares data container commonly used by the metmm sensor setting system.
#
# It has to be called before any sensor_*.arts using the metmm system.

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.ArrayOfArrayOfIndexCreate("met_mm_available_accuracies")
ws.ArrayOfIndexCreate("freq_number_tmp")
ws.VectorCreate("freq_spacing_tmp")
ws.ArrayOfIndexCreate("met_mm_freq_number")
ws.NumericCreate("current_spacing")
ws.VectorCreate("met_mm_freq_spacing")
ws.IndexCreate("met_mm_nchannels")
