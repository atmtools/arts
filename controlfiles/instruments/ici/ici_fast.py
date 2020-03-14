# DEFINITIONS:  -*-sh-*-
# ARTS setup file for fast ICI simulations
#
# Uses ici.arts as basic setup (observation geometry, channel definitions).
#
# This is a version with a reduced frequency grid (only 1-2
# frequencies per channel).
# The frequencies have been derived with atmlab's simulated annealing
# functionality from the Garand data set (42 atmospheric cases).
# Optimization was run by Antonin Verlet-Banide on 140609 for channel
# accuracies of ?K.
#
# Accuracy for the different channels:
# Channel RMS_err[K]  Max_err[K]
#    01   0.0036
#    02   0.0546
#    03   0.0897
#    04   0.0597
#    05   0.0708
#    06   0.0997
#    07   0.00722
#    08   0.0804
#    09   0.0763
#    10   0.0836
#    11   0.0755
#
# 2014-06-09 Jana Mendrok, Antonin Verlet-Banide

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("instruments/ici/ici.arts")
# Replace f_grid and sensor_response by optimized ones.
#  Also, sensor_response_f is reset to effective channel frequencies resulting
#  from the chosen f_grid. This as a safety fix: In case yApplyUnit is
#  used, gives largely wrong answers for the default sensor_response_f setting,
#  particularly for bands with large IF.
# (other parameters remain like initialized above).
ws.ReadXML(ws.f_grid, "instruments/ici/ici.f_grid_fast.xml")
ws.ReadXML(ws.sensor_response, "instruments/ici/ici.sensor_response_fast.xml")
ws.ReadXML(ws.sensor_response_f, "instruments/ici/ici.sensor_response_f_fast.xml")
