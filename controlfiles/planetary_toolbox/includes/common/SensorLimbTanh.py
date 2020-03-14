#
# include file for limb sensor with
#  - specified tangent altitudes
#  - above-atmosphere platform altitude (here: fixed to 5000km, which should
#     always be above TOA)
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# sensor position and LOS
ws.nelemGet(ws.indtmp, ws.tanh)
ws.MatrixSetConstant(ws.sensor_pos, ws.indtmp, 1, 5000000.0)
ws.VectorZtanToZa1D(ws.tanh, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.tanh)
ws.Matrix1ColFromVector(ws.sensor_los, ws.tanh)
ws.WriteXML("ascii", ws.sensor_los)
