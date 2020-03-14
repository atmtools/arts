# DEFINITIONS:  -*-sh-*-
#
# Demonstration of a DOIT scattering calculation
#
# Author: Oliver Lemke
#
# This controlfile does the same calculation as TestDOIT,
# except that the sensor is positioned inside the cloudbox.
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.IndexSet(ws.stokes_dim, 4)
ws.execute_controlfile("artscomponents/doit/doit_setup.arts")
# Sensor altitude from earth surface
ws.nelemGet(ws.nelem, ws.vector_1)
# VectorSetConstant( vector_2, nelem, 95000.1 )
ws.VectorSetConstant(ws.vector_2, ws.nelem, 10000.0)
ws.Matrix1ColFromVector(ws.sensor_pos, ws.vector_2)
ws.execute_controlfile("artscomponents/doit/doit_calc.arts")
# WriteXML( in=y )
# ==================check==========================
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "artscomponents/doit/yREFERENCE_DOITsensorInsideCloudbox.xml")
ws.Compare(ws.y, ws.yREFERENCE, 1e-06)
# End of Main
