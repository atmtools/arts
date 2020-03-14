# DEFINITIONS:  -*-sh-*-
#
# filename: TestDOIT.arts
#
# Demonstration of a DOIT scattering calculation
#
# Author: Claudia Emde
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.IndexSet(ws.stokes_dim, 4)
ws.execute_controlfile("artscomponents/doit/doit_setup.arts")
ws.execute_controlfile("artscomponents/doit/doit_calc.arts")
ws.WriteXML(ws.output_file_format, ws.y, "", 0)
# ==================check==========================
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "artscomponents/doit/yREFERENCE_DOIT.xml")
ws.Compare(ws.y, ws.yREFERENCE, 1e-06)
# End of Main
