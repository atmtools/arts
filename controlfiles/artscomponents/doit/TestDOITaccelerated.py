# DEFINITIONS:  -*-sh-*-
#
# filename: TestDOITaccelerated.arts
#
# Demonstration of an accelerated DOIT scattering calculation
#
# Authors: Jakob Doerr
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.IndexSet(ws.stokes_dim, 4)
ws.execute_controlfile("artscomponents/doit/doit_setup.arts")
ws.execute_controlfile("artscomponents/doit/doit_setup_accelerated.arts")
ws.execute_controlfile("artscomponents/doit/doit_calc.arts")
ws.WriteXML(ws.output_file_format, ws.y, "", 0)
# ==================check==========================
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "artscomponents/doit/yREFERENCE_DOITaccelerated.xml")
ws.Compare(ws.y, ws.yREFERENCE, 1e-06)
# End of Main
