import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.VectorCreate("data_f_grid")
ws.VectorCreate("data_t_grid")
ws.VectorSet(ws.data_f_grid, np.array([2.3e11, 2.4e11]))
ws.VectorSet(ws.data_t_grid, np.array([220.0, 250.0, 270.0]))
ws.complex_refr_indexIceMatzler06(
    data_f_grid=ws.data_f_grid, data_T_grid=ws.data_t_grid
)
# WriteXML("ascii", complex_refr_index)
ws.GriddedField3Create("ref")
ws.ReadXML(ws.ref, "TestRefice.complex_refr_indexREFERENCE.xml")
ws.Compare(ws.complex_refr_index, ws.ref, 0.001)
