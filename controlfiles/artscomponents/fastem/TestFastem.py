import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.MatrixCreate("emissivity")
ws.MatrixCreate("reflectivity")
ws.VectorSet(ws.f_grid, np.array([1.80e11, 1.83e11]))
ws.VectorCreate("transmit")
ws.VectorSet(ws.transmit, np.array([0.9, 0.9]))
ws.FastemStandAlone(
    ws.emissivity,
    ws.reflectivity,
    ws.f_grid,
    283.0,
    180.0,
    0.1,
    3.0,
    0.0,
    ws.transmit,
    6,
)
ws.MatrixCreate("REFemissivity")
ws.MatrixCreate("REFreflectivity")
ws.ReadXML(ws.REFemissivity, "TestFastem.emissivityREFERENCE.xml")
ws.ReadXML(ws.REFreflectivity, "TestFastem.reflectivityREFERENCE.xml")
ws.Compare(ws.emissivity, ws.REFemissivity, 1e-06)
ws.Compare(ws.reflectivity, ws.REFreflectivity, 1e-06)
