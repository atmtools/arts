import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.TessemNNReadAscii(ws.tessem_neth, "testdata/tessem_sav_net_H.txt")
ws.Print(ws.tessem_neth, 0)
ws.TessemNNReadAscii(ws.tessem_netv, "testdata/tessem_sav_net_V.txt")
ws.Print(ws.tessem_netv, 0)
ws.VectorCreate("tessem_in")
ws.VectorCreate("tessem_out")
ws.VectorCreate("tessem_ref")
ws.VectorSet(
    ws.tessem_in,
    np.array([1.0000000e10, 0.0000000e00, 0.0000000e00, 2.7314999e02, 3.0000000e-03]),
)
ws.VectorSet(ws.tessem_ref, np.array([0.395911]))
ws.TestTessem(ws.tessem_out, ws.tessem_neth, ws.tessem_in)
ws.Compare(ws.tessem_out, ws.tessem_ref, 1e-06)
ws.VectorSet(
    ws.tessem_in,
    np.array([1.0000000e10, 0.0000000e00, 0.0000000e00, 2.7314999e02, 3.0000000e-03]),
)
ws.VectorSet(ws.tessem_ref, np.array([0.374513]))
ws.TestTessem(ws.tessem_out, ws.tessem_netv, ws.tessem_in)
ws.Compare(ws.tessem_out, ws.tessem_ref, 1e-06)
