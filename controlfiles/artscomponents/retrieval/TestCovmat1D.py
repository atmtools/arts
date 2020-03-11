import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.VectorCreate("grid_1")
ws.VectorCreate("sigma_1")
ws.VectorCreate("cls_1")
ws.VectorLinSpace(ws.grid_1, 0.0, 30000.0, 1500.0)
ws.VectorSetConstant(out=ws.sigma_1, nelem=1, value=1.0)
ws.VectorLinSpace(ws.cls_1, 1000.0, 2000.0, 50.0)
ws.VectorCreate("grid_2")
ws.VectorCreate("sigma_2")
ws.VectorCreate("cls_2")
ws.VectorLinSpace(ws.grid_2, 0.0, 30000.0, 3000.0)
ws.VectorSetConstant(out=ws.sigma_2, nelem=1, value=2.0)
ws.VectorLinSpace(ws.cls_2, 1000.0, 2000.0, 100.0)
ws.SparseCreate("covmat_reference")
ws.SparseCreate("covmat")
ws.ReadXML(ws.covmat_reference, "artscomponents/retrieval/covmat1D_lin.xml")
ws.covmat1D(
    ws.covmat,
    ws.grid_1,
    ws.grid_2,
    ws.sigma_1,
    ws.sigma_2,
    ws.cls_1,
    ws.cls_2,
    0.001,
    "lin",
)
ws.Compare(
    ws.covmat_reference, ws.covmat, 1e-06, "Result of covmat1D deviates from reference."
)
ws.ReadXML(ws.covmat_reference, "artscomponents/retrieval/covmat1D_exp.xml")
ws.covmat1D(
    ws.covmat,
    ws.grid_1,
    ws.grid_2,
    ws.sigma_1,
    ws.sigma_2,
    ws.cls_1,
    ws.cls_2,
    0.001,
    "exp",
)
ws.Compare(
    ws.covmat_reference, ws.covmat, 1e-06, "Result of covmat1D deviates from reference."
)
ws.ReadXML(ws.covmat_reference, "artscomponents/retrieval/covmat1D_gau.xml")
ws.covmat1D(
    ws.covmat,
    ws.grid_1,
    ws.grid_2,
    ws.sigma_1,
    ws.sigma_2,
    ws.cls_1,
    ws.cls_2,
    0.001,
    "gau",
)
ws.Compare(
    ws.covmat_reference, ws.covmat, 1e-06, "Result of covmat1D deviates from reference."
)
