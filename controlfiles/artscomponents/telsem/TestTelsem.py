import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.telsem_atlasesReadAscii(directory="/home/simon/src/telsem2/Atlas")
ws.TelsemAtlasCreate("telsem_atlas")
ws.Extract(ws.telsem_atlas, ws.telsem_atlases, 0)
ws.MatrixCreate("emis")
ws.MatrixCreate("REFemis")
ws.MatrixSet(ws.REFemis, np.array([[0.95400876, 0.95400876]]))
ws.VectorCreate("fs")
ws.VectorSet(ws.fs, np.array([3.0e11]))
ws.telsemStandalone(ws.emis, -30.0, 302.0, 15.0, ws.fs, ws.telsem_atlas)
# Compare(REFemis,  emis, 1e-6)
ws.VectorCreate("emis_vector")
ws.telsemAtlasLookup(ws.emis_vector, 89.75, 76.0, ws.telsem_atlas)
ws.Print(ws.emis_vector, 0)
