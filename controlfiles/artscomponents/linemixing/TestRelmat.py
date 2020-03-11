import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.Wigner6Init()
ws.partition_functionsInitFromBuiltin()
ws.isotopologue_ratiosInitFromBuiltin()
ws.abs_speciesSet(species=["CO2-LM-626"])
ws.NumericSet(ws.rtp_pressure, 10.0)
ws.VectorNLinSpace(ws.abs_t, 50, 200.0, 350.0)
ws.ReadXML(ws.abs_lines_per_species, "testdata/abs_lines_per_band_relmat.xml")
ws.ReadXML(ws.band_identifiers, "testdata/co2band_relmat.xml")
ws.abs_lines_per_bandFromband_identifiers()
ws.SetRelaxationMatrixCalcType(type=[0])
ws.SetLineMixingCoefficinetsFromRelmat(order_of_linemixing=2)
# ReadXML(relmat_per_band, "test_relmat_per_band.xml")  # The reference test
ws.ArrayOfArrayOfMatrixCreate("test")
ws.ReadXML(ws.test, "relmat_per_bandREFERENCE.xml")
ws.CompareRelative(ws.test, ws.relmat_per_band, 1e-06)
