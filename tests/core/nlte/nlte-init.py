import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["H2O"])

ws.abs_bands.readxml("nlte_lines.xml")

ws.atm_field["t"] = pyarts.arts.GriddedField3.fromxml("t.xml")

ws.atm_fieldInitializeNonLTE(normalization=0.75)

v = 0

for x in ws.atm_field.nlte:
    v+=ws.atm_field.nlte[x].data

assert np.allclose(v, 0.75, rtol=1e-6)
