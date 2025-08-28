import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["H2O"])

ws.absorption_bands.readxml("nlte_lines.xml")

ws.atmospheric_field["t"] = pyarts.arts.GriddedField3.fromxml("t.xml")

ws.atmospheric_fieldInitializeNonLTE(normalization=0.75)

v = 0

for x in ws.atmospheric_field.nlte:
    v+=ws.atmospheric_field.nlte[x].data

assert np.allclose(v, 0.75, rtol=1e-6)
