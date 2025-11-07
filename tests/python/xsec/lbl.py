import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

atm = pyarts.arts.AtmPoint()
atm.pressure = 101325.0
atm.temperature = 273.15
atm["O2"] = 0.21
nd = atm.number_density("O2-66")

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
f = np.linspace(50e9, 70e9, 10)
x = ws.abs_bands.propagation_matrix(f=f, atm=atm) / nd

assert np.allclose(
    x[:, 0],
    [
        1.09249629e-29,
        2.79155620e-29,
        1.24134579e-28,
        4.10225824e-28,
        6.40578024e-28,
        6.96936469e-28,
        4.31517120e-28,
        1.19008426e-28,
        3.51438019e-29,
        1.81454497e-29,
    ],
    atol=1e-36,
)
