import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

atm = pyarts.arts.AtmPoint()
atm.pressure = 101325.0
atm.temperature = 273.15
atm["O2"] = 0.21
atm["H2O"] = 0.0
nd = atm.number_density("O2-66")

ws.abs_speciesSet(species=["O2-PWR98"])
ws.ReadCatalogData()
f = np.linspace(50e9, 70e9, 10)
x = ws.abs_predef_data.propagation_matrix(f=f, atm=atm) / nd

assert np.allclose(
    x[:, 0],
    [
        1.31538279e-29,
        3.13886678e-29,
        1.28291163e-28,
        4.16609430e-28,
        6.43587999e-28,
        7.09865626e-28,
        4.24233633e-28,
        1.13296982e-28,
        3.04226249e-29,
        1.43713835e-29,
    ],
    atol=1e-36,
)
