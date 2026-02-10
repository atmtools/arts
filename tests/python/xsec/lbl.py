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
x = ws.abs_bands.spectral_propmat(f=f, atm=atm) / nd

assert np.allclose(
    x[:, 0],
    [
        1.09244886e-29,
        2.79146590e-29,
        1.24133461e-28,
        4.10229972e-28,
        6.40580290e-28,
        6.96946879e-28,
        4.31520002e-28,
        1.19006194e-28,
        3.51422974e-29,
        1.81445133e-29,
    ],
    atol=1e-36,
)
