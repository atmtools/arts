import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

atm = pyarts.arts.AtmPoint()
atm.pressure = 101325.0
atm.temperature = 273.15
atm["O2"] = 0.21
atm["N2"] = 0.79
nd = atm.number_density("O2-66")

ws.abs_speciesSet(species=["O2-CIA-O2"])
ws.ReadCatalogData()
f = ws.abs_cia_data[0].data[0].grids[0]
x = ws.abs_cia_data.spectral_propmat(f=f, atm=atm) / nd

assert np.allclose(
    x[::300, 0],
    [
        0.00000000e00,
        0.00000000e00,
        8.28862173e-33,
        3.95506323e-32,
        8.94099309e-32,
        7.30562463e-31,
        2.05122205e-30,
        3.47065863e-30,
        2.61488493e-30,
        8.57954839e-31,
        1.08628874e-31,
        5.20303987e-32,
        5.68696897e-32,
        0.00000000e00,
    ],
    atol=1e-36,
)

ws.abs_speciesSet(species=["O2-CIA-N2"])
ws.ReadCatalogData()
f = ws.abs_cia_data[0].data[1].grids[0]
x = ws.abs_cia_data.spectral_propmat(f=f, atm=atm) / nd

assert np.allclose(
    x[::300, 0],
    [
        0.00000000e00,
        0.00000000e00,
        0.00000000e00,
        0.00000000e00,
        8.35363249e-32,
        4.48072042e-31,
        1.73662572e-30,
        5.76244958e-31,
        7.52829275e-32,
        6.40651315e-33,
        0.00000000e00,
        0.00000000e00,
        0.00000000e00,
        0.00000000e00,
        0.00000000e00,
    ],
    atol=1e-36,
)
