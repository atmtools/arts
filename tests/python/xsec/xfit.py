import pyarts
import numpy as np

ws = pyarts.Workspace()

atm = pyarts.arts.AtmPoint()
atm.pressure = 101325.0
atm.temperature = 273.15
atm["O3"] = 1e-8
nd = atm.number_density("O3-666")

ws.absorption_speciesSet(species=["O3-XFIT"])
ws.ReadCatalogData()
f = ws.absorption_xsec_fit_data[0].fitcoeffs[0].grids[0]
x = ws.absorption_xsec_fit_data.propagation_matrix(f=f, atm=atm) / nd

assert np.allclose(
    x[::60000, 0],
    [
        5.48830541e-29,
        7.01126943e-26,
        2.25566466e-25,
        3.67157516e-27,
        5.42679027e-26,
        7.56822053e-23,
        1.11608430e-21,
        3.90151892e-22,
    ],
    atol=1e-36,
)
