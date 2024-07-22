import pyarts
import numpy as np


def rem(x, v, start, length):
    end = start + length
    x[start:end] -= v
    return x, end


toa = 100e3
nalt = 20
nlat = 10
nlon = 5


ws = pyarts.Workspace()

ws.jacobian_targetsInit()

ws.atmospheric_fieldInit(toa=toa)
ws.surface_fieldEarth()

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()

# %% Temperature and pressure

ws.atmospheric_field[pyarts.arts.AtmKey.t] = 300.0
ws.jacobian_targetsAddTemperature()

orig_p = np.random.random((nalt, nlat, nlon))
ws.atmospheric_field[pyarts.arts.AtmKey.p] = pyarts.arts.GriddedField3(
    data=orig_p,
    name="Pressure",
    grids=(
        np.linspace(0, toa, nalt),
        np.linspace(0, 90, nlat),
        np.linspace(0, 180, nlon),
    ),
    grid_names=["Altitude", "Latitude", "Longitude"],
)

ws.jacobian_targetsAddPressure()


# %% VMR

spec = pyarts.arts.SpeciesEnum("H2O")
ws.atmospheric_field[spec] = 0.1
ws.jacobian_targetsAddSpeciesVMR(species="H2O")


# %% Isotopologue

isot = pyarts.arts.SpeciesIsotope("H2O-161")
ws.jacobian_targetsAddSpeciesIsotopologueRatio(species="H2O-161")


# %% Finalize

ws.jacobian_targetsFinalize()

# %% Set the model_state_vector

ws.model_state_vectorFromData()

x = ws.model_state_vector * 1.0

n = 0

x, n = rem(x, 300.0, n, 1)
x, n = rem(x, orig_p.flatten(), n, nalt * nlat * nlon)
x, n = rem(x, ws.atmospheric_field[spec].data, n, 1)
x, n = rem(x, ws.atmospheric_field[isot].data, n, 1)

assert np.allclose(x, 0)
