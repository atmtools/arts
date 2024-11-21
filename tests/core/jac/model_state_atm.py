import pyarts
import numpy as np


toa = 100e3
nalt = 20
nlat = 10
nlon = 5


# Remove value from x
def rem(x, v, start, length):
    end = start + length
    x[start:end] -= v
    return x


# Scale x by 1.1
def scl110(x, v, start, length):
    end = start + length
    x[start:end] = 1.1 * v
    return x


# A no-op call
def noop(x, v, start, length):
    return x


# Invoke some operation on all data
def ops(ws, op):
    x = ws.model_state_vector * 1.0

    for key in ws.jacobian_targets.atm:
        data = np.asarray(ws.atmospheric_field[key.type].data)
        x = op(x, data.flatten(), key.x_start, key.x_size)

    for key in ws.jacobian_targets.surf:
        data = np.asarray(ws.surface_field[key.type].data)
        x = op(x, data.flatten(), key.x_start, key.x_size)

    for key in ws.jacobian_targets.line:
        x = op(x, ws.absorption_bands[key.type], key.x_start, key.x_size)

    return x


ws = pyarts.Workspace()

ws.measurement_sensor = []

ws.jacobian_targetsInit()

ws.atmospheric_fieldInit(toa=toa)
ws.surface_fieldEarth()
ws.surface_field["t"] = 295
ws.surface_field["h"] = pyarts.arts.GriddedField2(
    data=2 + np.random.random((nlat, nlon)),
    name="Elevation",
    grids=(
        np.linspace(0, 90, nlat),
        np.linspace(0, 180, nlon),
    ),
    grid_names=["Latitude", "Longitude"],
)

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=100e9, fmax=120e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.absorption_bands = {bandkey: ws.absorption_bands[bandkey]}

# %% Temperature and pressure and VMR and ratios

ws.atmospheric_field[pyarts.arts.AtmKey.t] = 300.0
ws.jacobian_targetsAddTemperature()

ws.atmospheric_field[pyarts.arts.AtmKey.p] = pyarts.arts.GriddedField3(
    data=2 + np.random.random((nalt, nlat, nlon)),
    name="Pressure",
    grids=(
        np.linspace(0, toa, nalt),
        np.linspace(0, 90, nlat),
        np.linspace(0, 180, nlon),
    ),
    grid_names=["Altitude", "Latitude", "Longitude"],
)

ws.jacobian_targetsAddPressure()

spec = pyarts.arts.SpeciesEnum("H2O")
ws.atmospheric_field[spec] = 0.1
ws.jacobian_targetsAddSpeciesVMR(species=spec)

isot = pyarts.arts.SpeciesIsotope("H2O-161")
ws.jacobian_targetsAddSpeciesIsotopologueRatio(species=isot)

# %% Line center and Einstein and G0 and Y

# FIXME: Make a better interface!
# key = ws.absorption_bands[0].key
# ws.jacobian_targetsAddLineParameter(id=key, line_index=0, parameter="f0")
# ws.jacobian_targetsAddLineParameter(id=key, line_index=0, parameter="a")
# ws.jacobian_targetsAddLineParameter(
#     id=key, line_index=0, parameter="G0", species="O2"
# )
# ws.jacobian_targetsAddLineParameter(
#     id=key, line_index=0, parameter="Y", species="O2"
# )
# ws.jacobian_targetsAddLineParameter(
#     id=key, line_index=0, parameter="G0", species="AIR"
# )
# ws.jacobian_targetsAddLineParameter(
#     id=key, line_index=0, parameter="Y", species="AIR"
# )


# %% Surface temperature and altitude

ws.jacobian_targetsAddSurface(target="t")
ws.jacobian_targetsAddSurface(target="h")


# %% Finalize

ws.jacobian_targetsFinalize()

# %% Set the model_state_vector

ws.model_state_vectorFromData()

# %% Check equality

x = ops(ws, rem)
assert np.allclose(x, 0)

# %% Check scaling

x1 = ops(ws, scl110)
x2 = ops(ws, noop)
assert np.allclose(x1 / x2, 1.1)

# %% Check update

ws.UpdateModelStates(model_state_vector=x1)

ws.model_state_vectorFromData()
x2 = ops(ws, noop)

assert np.allclose(x1 / x2, 1.0)
