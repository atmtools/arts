import numpy as np
import pyarts3 as pyarts

# The line-by-line dispersion is the refractivity, n - 1, i.e., -c Im(pm) / (4 pi f)

ws = pyarts.Workspace()
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=130e9)
ws.WignerInit()

ws.atm_point["t"] = 290.0
ws.atm_point["p"] = 1013.25e2
ws.atm_point["O2"] = 0.2095
ws.ray_point.los = [180.0, 0.0]
ws.select_species = "AIR"

ws.jac_targetsInit()
ws.jac_targetsAddWindField(component="u")
ws.jac_targetsAddTemperature()

f = np.array([50e9, 58e9, 62e9, 70e9, 130e9])
ws.freq_grid = f


def spectral(df=0.0):
    ws.freq_grid = f + df
    ws.spectral_propmatInit()
    ws.spectral_propmatAddVoigtLTE()
    ws.freq_grid = f
    return ws.spectral_dispersion * 1.0, ws.spectral_dispersion_jac * 1.0


d, dd = spectral()

# The O2 60 GHz band at sea level gives refractivity of order 1e-6 that
# is positive below and negative above the band
assert np.all(np.abs(d) < 1e-5), d
assert np.all(np.abs(d) > 1e-8), d
assert d[1] > 0 and d[2] < 0, d

# Same result from the memory intensive method, also when called twice
ws.spectral_propmatInit()
ws.spectral_propmatMemoryIntensiveAddVoigtLTE()
ws.spectral_propmatInit()
ws.spectral_propmatMemoryIntensiveAddVoigtLTE()
assert np.allclose(ws.spectral_dispersion, d, rtol=1e-6, atol=0), (ws.spectral_dispersion, d)
assert np.allclose(ws.spectral_dispersion_jac, dd, rtol=1e-6, atol=0)

# Line-by-line wind derivatives are frequency derivatives, including that of 1 / f.
# The line shape derivative is numerical, so allow 1e-3.  Ignoring 1 / f gives > 3e-2.
df = 1e5
d2, _ = spectral(df)
assert np.allclose(dd[0], (d2 - d) / df, rtol=1e-3, atol=0), (dd[0], (d2 - d) / df)

# Single frequency points agree and accumulate onto what is already there
ws.freq_wind_shift_jac = [0.0, 0.0, 0.0]
for i, fi in enumerate(f):
    ws.freq = fi
    ws.single_propmatInit()
    ws.single_propmatAddVoigtLTE()
    assert np.isclose(ws.single_dispersion, d[i], rtol=1e-12, atol=0)
    assert np.allclose(ws.single_dispersion_jac, dd[:, i], rtol=1e-12, atol=0)

    ws.single_propmatInit()
    ws.single_dispersionAddGasMicrowavesEarth()
    earth = ws.single_dispersion * 1.0
    ws.single_propmatAddVoigtLTE()
    assert np.isclose(ws.single_dispersion, earth + d[i], rtol=1e-12, atol=0)
