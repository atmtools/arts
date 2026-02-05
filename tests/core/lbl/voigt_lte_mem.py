import pyarts3 as pyarts
import numpy as np
import os

ws = pyarts.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-5e6, 5e6, 1001) + line_f0

# %% Species and line absorption
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=1200e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Set the atmosphere
ws.atm_point["t"] = 300.0
ws.atm_point["p"] = 1013.25e-2
ws.atm_point["O2"] = 0.2095
ws.atm_point.mag = [10e-6, 50e-6, 1e-6]
ws.ray_point.los = [180.0, 0.0]

# %% Set the jacobian
ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
ws.jac_targetsAddSpeciesVMR(species="O2")
ws.jac_targetsAddMagneticField(component="u")
ws.jac_targetsAddMagneticField(component="v")
ws.jac_targetsAddMagneticField(component="w")
ws.jac_targetsAddWindField(component="u")
ws.jac_targetsAddWindField(component="v")
ws.jac_targetsAddWindField(component="w")
x = ["Temperature", "VMR_O2", "Mag_u", "Mag_v", "Mag_w", "Wind_u", "Wind_v", "Wind_w"]

# %% Calculations

ws.spectral_propmatInit()
ws.spectral_propmatMemoryIntensiveAddVoigtLTE()

d = ws.spectral_dispersion * 1.0
pm = ws.spectral_propmat * 1.0
dpm = ws.spectral_propmat_jac * 1.0
dd = ws.spectral_dispersion_jac * 1.0
f = ws.freq_grid * 1.0
ws.freq_grid = f + 10
ws.spectral_propmatInit()
ws.spectral_propmatMemoryIntensiveAddVoigtLTE()
d2 = ws.spectral_dispersion * 1.0
pm2 = ws.spectral_propmat * 1.0
dpm2 = (pm2 - pm) / 10
ws.freq_grid = f

# %% Standard code

ws.spectral_propmat = []
ws.spectral_propmat_jac = [[]]
ws.spectral_propmatInit()
ws.spectral_propmatAddLines()

if "ARTS_HEADLESS" not in os.environ:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(8, 12))
    ax[0].plot(ws.freq_grid/1e9, d)
    ax[1].plot(ws.freq_grid/1e9, pm)
    ax[2].plot(ws.freq_grid/1e9, ws.spectral_propmat)

    for i in range(dd.shape[0]):
        fig, ax = plt.subplots(
            3 + ("Wind" in x[i]), 1, figsize=(8, 12 + ("Wind" in x[i])*4))
        ax[0].plot(ws.freq_grid/1e9, dd[i])
        ax[1].plot(ws.freq_grid/1e9, dpm[i])
        if "Wind" in x[i]:
            ax[3].plot(ws.freq_grid/1e9, dpm2)
        ax[2].plot(ws.freq_grid/1e9, ws.spectral_propmat_jac[i])
        ax[0].set_title(f"Jacobian for {x[i]}")

assert np.allclose(pm, ws.spectral_propmat)
# assert np.allclose(dpm, ws.spectral_propmat_jac)  # Disabled due zero-crossings

