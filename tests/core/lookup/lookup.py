import pyarts
import numpy as np
import matplotlib.pyplot as plt
from time import time


toa = 100e3

# %% Setup workspace

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["CO2-626"])

ws.ReadCatalogData()
for key in ws.absorption_bands:
    ws.absorption_bands[key].cutoff = "ByLine"
    ws.absorption_bands[key].cutoff_value = 750e9

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field["t"] = 295.0

ws.atmospheric_fieldRead(
    toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
data = pyarts.arts.Tensor3(ws.atmospheric_field["t"].data.data)

v = np.linspace(400, 2500, 1001)
ws.frequency_grid = pyarts.arts.convert.kaycm2freq(v)

t = time()
ws.absorption_lookup_tableFromProfiles(
    pressure_profile=np.array(ws.atmospheric_field["p"].data.flatten()),
    temperature_profile=np.array(ws.atmospheric_field["t"].data.flatten()),
    vmr_profiles={"CO2": np.array(ws.atmospheric_field["CO2"].data.flatten())},
    temperature_perturbation=np.linspace(-40, 40, 8),
)
print(time() - t, "s to train the LUT")

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)

# %% Checks and settings for LBL
for temperature_offset in np.linspace(-20, 20, 5):
    ws.atmospheric_field["t"].data.data = data + temperature_offset

    t = time()
    ws.propagation_matrix_agendaAuto(use_absorption_lookup_table=0)
    ws.spectral_radianceClearskyEmission()
    ws.spectral_radianceApplyUnitFromSpectralRadiance()
    lbl = ws.spectral_radiance[:, 0] * 1.0
    print(time() - t, "s to compute the LBL spectral radiance")

    t = time()
    ws.propagation_matrix_agendaAuto(
        f_interp_order=0,
        p_interp_order=5,
        t_interp_order=7,
        use_absorption_lookup_table=1,
    )
    ws.spectral_radianceClearskyEmission()
    ws.spectral_radianceApplyUnitFromSpectralRadiance()
    lut = ws.spectral_radiance[:, 0] * 1.0
    print(time() - t, "s to compute the LUT spectral radiance")

    plt.figure()
    plt.plot(v, lbl - lut, label="LBL - LUT")
    plt.legend()
    plt.xlabel("Frequency [cm$^{-1}$]")
    plt.ylabel("Radiance difference [K]")

    assert np.allclose(lbl, lut, atol=1e-3)

# %%
