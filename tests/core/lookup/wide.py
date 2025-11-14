import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from time import time


toa = 100e3

# %% Setup workspace

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["CO2-626", "H2O-161"])

ws.ReadCatalogData()
for key in ws.abs_bands:
    ws.abs_bands[key].cutoff = "ByLine"
    ws.abs_bands[key].cutoff_value = 750e9

ws.abs_bands.keep_hitran_s(70)

ws.surface_fieldPlanet(option="Earth")
ws.surface_field["t"] = 295.0

ws.atm_fieldRead(
    toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
tdata = pyarts.arts.Tensor3(ws.atm_field["t"].data.data)
wdata = pyarts.arts.Tensor3(ws.atm_field["H2O"].data.data)

v = np.linspace(400, 2500, 101)
ws.freq_grid = pyarts.arts.convert.kaycm2freq(v)

t = time()
ws.abs_lookup_dataSimpleWide(water_affected_species=["H2O"],
                                     pressure_range=[1e-2, 1100e2])
print(round(1000 * (time() - t)), "ms to train the LUT")

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)

# %% Checks and settings for LBL
for water_ratio in [5e-1, 5]:
    for temperature_offset in np.linspace(-20, 20, 5):
        print(water_ratio, temperature_offset)
        ws.atm_field["t"].data.data = tdata + temperature_offset
        ws.atm_field["H2O"].data.data = wdata * water_ratio

        t = time()
        ws.propagation_matrix_agendaAuto(use_abs_lookup_data=0)
        ws.spectral_radianceClearskyEmission()
        ws.spectral_radianceApplyUnitFromSpectralRadiance()
        lbl = ws.spectral_radiance[:, 0] * 1.0
        print(round(1000 * (time() - t)), "ms to compute the LBL spectral radiance")

        t = time()
        ws.propagation_matrix_agendaAuto(
            f_interp_order=0,
            p_interp_order=5,
            t_interp_order=7,
            use_abs_lookup_data=1,
        )
        ws.spectral_radianceClearskyEmission()
        ws.spectral_radianceApplyUnitFromSpectralRadiance()
        lut = ws.spectral_radiance[:, 0] * 1.0
        print(round(1000 * (time() - t)), "ms to compute the LUT spectral radiance")

        plt.figure()
        plt.plot(v, lbl - lut, label="LBL - LUT")
        plt.legend()
        plt.xlabel("Frequency [cm$^{-1}$]")
        plt.ylabel("Radiance difference [K]")
        plt.title(
            f"Water ratio: {100*water_ratio}%, Temperature offset: {temperature_offset} K"
        )

        assert np.allclose(lbl, lut, atol=1e-3)
