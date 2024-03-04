import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-50e6, 50e6, 1001) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=40e9, fmax=120e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core calculations

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.propagation_pathGeometric(pos=pos, los=los, max_step=1000.)
ws.spectral_radianceStandardEmission()

# %% Show results

plt.plot((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line")

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [2.27786383e02, -1.99350269e-05, -4.38513166e-04, 5.61378144e-02],
            [2.30865383e02, -3.27976598e-05, -6.79061901e-04, 6.94379197e-02],
            [2.34808036e02, -6.39244866e-05, -1.20203107e-03, 9.20893766e-02],
            [2.40362309e02, -1.71072539e-04, -2.69107006e-03, 1.38088979e-01],
            [2.49783872e02, -1.00989029e-03, -1.00100939e-02, 2.65991699e-01],
            [2.07223843e02, 6.27903175e00, -2.10085556e01, 1.22487429e-05],
            [2.49783340e02, -1.01053458e-03, -1.00136877e-02, -2.66067280e-01],
            [2.40361289e02, -1.71259219e-04, -2.69291112e-03, -1.38163753e-01],
            [2.34806566e02, -6.40205116e-05, -1.20325640e-03, -9.21636026e-02],
            [2.30863475e02, -3.28598660e-05, -6.79985804e-04, -6.95121708e-02],
            [2.27784043e02, -1.99804947e-05, -4.39260082e-04, -5.62120412e-02],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
