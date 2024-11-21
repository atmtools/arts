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
ws.ray_path_observer_agendaSet(option="Geometric")
ws.spectral_radiance_observer_agendaSet(option="Emission")

# %% Set up a sensor with Gaussian standard deviation channel widths on individual frequency ranges

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.measurement_sensorSimpleGaussian(std=1e5, pos=pos, los=los, pol="RC")

# %% Core calculations

ws.measurement_vectorFromSensor()

# %% Show results

plt.plot((ws.frequency_grid - line_f0) / 1e6, ws.measurement_vector)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(
    f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line with Gaussian channels on individual grids"
)

# %% Test

assert np.allclose(
    ws.measurement_vector[::100],
    np.array(
        [
            227.85626271,
            230.93430882,
            234.89998118,
            240.50100578,
            250.05272664,
            209.9140708,
            249.51258095,
            240.21976958,
            234.71155853,
            230.79135297,
            227.73970752,
        ]
    ),
)
