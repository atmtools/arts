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
ws.spectral_radiance_observer_agendaSet(option="EmissionUnits")

# %% Set up a sensor with Gaussian channel widths on individual frequency ranges

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.measurement_vector_sensorGaussian(
    f0_fwhm_df=[[f, 1e5, 1e4] for f in ws.frequency_grid],
    pos=pos,
    los=los,
    pol="RC",
)

# %% Core calculations

result = pyarts.arts.Vector()
result_jac = pyarts.arts.Matrix()
ws.measurement_vectorFromSensor(result, result_jac)

# %% Show results

plt.plot((ws.frequency_grid - line_f0) / 1e6, result)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(
    f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line with Gaussian channels on individual grids"
)


# %% Test

assert np.allclose(
    result[::100],
    np.array(
        [
            227.84333794,
            230.93584277,
            234.90150277,
            240.50252995,
            250.05423728,
            209.93307662,
            249.51413112,
            240.2213661,
            234.713189,
            230.79303123,
            227.7270699,
        ]
    ),
)


# %% Set up a sensor with a fixed grid size

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.measurement_vector_sensorGaussianFrequencyGrid(
    f0_fwhm=[[f, 1e5] for f in ws.frequency_grid],
    pos=pos,
    los=los,
    pol="RC",
)

# %% Core calculations

result = pyarts.arts.Vector()
result_jac = pyarts.arts.Matrix()
ws.measurement_vectorFromSensor(result, result_jac)

# %% Show results

plt.plot((ws.frequency_grid - line_f0) / 1e6, result)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(
    f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line with Gaussian channels on a single grid"
)


# %% Test

assert np.allclose(
    result[::100],
    np.array(
        [
            227.85779479,
            230.93584267,
            234.9015026,
            240.50252956,
            250.05423597,
            209.91116972,
            249.51412988,
            240.22136572,
            234.71318884,
            230.79303114,
            227.74141792,
        ]
    ),
)
