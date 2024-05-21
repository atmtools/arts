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
            227.84254707,
            230.93486349,
            234.90020182,
            240.50057327,
            250.05044696,
            209.877013,
            249.51782475,
            240.22329289,
            234.71447726,
            230.79400404,
            227.7278571,
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
            227.85700315,
            230.93486339,
            234.90020165,
            240.50057288,
            250.05044565,
            209.854878,
            249.51782351,
            240.22329251,
            234.71447709,
            230.79400394,
            227.74220588,
        ]
    ),
)
