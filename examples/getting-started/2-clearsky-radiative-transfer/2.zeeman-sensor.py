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
ws.measurement_sensorSimpleGaussian(fwhm=1e5, pos=pos, los=los, pol="RC")

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
            227.78791323,
            230.8653163,
            234.80800379,
            240.3623207,
            249.78397782,
            207.61855855,
            249.78344508,
            240.36129972,
            234.80653428,
            230.86340781,
            227.78557379,
        ]
    ),
)
