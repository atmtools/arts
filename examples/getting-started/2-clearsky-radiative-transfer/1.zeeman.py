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
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()

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
            [2.27786385e02, -2.23642454e-04, -3.77012348e-04, 5.69263236e-02],
            [2.30865387e02, -3.48094282e-04, -5.82869673e-04, 7.04134321e-02],
            [2.34808042e02, -6.21526538e-04, -1.02887768e-03, 9.33834409e-02],
            [2.40362325e02, -1.41676037e-03, -2.28982262e-03, 1.40030136e-01],
            [2.49783921e02, -5.60527989e-03, -8.33754441e-03, 2.69732924e-01],
            [2.07245890e02, -4.30794421e00, -2.14467490e01, 1.24566607e-05],
            [2.49783389e02, -5.60754335e-03, -8.34040313e-03, -2.69809568e-01],
            [2.40361304e02, -1.41779162e-03, -2.29135603e-03, -1.40105962e-01],
            [2.34806573e02, -6.22187547e-04, -1.02991185e-03, -9.34587096e-02],
            [2.30863479e02, -3.48583496e-04, -5.83654375e-04, -7.04887260e-02],
            [2.27784046e02, -2.24033591e-04, -3.77649070e-04, -5.70015926e-02],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
