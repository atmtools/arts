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
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
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

plt.plot((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance + 0)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line")

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [2.27784834e02, -2.23635319e-04, -3.77001634e-04, 5.69266632e-02],
            [2.30863853e02, -3.48081919e-04, -5.82851478e-04, 7.04138002e-02],
            [2.34806521e02, -6.21501805e-04, -1.02884266e-03, 9.33836080e-02],
            [2.40360801e02, -1.41669913e-03, -2.28974401e-03, 1.40029968e-01],
            [2.49782407e02, -5.60506178e-03, -8.33739163e-03, 2.69735122e-01],
            [2.07248447e02, -4.30842579e00, -2.14475851e01, 1.25057382e-05],
            [2.49781840e02, -5.60734272e-03, -8.34027514e-03, -2.69812575e-01],
            [2.40359708e02, -1.41773814e-03, -2.29128954e-03, -1.40106546e-01],
            [2.34804943e02, -6.22167771e-04, -1.02988483e-03, -9.34596090e-02],
            [2.30861802e02, -3.48574780e-04, -5.83642167e-04, -7.04898114e-02],
            [2.27782317e02, -2.24029336e-04, -3.77643130e-04, -5.70026253e-02],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
