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

plt.plot((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance + 0)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line")

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [2.27786279e02, -2.23644583e-04, -3.77015913e-04, 5.69269522e-02],
            [2.30865312e02, -3.48096654e-04, -5.82873607e-04, 7.04140485e-02],
            [2.34807995e02, -6.21529170e-04, -1.02888198e-03, 9.33839855e-02],
            [2.40362302e02, -1.41676334e-03, -2.28982732e-03, 1.40030558e-01],
            [2.49783915e02, -5.60528357e-03, -8.33754971e-03, 2.69733178e-01],
            [2.07245890e02, -4.30794421e00, -2.14467490e01, 1.24566607e-05],
            [2.49783382e02, -5.60754703e-03, -8.34040843e-03, -2.69809821e-01],
            [2.40361281e02, -1.41779460e-03, -2.29136074e-03, -1.40106384e-01],
            [2.34806526e02, -6.22190182e-04, -1.02991615e-03, -9.34592546e-02],
            [2.30863403e02, -3.48585873e-04, -5.83658317e-04, -7.04893428e-02],
            [2.27783939e02, -2.24035725e-04, -3.77652642e-04, -5.70022211e-02],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
