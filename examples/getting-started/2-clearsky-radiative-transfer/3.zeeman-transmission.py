import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
nf = 1001
ws.frequency_grid = np.linspace(-50e6, 50e6, nf) + line_f0

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

ws.spectral_radiance_unit = "1"
ws.spectral_radiance_space_agendaSet(option="Transmission")
ws.spectral_radiance_surface_agendaSet(option="Transmission")

# %% Core calculations

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.propagation_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyTransmission()

# %% Show results

plt.semilogy((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line")

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [3.41938759e-06, 7.73629681e-12, 1.48466689e-10, -4.44007326e-08],
            [1.68296779e-06, 6.36563685e-12, 1.12515548e-10, -2.69027168e-08],
            [6.84951553e-07, 5.12567520e-12, 7.99195217e-11, -1.43179215e-08],
            [1.98252387e-07, 4.04809327e-12, 5.08300685e-11, -6.05374848e-09],
            [2.54316071e-08, 3.27476721e-12, 2.52886453e-11, -1.49105559e-09],
            [3.22328343e-13, 1.81363558e-13, -2.57708670e-13, 7.96334058e-18],
            [2.54531426e-08, 3.28073742e-12, 2.53293549e-11, 1.49570802e-09],
            [1.98544497e-07, 4.06121370e-12, 5.09832834e-11, 6.08600656e-09],
            [6.86339664e-07, 5.14925742e-12, 8.02667321e-11, 1.44210546e-08],
            [1.68722782e-06, 6.40330448e-12, 1.13148512e-10, 2.71412968e-08],
            [3.42964480e-06, 7.79200472e-12, 1.49486896e-10, 4.48616704e-08],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
