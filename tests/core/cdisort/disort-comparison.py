#!/usr/bin/env python3

import os
import sys

import pyarts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from datetime import datetime

# pyarts.data.parameters.datapath.append("/home/user/arts-xml-data-git")

if not pyarts.arts.globals.data.has_cdisort:
    print("cdisort is not available, skipping disort-comparison test")
    sys.exit(0)

# =============================================================================
# paths/constants
# =============================================================================

PLOT = False
toa = 100e3
lat = 0
lon = 0
N_quad = 40
N_freq = 50

# =============================================================================
# set up simulation
# =============================================================================

ws = pyarts.Workspace()

# Sampled frequency range
ws.frequency_grid = np.linspace(170, 194, N_freq) * 1e9
frequency_grid = np.array(ws.frequency_grid)

# Species and line absorption
ws.absorption_speciesSet(species=["H2O-PWR2022"])
ws.ReadCatalogData()

# Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# Grids and planet
ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=toa,
    basename="planets/Earth/afgl/tropical/",
)

# Checks and settings
ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")


# =============================================================================
# Disort calculations
# =============================================================================
ws.disort_settings_agendaSetup()

width_in_cm = 29.7 * 1.5
height_in_cm = width_in_cm / 1.421

fig, ax = plt.subplots(2, 4, sharex=True, layout="constrained")
fig.set_size_inches(width_in_cm / 2.54, h=height_in_cm / 2.54)

# cdisort calculation
dt = datetime.now()
ws.disort_spectral_radiance_fieldProfileCdisort(
    longitude=lon,
    latitude=lat,
    disort_quadrature_dimension=N_quad,
    disort_legendre_polynomial_dimension=1,
    disort_fourier_mode_dimension=1,
)
print("cdisort:   ", datetime.now() - dt)

cdisort_spectral_radiance_field = ws.disort_spectral_radiance_field * 1.0

# cppdisort calculation
dt = datetime.now()
ws.disort_spectral_radiance_fieldProfile(
    longitude=lon,
    latitude=lat,
    disort_quadrature_dimension=N_quad,
    disort_legendre_polynomial_dimension=1,
    disort_fourier_mode_dimension=1,
)
print("cppdisort: ", datetime.now() - dt)

radiance_field_absdiff = (
    cdisort_spectral_radiance_field - ws.disort_spectral_radiance_field
)

radiance_field_reldiff = (
    (cdisort_spectral_radiance_field - ws.disort_spectral_radiance_field)
    / ws.disort_spectral_radiance_field
    * 100
)

print(f"max reldiff: {np.max(np.abs(radiance_field_reldiff))}%")
imax = np.unravel_index(
    np.argmax(np.abs(radiance_field_reldiff)), radiance_field_reldiff.shape
)
print(
    f"max location: {ws.ray_path[imax[1]].pos[0] / 1000:.1f} km, "
    f"{ws.frequency_grid[imax[0]] / 1e9:.1f} GHz, "
    f"{ws.disort_quadrature_angles[imax[-1]]:.2f}°"
)


# =============================================================================
# Plotting
# =============================================================================
def format_plot(ax):
    ax.grid(which="both", linestyle=":", linewidth=0.25)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.axes.tick_params(direction="in", which="both")


def plot_field(fig, I, ax, title, reldiff=False):
    # =============================================================================
    # plot disort_cpp
    # =============================================================================
    cmap = cm["managua"]
    colors = cmap(np.linspace(0, 1, N_quad))
    angles = ws.disort_quadrature_angles[:]
    idx = np.argsort(angles)
    angles = angles[idx]
    I = I[:, :, idx] * 1.0

    # plot spectral radiance TOA and SFC
    for i in range(np.size(I, -1)):
        ax[0].plot(
            frequency_grid / 1e9,
            I[:, 0, i],
            label=f"{angles[i]:.2f}$^\\circ$",
            color=colors[i, :],
        )
        ax[1].plot(
            frequency_grid / 1e9,
            I[:, -1, i],
            label=f"{angles[i]:.2f}$^\\circ$",
            color=colors[i, :],
        )

    ylabel = "relative diff / %" if reldiff else "I / Wm$^{-2}$sr$^{-1}$Hz$^{-1}$"
    ax[1].legend(fontsize=6, ncols=2)
    ax[1].set_xlabel("frequency / GHz")
    ax[0].set_ylabel(ylabel)
    ax[0].set_title(f"TOA {title}")
    ax[1].set_ylabel(ylabel)
    ax[1].set_title(f"SFC {title}")

    format_plot(ax[0])
    format_plot(ax[1])


if PLOT:
    plot_field(fig, cdisort_spectral_radiance_field[:, :, 0, :], ax[0:2, 0], "cdisort")
    plot_field(
        fig, ws.disort_spectral_radiance_field[:, :, 0, :], ax[0:2, 1], "cppdisort"
    )

    plot_field(
        fig, radiance_field_absdiff[:, :, 0, :], ax[0:2, 2], "cdisort - cppdisort"
    )

    plot_field(
        fig,
        radiance_field_reldiff[:, :, 0, :],
        ax[0:2, 3],
        "(cdisort - cppdisort) / cppdisort",
        reldiff=True,
    )

    # fig.savefig(f"disort-toa{toa / 1e3:.0f}km.pdf")
    plt.show()

assert np.allclose(cdisort_spectral_radiance_field, ws.disort_spectral_radiance_field)
