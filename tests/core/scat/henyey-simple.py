import pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False
toa = 100e3
lat = 0
lon = 0
NQuad = 40

if pyarts.arts.globals.has_sht:
    ws = pyarts.Workspace()

    line_f0 = 118750348044.712
    ws.frequency_grid = [line_f0]
    ws.frequency_grid = np.linspace(-20e9, 2e6, 101) + line_f0

    # %% Species and line absorption

    ws.absorption_speciesSet(species=["O2-66"])
    ws.ReadCatalogData()

    # %% Use the automatic agenda setter for propagation matrix calculations
    ws.propagation_matrix_agendaAuto()

    # %% Grids and planet

    ws.surface_fieldSetPlanetEllipsoid(option="Earth")
    ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
    ws.atmospheric_fieldRead(
        toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
    )
    ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")
    # ws.atmospheric_field[pyarts.arts.AtmKey.t] = 300.0

    # %% Checks and settings

    ws.spectral_radiance_unit = "Tb"
    ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
    ws.spectral_radiance_surface_agendaSet(option="Blackbody")

    # %% Core Disort calculations

    ws.disort_spectral_radiance_fieldClearsky(
        longitude=lon,
        latitude=lat,
        disort_quadrature_dimension=NQuad,
        disort_legendre_polynomial_dimension=1,
        disort_fourier_mode_dimension=1,
    )

    # %% Equivalent ARTS calculations

    ws.ray_pathGeometricDownlooking(
        latitude=lat,
        longitude=lon,
        max_step=1000.0,
    )
    ws.spectral_radianceClearskyEmission()

    # %% Plot results
    f = (ws.frequency_grid - line_f0) / 1e9

    plt.figure()
    plt.semilogy(f, ws.spectral_radiance[:, 0], "k--", lw=3)
    plt.semilogy(
        f,
        ws.disort_spectral_radiance_field[:, 0, 0, 0],
        "g:",
        lw=3,
    )
    plt.semilogy(
        f,
        ws.disort_spectral_radiance_field[:, 0, 0, (NQuad // 2) - 1],
        "m:",
        lw=3,
    )
    plt.ylabel("Spectral radiance [W sr$^{-1}$ m$^{-2}$ Hz$^{-1}$]")
    plt.xlabel("Dirac frequency [GHz]")
    plt.title("Downlooking")

    ws.propagation_matrix_scattering_totally_random_orientation_spectral_agendaSet(
        option="FromSpecies"
    )

    for EXT in np.logspace(-9, -1, 3):
        for SSA in np.logspace(-5, -1, 5):
            for g in [-1.0, -0.8, 0.0, 0.2, 1.0]:
                print(EXT, SSA, g)
                ws.scattering_species = [
                    pyarts.arts.HenyeyGreensteinScatterer(
                        lambda f, atm: (
                            np.array([EXT, SSA])
                            if atm.pressure > 1e4 and atm.pressure < 9e4
                            else [0, 0]
                        ),
                        g,
                    )
                ]

                ws.disort_spectral_radiance_fieldScatteringSpecies(
                    longitude=lon,
                    latitude=lat,
                    disort_quadrature_dimension=NQuad,
                    disort_legendre_polynomial_dimension=15,
                    disort_fourier_mode_dimension=15,
                )

                plt.figure()
                plt.semilogy(f, ws.spectral_radiance[:, 0], "k--", lw=3)
                plt.semilogy(
                    f,
                    ws.disort_spectral_radiance_field[:, 0, 0, 0],
                    "g:",
                    lw=3,
                )
                plt.semilogy(
                    f,
                    ws.disort_spectral_radiance_field[:, 0, 0, (NQuad // 2) - 1],
                    "m:",
                    lw=3,
                )
                plt.ylabel("Spectral radiance [W sr$^{-1}$ m$^{-2}$ Hz$^{-1}$]")
                plt.xlabel("Dirac frequency [GHz]")
                plt.title(f"EXT {EXT}; SSA {SSA}; g {g}")
