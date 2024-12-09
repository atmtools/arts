import pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False
toa = 100e3
lat = 0
lon = 0
NQuad = 40

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

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core Disort calculations

ws.propagation_matrix_scattering_spectral_agendaSet(option="FromSpeciesTRO")

ws.disort_settings_agendaSetup()

ws.disort_spectral_radiance_fieldProfile(
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

for manual in [False, True]:
    for EXT in [1e-4]:
        for SSA in [0.01]:
            for g in [0.9]:
                print(EXT, SSA, g)
                hspec = pyarts.arts.HenyeyGreensteinScatterer(
                    lambda f, atm: (
                        np.array([EXT, SSA])
                        if atm.pressure > 1e4 and atm.pressure < 9e4
                        else [0, 0]
                    ),
                    g,
                )

                if manual:
                    ws.scattering_species = [hspec]
                else:
                    def python_func(atm, f, index):
                        return hspec.get_bulk_scattering_properties_tro_spectral(
                            atm, f, index
                        )

                    ws.scattering_species = [
                        pyarts.arts.ScatteringGeneralSpectralTRO(python_func)
                    ]

                
                ws.disort_settings_agendaSetup(
                    scattering_setting="ScatteringSpecies",
                )

                ws.disort_spectral_radiance_fieldProfile(
                    longitude=lon,
                    latitude=lat,
                    disort_quadrature_dimension=NQuad,
                    disort_legendre_polynomial_dimension=1,
                    disort_fourier_mode_dimension=1,
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
                plt.title(f"EXT {EXT}; SSA {SSA}; g {g}; manual {manual}")
