"""
Scattering calculations with DISORT

This script demonstrates how to simulate brightness temperatures of cloudy
skies using the DISORT scattering solver.
"""

import numpy as np
import pyarts3 as pyarts
from pyarts3.arts import (
    AtmPoint,
    GriddedField3,
    MGDSingleMoment,
    ParticleHabit,
    ScatteringHabit,
)
from pyarts3.xml import load

toa = 100e3
lat = 0
lon = 0
NQuad = 40
ws = pyarts.Workspace()

pyarts.arts.globals.omp_set_num_threads(1)

ws.frequency_grid = [31.5e9, 165e9, 666e9]

# %% Species and line absorption
ws.absorption_speciesSet(species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"])
ws.ReadCatalogData()
ws.propagation_matrix_agendaAuto()

# %% Defining scattering species
#
# ### Loading scattering data
#
# We load scattering and scattering meta data from the ARTS .xml format. The
# scattering and scattering meta data arrays contain scattering data for rain
# particles (index 0) and ice particle (index 1).
scat_data_raw = load("scat_data.xml")
scat_data_meta = load("scat_meta.xml")

t_grid = scat_data_raw[0][0].T_grid
f_grid = scat_data_raw[0][0].f_grid

# Next we define a rain habit that holds the scattering data for liquid cloud and rain drops of different sizes.
#
# > **Note**: We transform the rain habit to the spectral TRO representation that is expected by DISORT already here.
rain_habit = ParticleHabit.from_legacy_tro(scat_data_raw[0], scat_data_meta[0])
rain_habit = rain_habit.to_tro_spectral(t_grid, f_grid, 39)


# The particle habit, which holds the scattering data, needs to be combined
# with a PSD so that it can be used as a scattering species in ARTS. The PSD
# needs to be linked to an atmospheric field through a particule property. The
# particulate property represents the moments of the PSD that vary throughout
# the atmosphere.
#
# For the rain example considered here, we use a single moment corresponding
# to the mass density of the rain.
rain_first_moment = pyarts.arts.ScatteringSpeciesProperty(
    "rain", pyarts.arts.ParticulateProperty("MassDensity")
)
psd = MGDSingleMoment(rain_first_moment, "Wang16", 270, 300, False)
rain = ScatteringHabit(rain_habit, psd)

ws.scattering_species = [rain]

# %% Demonstration
#
# The scattering habit can now be used to calculate the bulk properties of our
# rain scattering species for any point in the atmosphere.

point = AtmPoint()
point["t"] = 280
point[rain_first_moment] = 1e-4

bulk_props = rain.get_bulk_scattering_properties_tro_spectral(point, f_grid, 1)

# %% Grids and Planet
p_grid = load("p_grid.xml")
t_field = load("t_field.xml")
z_field = load("z_field.xml")
vmr_field = load("vmr_field.xml")
pbf_field = load("particle_bulkprop_field.xml")
pbf_names = load("particle_bulkprop_names.xml")
ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = t_field[0, 0, 0]

lat_grid = np.array([0.0])
lon_grid = np.array([0.0])
z_grid = z_field[..., 0, 0]

pressure = GriddedField3(
    "p", p_grid[..., None, None], ["alt", "lon", "lat"], (z_grid, lon_grid, lat_grid)
)
temperature = GriddedField3(
    "t", t_field, ["alt", "lon", "lat"], (z_grid, lon_grid, lat_grid)
)
n2 = GriddedField3(
    "N2", vmr_field[0], ["alt", "lon", "lat"], (z_grid, lon_grid, lat_grid)
)
o2 = GriddedField3(
    "O2", vmr_field[1], ["alt", "lon", "lat"], (z_grid, lon_grid, lat_grid)
)
h2o = GriddedField3(
    "H2O", vmr_field[2], ["alt", "lon", "lat"], (z_grid, lon_grid, lat_grid)
)
rwc = GriddedField3(
    "RWC", 1.0 * pbf_field[0], ["alt", "lon", "lat"], (z_grid, lon_grid, lat_grid)
)

ws.atmospheric_field["p"] = pressure
ws.atmospheric_field["t"] = temperature
ws.atmospheric_field["N2"] = n2
ws.atmospheric_field["O2"] = o2
ws.atmospheric_field["H2O"] = h2o
ws.atmospheric_field[rain_first_moment] = rwc
ws.atmospheric_field.top_of_atmosphere = 12.0e3

ws.propagation_matrix_scattering_spectral_agenda

# %% Checks and settings
ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

ws.disort_settings_agendaSetup(scattering_setting="ScatteringSpecies")
ws.disort_quadrature_dimension = 40
ws.disort_fourier_mode_dimension = 1
ws.propagation_matrix_scattering_spectral_agendaSet()
ws.disort_legendre_polynomial_dimension = 40
ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")


def calculate_tbs_disort():
    ws.disort_settings_agendaSetup(scattering_setting="ScatteringSpecies")
    ws.disort_spectral_radiance_fieldProfile(
        longitude=lon,
        latitude=lat,
        disort_quadrature_dimension=NQuad,
        disort_legendre_polynomial_dimension=40,
        disort_fourier_mode_dimension=1,
        max_step=100,
    )
    disort_stokes = [
        [ws.disort_spectral_radiance_field.data[f_ind, 0, 0, 0], 0.0, 0.0, 0.0]
        for f_ind in range(3)
    ]
    ws.spectral_radiance = disort_stokes
    ws.spectral_radianceApplyForwardUnit(ray_path_point=ws.ray_path[0])
    return ws.spectral_radiance.value.copy()[:, 0]


tbs_cloudy = calculate_tbs_disort()

# %% Cloudy-sky brightness temperatures
#
# The cloudy-sky brightness temperature are obviously off. I checked the
# single-scattering albedo and extinction matrix between ARTS 2.6 and the new
# implementation and they agree. Currently, I am assuming there is an issue in
# DISORT.

# ARTS 2.6 results: ``271.694859567588, 272.601957925916, 251.643215266136``
tbs_cloudy

# ## Clear-sky brightness temperatures
#
# The clearsky brightness temperatures, however, agree well with the ARTS 2.6 results.

# ARTS 2.6 results: ``298.566120236439 283.35611518369 251.643322551348``
ws.atmospheric_field[rain_first_moment] = 0.0
tbs_clear = calculate_tbs_disort()
tbs_clear
