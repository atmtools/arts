import pyarts
import numpy as np
import xarray as xa
from dataclasses import dataclass
import matplotlib.pyplot as plt

NQuad = 16
max_level_step = 1e3
atm_latitude = 0.0
atm_longitude = 0.0
solar_latitude = 0.0
solar_longitude = 0.0
surface_temperature = 293.0
surface_reflectivity = 0.05
cutoff = ["ByLine", 750e9]
remove_lines_percentile = 70
sunfile = "star/Sun/solar_spectrum_QUIET.xml"
planet = "Earth"

xarr = pyarts.data.xarray_open_dataset("atm.nc")

ws = pyarts.Workspace()

ws.frequency_grid = pyarts.arts.convert.kaycm2freq(np.linspace(500, 2500, 1001))
ws.atmospheric_field = pyarts.data.to_atmospheric_field(xarr)

v = pyarts.data.to_absorption_species(ws.atmospheric_field)

ws.absorption_species = v
ws.ReadCatalogData(ignore_missing=True)
ws.propagation_matrix_agendaAuto(T_extrapolfac=1e9)

for band in ws.absorption_bands:
    ws.absorption_bands[band].cutoff = cutoff[0]
    ws.absorption_bands[band].cutoff_value = cutoff[1]

ws.absorption_bands.keep_hitran_s(remove_lines_percentile)

ws.surface_fieldSetPlanetEllipsoid(option=planet)

sun = pyarts.arts.GriddedField2.fromxml(sunfile)

ws.surface_field["t"] = surface_temperature

ws.sunFromGrid(
    sun_spectrum_raw=sun,
    latitude=solar_latitude,
    longitude=solar_longitude,
)

ws.disort_quadrature_dimension = NQuad
ws.disort_fourier_mode_dimension = 1
ws.disort_legendre_polynomial_dimension = 1

ws.ray_pathGeometricDownlooking(
    latitude=atm_latitude,
    longitude=atm_longitude,
    max_step=max_level_step,
)
ws.ray_path_atmospheric_pointFromPath()

ws.ray_path_frequency_gridFromPath()
ws.ray_path_propagation_matrixFromPath()
ws.disort_settingsInit()
ws.disort_settingsOpticalThicknessFromPath()
ws.disort_settingsLayerThermalEmissionLinearInTau()
ws.disort_settingsSurfaceEmissionByTemperature(ray_path_point=ws.ray_path[0])
ws.disort_settingsCosmicMicrowaveBackgroundRadiation()
ws.disort_settingsSurfaceLambertian(value=surface_reflectivity)
ws.disort_settingsNoSingleScatteringAlbedo()
ws.disort_settingsNoFractionalScattering()
ws.disort_settingsNoLegendre()
ws.disort_settingsSetSun(ray_path_point=ws.ray_path[-1])
ws.disort_spectral_flux_fieldCalc()

plt.semilogy(pyarts.arts.convert.freq2kaycm(ws.frequency_grid),
             ws.disort_spectral_flux_field[:, 1])

f, s = pyarts.plots.AtmField.plot(ws.atmospheric_field,
                           alts=np.linspace(0, ws.atmospheric_field.top_of_atmosphere))
