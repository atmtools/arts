import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

if pyarts.arts.globals.data.is_lgpl:
    print(
        "CKDMT models are not available in LGPL mode, compile with -DENABLE_ARTS_LGPL=0"
    )
    exit(0)

# Download catalogs
pyarts.data.download()

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

# Input file names are prefixed with the example number
prefix = "3-"

xarr = pyarts.data.xarray_open_dataset(prefix + "atmosphere.nc")

ws = pyarts.Workspace()

ws.freq_grid = pyarts.arts.convert.kaycm2freq(np.linspace(500, 2500, 1001))
ws.atm_field = pyarts.data.to_atm_field(xarr)

v = pyarts.data.to_abs_species(ws.atm_field)

ws.abs_species = v
ws.ReadCatalogData(ignore_missing=True)
ws.propagation_matrix_agendaAuto(T_extrapolfac=1e9)

for band in ws.abs_bands:
    ws.abs_bands[band].cutoff = cutoff[0]
    ws.abs_bands[band].cutoff_value = cutoff[1]

ws.abs_bands.keep_hitran_s(remove_lines_percentile)

ws.surface_fieldPlanet(option=planet)

sun = pyarts.arts.GriddedField2.fromxml(sunfile)

ws.surface_field["t"] = surface_temperature

ws.sunFromGrid(
    sun_spectrum_raw=sun,
    lat=solar_latitude,
    lon=solar_longitude,
)

ws.disort_quadrature_dimension = NQuad
ws.disort_fourier_mode_dimension = 1
ws.disort_legendre_polynomial_dimension = 1

ws.ray_pathGeometricDownlooking(
    lat=atm_latitude,
    lon=atm_longitude,
    max_stepsize=max_level_step,
)

ws.disort_settings_agendaSetup(
    sun_setting="Sun",
    surface_setting="ThermalLambertian",
    surface_lambertian_value=surface_reflectivity * np.ones_like(ws.freq_grid),
)

ws.disort_spectral_flux_fieldFromAgenda()

f, s = pyarts.plot(
    ws.disort_spectral_flux_field,
    freqs=pyarts.arts.convert.freq2kaycm(ws.freq_grid),
    alts=ws.disort_spectral_flux_field.alt_grid/1e3,
    select='down_diffuse',
    plotstyle='plot',
)
s.set_xlabel("Kaysers [cm$^{-1}$]")
s.set_ylabel("Altitude [km]")
s.set_title("Downward diffuse spectral flux from DISORT")

ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

alts = np.linspace(0, ws.atm_field.top_of_atmosphere)
f, s = pyarts.plot(
    ws.atm_field,
    alts=alts,
    ygrid=alts/1e3,
    apply_natural_scale=True,
)
for a in s.flatten():
    a.set_ylabel("Altitude [km]")
    lines = a.get_lines()
    if len(lines) > 0:
        a.set_xlabel(lines[0].get_label())
f.suptitle("Atmospheric field")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()
