import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["H2O"])

ws.absorption_bands.readxml("nlte_lines.xml")

toa = 4.4825000e05
ws.atmospheric_fieldInit(toa=toa)

ws.atmospheric_field["t"] = pyarts.arts.GriddedField3.fromxml("t.xml")
ws.atmospheric_field["p"] = pyarts.arts.GriddedField3.fromxml("p.xml")
ws.atmospheric_field["N2"] = 0.0
ws.atmospheric_field["O2"] = 0.0
ws.atmospheric_field["H2O"] = 1.0
ws.atmospheric_field["CO2"] = 0.0
ws.atmospheric_field["H2"] = 0.0
ws.atmospheric_field["He"] = 0.0

ws.surface_fieldGanymede()
ws.surface_field["t"] = ws.atmospheric_field["t"].data[0, 0, 0]

line_f0 = 556936000000.0
ws.frequency_grid = np.linspace(-5e6, 5e6, 11) + line_f0

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.propagation_matrix_agendaAuto()

pos = [toa, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)

ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()
lte = ws.spectral_radiance[:, 0] * 1.0

ws.absorption_bandsSetNonLTE()
ws.atmospheric_fieldInitializeNonLTE()
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()
nlte = ws.spectral_radiance[:, 0] * 1.0

assert np.allclose(nlte, lte, rtol=1e-3)

