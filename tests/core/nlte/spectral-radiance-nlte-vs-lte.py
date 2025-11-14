import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["H2O"])

ws.abs_bands.readxml("nlte_lines.xml")

toa = 4.4825000e05
ws.atm_fieldInit(toa=toa)

ws.atm_field["t"] = pyarts.arts.GriddedField3.fromxml("t.xml")
ws.atm_field["p"] = pyarts.arts.GriddedField3.fromxml("p.xml")
ws.atm_field["N2"] = 0.0
ws.atm_field["O2"] = 0.0
ws.atm_field["H2O"] = 1.0
ws.atm_field["CO2"] = 0.0
ws.atm_field["H2"] = 0.0
ws.atm_field["He"] = 0.0

ws.surf_fieldGanymede()
ws.surf_field["t"] = ws.atm_field["t"].data[0, 0, 0]

line_f0 = 556936000000.0
ws.freq_grid = np.linspace(-5e6, 5e6, 11) + line_f0

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.spectral_propmat_agendaAuto()

pos = [toa, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)

ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()
lte = ws.spectral_radiance[:, 0] * 1.0

ws.abs_bandsSetNonLTE()
ws.atm_fieldInitializeNonLTE()
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()
nlte = ws.spectral_radiance[:, 0] * 1.0

assert np.allclose(nlte, lte, rtol=1e-3)

