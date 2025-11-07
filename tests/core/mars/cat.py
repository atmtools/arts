import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["H2O-161"])

ws.surface_fieldMars()
ws.surface_field["t"] = pyarts.arts.GriddedField2.fromxml(
    "planets/Mars/Ls0.day.dust-medium/surface_temperature.xml"
).make_geodetic()
ws.surface_field["h"] = pyarts.arts.GriddedField2.fromxml(
    "planets/Mars//surface_elevation.xml"
).make_geodetic()
ws.surface_field["t"].set_extrapolation("Nearest")
ws.surface_field["h"].set_extrapolation("Nearest")

pos = [200e3, 0, 0]
los = [180, 0]
f0 = 556.935985e9
ws.frequency_grid = f0 + np.linspace(-5e9, 5e9, 101)

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

out = []
for abs_scenario in ["FullMars", "IsotEarth", "FullEarth"]:
    print(abs_scenario)
    if abs_scenario == "FullMars":
        ws.abs_bandsReadSpeciesSplitARTSCAT(basename="spectroscopy/Perrin/")
        ws.abs_bands = ws.abs_bands.extract_species("H2O-161")
        ws.propagation_matrix_agendaAuto()
        ws.atmospheric_fieldRead(
            toa=150e3,
            basename="planets/Mars/Ls0.day.dust-medium/Ls0.day.dust-medium.sol-avg/",
            missing_is_zero=True,
            extrapolation="Nearest",
        )
        ws.atmospheric_fieldAppendLineIsotopologueData(
            basename="planets/Mars/isotopologue_ratios/", replace_existing=True
        )
        assert 6 == len(
            ws.atmospheric_field.specs
        ), "Reading ARTSCAT-4 but not getting 6 species"
    elif abs_scenario == "IsotEarth":
        ws.abs_bandsReadSpeciesSplitARTSCAT(basename="spectroscopy/Perrin/")
        ws.abs_bands = ws.abs_bands.extract_species("H2O-161")
        ws.propagation_matrix_agendaAuto()
        ws.atmospheric_fieldRead(
            toa=150e3,
            basename="planets/Mars/Ls0.day.dust-medium/Ls0.day.dust-medium.sol-avg/",
            missing_is_zero=True,
            extrapolation="Nearest",
        )
    elif abs_scenario == "FullEarth":
        ws.ReadCatalogData()
        ws.abs_bandsSelectFrequencyByLine(fmin=550e9, fmax=560e9)
        ws.propagation_matrix_agendaAuto()
        ws.atmospheric_fieldRead(
            toa=150e3,
            basename="planets/Mars/Ls0.day.dust-medium/Ls0.day.dust-medium.sol-avg/",
            missing_is_zero=True,
            extrapolation="Nearest",
        )
    ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
    ws.spectral_radianceClearskyEmission()
    ws.spectral_radianceApplyUnitFromSpectralRadiance()

    out.append(1.0*ws.spectral_radiance[:, 0])

assert np.all(out[0] != out[1]), "FullMars and IsotEarth should be different, changed isotopologue ratios"
assert np.all(out[0] != out[2]), "FullMars and FullEarth should be different, changed absorption bands"
assert np.all(out[1] != out[2]), "IsotEarth and FullEarth should be different, changed absorption bands"
