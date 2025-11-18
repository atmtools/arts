import pyarts3 as pyarts

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
ws.atm_fieldInitializeNonLTE(normalization=0.75)
ws.abs_bandsSetNonLTE()

ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")
ws.spectral_propmat_agendaAuto()

collision_data = pyarts.arts.QuantumIdentifierGriddedField1Map.fromxml("Cij.xml")

ws.freq_gridFitNonLTE(nf=401, df=1e-4)

levels = pyarts.arts.ArrayOfQuantumIdentifier(
    [
        "H2O-161 J 1 1 Ka 0 0 Kc 1 1",
        "H2O-161 J 1 1 Ka 1 1 Kc 0 0",
        "H2O-161 J 2 2 Ka 1 1 Kc 2 2",
        "H2O-161 J 2 2 Ka 2 2 Kc 1 1",
        "H2O-161 J 3 3 Ka 0 0 Kc 3 3",
        "H2O-161 J 3 3 Ka 1 1 Kc 2 2",
        "H2O-161 J 3 3 Ka 2 2 Kc 1 1",
    ]
)

ws.ray_pathInit(pos = [ws.atm_field.top_of_atmosphere, 0, 0], los=[180, 0])
ws.ray_pathSetGeometricExtremes()
ws.ray_pathAddGeometricGridCrossings()
ws.ray_path.pop(0)

ws.disort_quadrature_dimension = 40
ws.disort_fourier_mode_dimension = 1
ws.disort_legendre_polynomial_dimension = 1
ws.disort_settings_agendaSetup()
ws.disort_settings_agendaExecute()
