import pyarts
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
ws.atmospheric_fieldInitializeNonLTE(normalization=0.75)
ws.absorption_bandsSetNonLTE()

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric(add_crossings=True, remove_non_crossings=True)
ws.propagation_matrix_agendaAuto()

ws.ray_path_fieldFluxProfile(dza=10)

collision_data = pyarts.arts.QuantumIdentifierGriddedField1Map.fromxml("Cij.xml")

ws.frequency_gridFitNonLTE(nf=401, df=1e-4)

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

ws.atmospheric_fieldFitNonLTE(
    collision_data=collision_data, levels=levels, convergence_criterion=5e-1
)

ref = [0.14906513, 0.16333853, 0.17919964, 0.19362265, 0.21227294,
       0.23620497, 0.26243519, 0.28530674, 0.30159183, 0.3133241 ,
       0.32257162]

assert np.allclose(
    ws.atmospheric_field.nlte[levels[0]].data.flatten()[::10], ref
), f"{ws.atmospheric_field.nlte[levels[0]].data.flatten()[::10]} vs {ref}"
