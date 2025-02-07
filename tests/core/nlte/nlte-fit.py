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

ws.ray_path_observersFluxProfile(n=11)
ws.ray_path_fieldFromObserverAgenda()

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

# import matplotlib.pyplot as plt
# for x in levels:
#     plt.plot(ws.atmospheric_field[x].data.flatten())
# plt.show()

ws.atmospheric_fieldFitNonLTE(
    collision_data=collision_data, levels=levels, convergence_criterion=1e-2
)

# for x in levels:
#     plt.plot(ws.atmospheric_field[x].data.flatten())
# plt.show()

ref = [
    0.15155419,
    0.16352512,
    0.18911686,
    0.21319294,
    0.24843466,
    0.29213797,
    0.33774595,
    0.37876021,
    0.40840548,
    0.42642776,
    0.43561727,
]

# FIXME: The values are not close enough, there's a difference between compilers
assert np.allclose(
    ws.atmospheric_field.nlte[levels[0]].data.flatten()[::10], ref, rtol=1e-1
), f"{ws.atmospheric_field.nlte[levels[0]].data.flatten()[::10]} vs {ref}"
