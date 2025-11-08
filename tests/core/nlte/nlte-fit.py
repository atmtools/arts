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

ws.surface_fieldGanymede()
ws.surface_field["t"] = ws.atm_field["t"].data[0, 0, 0]
ws.atm_fieldInitializeNonLTE(normalization=0.75)
ws.abs_bandsSetNonLTE()

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric(add_crossings=True, remove_non_crossings=True)
ws.propagation_matrix_agendaAuto()

collision_data = pyarts.arts.QuantumIdentifierGriddedField1Map.fromxml("Cij.xml")

ws.frequency_gridFitNonLTE(nf=301, df=1e-4)

levels = pyarts.arts.ArrayOfQuantumLevelIdentifier(
    [
        "H2O-161 J 1 Ka 0 Kc 1",
        "H2O-161 J 1 Ka 1 Kc 0",
        "H2O-161 J 2 Ka 1 Kc 2",
        "H2O-161 J 2 Ka 2 Kc 1",
        "H2O-161 J 3 Ka 0 Kc 3",
        "H2O-161 J 3 Ka 1 Kc 2",
        "H2O-161 J 3 Ka 2 Kc 1",
    ]
)

ws.atm_profileFromGrid()
ws.atm_profileFitNonLTE(
    collision_data=collision_data,
    levels=levels,
    dza=15,
    consider_limb=0,
    convergence_limit=1e-2,
)
ws.atm_fieldFromProfile()

assert np.allclose(
    ws.atm_field.nlte[levels[0]].data.flatten(),
    np.array(
        [
            0.1471907,
            0.14725853,
            0.14820754,
            0.14922162,
            0.15026332,
            0.15132991,
            0.1524223,
            0.1535408,
            0.1546861,
            0.15585901,
            0.15705965,
            0.15828437,
            0.15947635,
            0.16076186,
            0.16208209,
            0.16342939,
            0.1647907,
            0.16611868,
            0.16701486,
            0.16768155,
            0.16839643,
            0.16906924,
            0.16981641,
            0.17053676,
            0.17133858,
            0.17212494,
            0.17299928,
            0.17387131,
            0.17483953,
            0.17582054,
            0.17690861,
            0.17802791,
            0.17926956,
            0.18054853,
            0.18187759,
            0.18330887,
            0.18485518,
            0.18652491,
            0.18832625,
            0.19026404,
            0.19234227,
            0.19456418,
            0.19692697,
            0.19941492,
            0.20204203,
            0.20482096,
            0.20774923,
            0.21082012,
            0.21403305,
            0.21738219,
            0.22085925,
            0.22445576,
            0.22815107,
            0.23192582,
            0.23579137,
            0.23975452,
            0.24380735,
            0.24794062,
            0.25214388,
            0.2564086,
            0.26072478,
            0.2650811,
            0.26946922,
            0.27388121,
            0.27830467,
            0.28273146,
            0.28715527,
            0.29158426,
            0.29598041,
            0.30033719,
            0.30464426,
            0.30889241,
            0.31307362,
            0.31718234,
            0.32120745,
            0.32513935,
            0.32892355,
            0.33264331,
            0.33630525,
            0.33990629,
            0.34344474,
            0.34691791,
            0.3503257,
            0.35366721,
            0.35694104,
            0.36014657,
            0.36328462,
            0.36635417,
            0.36935689,
            0.37229075,
            0.37516023,
            0.37797271,
            0.38073064,
            0.38344661,
            0.38613017,
            0.38880416,
            0.39150159,
            0.39427452,
            0.39720971,
            0.40044572,
            0.40420514,
        ]
    ),
)
