import pyarts3 as pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(118.5e9 - 30e9, 119e9 + 30e9, 14001)

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=120e9)
del ws.abs_bands["O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 1 1"]

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Settings

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric(add_crossings=1, remove_non_crossings=1)

nlimb = 15
nup = 30
ndown = 120

# %% Checks

# Configuration for the ray path and flux calculation
ws.ray_path_observersFieldProfilePseudo2D(nup=nup, ndown=ndown, nlimb=nlimb)
# Generate the ray path field based on observer's agenda
ws.ray_path_fieldFromObserverAgenda()
# Calculate spectral flux profile using the ray path and atmospheric data
ws.spectral_flux_profileFromPathField(
    alt_grid=ws.atm_field["t"].data.grids[0]
)

# Atmospheric point definition along the ray path
ws.atm_profile = [
    ws.atm_field(x, 0, 0) for x in ws.atm_field["t"].data.grids[0]
]

# Integrate line flux profile for NLTE calculation
ws.nlte_line_flux_profileIntegrate()

# We are only close

assert np.allclose(
    ws.nlte_line_flux_profile["O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"],
    [
        3.81600630971944e-15,
        3.767359969177587e-15,
        3.7431764062592374e-15,
        3.720047296554584e-15,
        3.693869872354687e-15,
        3.662728223602894e-15,
        3.628250526592601e-15,
        3.591223421186791e-15,
        3.5520819474662738e-15,
        3.5110755656352577e-15,
        3.4685814048628366e-15,
        3.4245359385658345e-15,
        3.3793899037737243e-15,
        3.3336123751690905e-15,
        3.2871024717625563e-15,
        3.2395424614287704e-15,
        3.19086716538334e-15,
        3.1486698604832502e-15,
        3.1248240147012255e-15,
        3.118128095994502e-15,
        3.118900969333555e-15,
        3.1240681398056813e-15,
        3.1319375633920453e-15,
        3.1396352860671e-15,
        3.1457404539815727e-15,
        3.151385891297707e-15,
        3.1547355508761993e-15,
        3.163724787455442e-15,
        3.1757556971694315e-15,
        3.189746791011155e-15,
        3.205024184701511e-15,
        3.222190177679822e-15,
        3.250501708915215e-15,
        3.33028097112648e-15,
        3.534658071012364e-15,
        3.90431042796717e-15,
        5.163801133874627e-15,
        6.575174943960289e-15,
        6.685163771122741e-15,
        4.993503540873373e-15,
        2.8167350310086086e-15,
        1.3346612442515194e-15,
        5.606394939077885e-16,
        2.27424025731406e-16,
        1.0187432722431091e-16,
        5.822386490600417e-17,
        6.433938582408049e-17,
        1.1661128165526783e-16,
        3.2692276684154534e-16,
        8.170438260755405e-16,
    ],
)
