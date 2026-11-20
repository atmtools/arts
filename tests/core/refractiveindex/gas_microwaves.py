import numpy as np
import pyarts3 as pyarts


def workspace(atm_point):
    ws = pyarts.Workspace()
    ws.jac_targetsInit()
    ws.freq = 10e9
    ws.freq_wind_shift_jac = [0.0, 0.0, 0.0]
    ws.select_species = "AIR"
    ws.atm_point = atm_point
    ws.ray_point = pyarts.arts.PropagationPathPoint()
    return ws


# Earth model, including the optional-water dry-air case.
earth = pyarts.arts.AtmPoint()
earth.pressure = 1e5
earth.temperature = 300.0
earth["H2O"] = 0.01

ws = workspace(earth)
ws.single_propmat_agendaSetGasMicrowavesEarth()
ws.single_propmat_agendaExecute()

# Reference values of the Earth model, confirmed with independent code.
assert np.isclose(ws.single_dispersion, 2.999711111111111e-4, rtol=0.0, atol=1e-15)

dry = pyarts.arts.AtmPoint()
dry.pressure = earth.pressure
dry.temperature = earth.temperature
ws.atm_point = dry
ws.single_propmat_agendaExecute()
assert np.isclose(ws.single_dispersion, 2.5866666666666665e-4, rtol=0.0, atol=1e-15)


# General model: individual species are optional.  Supplying only CO2 at the
# reference pressure and temperature must give the pure-CO2 reference value.
planet = pyarts.arts.AtmPoint()
planet.pressure = 101325.0
planet.temperature = 273.15
planet["CO2"] = 0.4

ws.atm_point = planet
ws.single_propmat_agendaSetGasMicrowavesGeneral()
ws.single_propmat_agendaExecute()
assert np.isclose(ws.single_dispersion, 495.16e-6, rtol=0.0, atol=1e-15)

# Missing O2, H2, He, and H2O remain valid.  The supplied supported VMRs are
# normalized, so these deliberately need not sum to one.
planet["N2"] = 0.1
ws.atm_point = planet
ws.single_propmat_agendaExecute()
assert np.isclose(ws.single_dispersion, 454.89e-6, rtol=0.0, atol=1e-15)

# With no supported species the general model contributes zero.
empty = pyarts.arts.AtmPoint()
empty.pressure = planet.pressure
empty.temperature = planet.temperature
ws.atm_point = empty
ws.single_propmat_agendaExecute()
assert ws.single_dispersion == 0.0


# Confirm that the Earth model feeds the existing refractive stepper.
ws.atm_point = earth
ws.single_propmat_agendaSetGasMicrowavesEarth()
ws.atm_fieldInit(toa=100e3)
ws.surf_fieldEarth()
ws.max_stepsize = 1000.0

start = pyarts.arts.PropagationPathPoint()
start.pos_type = "atm"
start.los_type = "atm"
start.pos = [50e3, 0.0, 0.0]
start.los = [100.0, 0.0]
ws.ray_path = [start]

ws.single_propmat_agendaExecute()
expected_index = 1.0 + ws.single_dispersion
ws.ray_point_back_propagation_agendaSet(option="RefractiveStepwise")
ws.ray_point_back_propagation_agendaExecute()
assert np.isclose(ws.ray_point.nreal, expected_index, rtol=0.0, atol=1e-15)
