import numpy as np
import pyarts3 as pyarts

# ARTS2 reference case: 0.488 micrometres, 283.15 K, liquid-water density.
frequency = pyarts.arts.convert.wavelen2freq(0.488e-6)
reference_index = 1.33821937010893

n = pyarts.arts.physics.refractive_index_water_visible_nir_harvey98(
    frequency=frequency,
    temperature=283.15,
    density=1000.0,
)
assert np.isclose(n, reference_index, rtol=0.0, atol=1e-13)


# The Python interface must broadcast all physical inputs.
wavelengths = np.array([0.488, 0.589, 1.0]) * 1e-6
temperatures = np.array([[283.15], [300.0]])
densities = np.array([[1000.0], [997.0]])
indices = pyarts.arts.physics.refractive_index_water_visible_nir_harvey98(
    frequency=pyarts.arts.convert.wavelen2freq(wavelengths),
    temperature=temperatures,
    density=densities,
)
assert indices.shape == (2, 3)
assert np.all(np.isfinite(indices))

try:
    pyarts.arts.physics.refractive_index_water_visible_nir_harvey98(
        frequency=pyarts.arts.convert.wavelen2freq(2.0e-6),
        temperature=283.15,
        density=1000.0,
    )
except RuntimeError:
    pass
else:
    raise AssertionError("Harvey98 accepted a wavelength outside its validity range")


# The single-propagation-matrix agenda supplies n - 1 to refractive path
# stepping.  Test both agendas together, independently of any absorption model.
ws = pyarts.Workspace()
ws.jac_targetsInit()
ws.freq = frequency
ws.freq_wind_shift_jac = [0.0, 0.0, 0.0]
ws.select_species = "AIR"
ws.atm_point = pyarts.arts.AtmPoint()
ws.atm_point.temperature = 283.15
ws.atm_point.pressure = 1e5
ws.atm_point["H2O"] = 0.01
ws.ray_point = pyarts.arts.PropagationPathPoint()

# Default agenda setup derives water density from the atmospheric H2O VMR.
ws.single_propmat_agendaSetWaterVisibleNIRHarvey98()
ws.single_propmat_agendaExecute()
assert 0.0 < ws.single_dispersion < 1e-3

ws.single_propmat_agendaSetWaterVisibleNIRHarvey98(
    water_mass_density=1000.0
)
ws.single_propmat_agendaExecute()
assert np.isclose(ws.single_dispersion, reference_index - 1.0, rtol=0.0, atol=1e-13)

ws.atm_fieldInit(toa=100e3)
ws.surf_fieldEarth()
ws.max_stepsize = 1000.0

start = pyarts.arts.PropagationPathPoint()
start.pos_type = "atm"
start.los_type = "atm"
start.pos = [50e3, 0.0, 0.0]
start.los = [100.0, 0.0]
ws.ray_path = [start]

ws.ray_point_back_propagation_agendaSet(option="RefractiveStepwise")
ws.ray_point_back_propagation_agendaExecute()
refractive_pos = np.array(ws.ray_point.pos)
refractive_los = np.array(ws.ray_point.los)
refractive_nreal = ws.ray_point.nreal

assert np.isclose(refractive_nreal, reference_index, rtol=0.0, atol=1e-13)

ws.ray_path = [start]
ws.ray_point_back_propagation_agendaSet(option="GeometricStepwise")
ws.ray_point_back_propagation_agendaExecute()
geometric_point = ws.ray_point

assert not np.allclose(refractive_pos, geometric_point.pos)
assert not np.isclose(refractive_los[0], geometric_point.los[0])
