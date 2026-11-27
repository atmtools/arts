"""Gas-scattering species and generic first-order solar scattering."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyarts3 as pyarts


A = pyarts.arts
frequency = 6e14
cross_section = 1e-32

# The coefficient model is per molecule and therefore follows number density.
gas = A.GasScatterer.constant_isotropic(cross_section)
point = A.AtmPoint()
point.pressure = 101325.0
point.temperature = 288.15
angles = A.IrregularZenithAngleGrid([0.0, 90.0, 180.0])
bulk = gas.tro_gridded(point, A.AtmSpectralPoint(), [frequency], angles)

point_twice_pressure = A.AtmPoint()
point_twice_pressure.pressure = 2.0 * point.pressure
point_twice_pressure.temperature = point.temperature
bulk_twice_pressure = gas.tro_gridded(
    point_twice_pressure, A.AtmSpectralPoint(), [frequency], angles
)
assert np.allclose(
    np.asarray(bulk_twice_pressure.extinction_matrix),
    2.0 * np.asarray(bulk.extinction_matrix),
)
assert np.all(np.asarray(bulk.absorption_vector) == 0.0)

with TemporaryDirectory() as directory:
    filename = Path(directory) / "gas_species.xml"
    pyarts.xml.save(A.ArrayOfScatteringSpecies([gas]), str(filename))
    restored = pyarts.xml.load(str(filename))
    restored_bulk = restored.tro_gridded(
        point, A.AtmSpectralPoint(), [frequency], angles
    )
    assert np.allclose(
        np.asarray(restored_bulk.extinction_matrix),
        np.asarray(bulk.extinction_matrix),
    )

# Exercise the same species through the generic first-order solar path.  For
# exact backscattering and zero depolarization, Rayleigh F11 is 3/2 while the
# isotropic phase function is one.  Extinction and transmission are identical.
ws = pyarts.Workspace()
ws.freq_grid = [frequency]
ws.atm_field.top_of_atmosphere = 80e3
ws.atm_field["t"] = point.temperature
ws.atm_field["p"] = point.pressure
ws.abs_speciesSet(species=[])
ws.spectral_propmat_agendaAuto()

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("t")] = point.temperature
ws.sunBlackbody()
ws.suns = [ws.sun]
ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric()
ws.spectral_propmat_scat_spectral_agendaSet(option="FromSpeciesTRO")
ws.legendre_degree = 2

ws.ray_pathGeometric(pos=[90e3, 0.0, 0.0], los=[180.0, 0.0], max_stepsize=2e3)
ws.ray_path_suns_pathFromPathObserver(just_hit=1)


def calculate(species):
    ws.atm_scat_species = [species]
    ws.spectral_radClearskyScattering()
    return np.asarray(ws.spectral_rad).copy()[0]


isotropic = calculate(A.GasScatterer.constant_isotropic(cross_section))
rayleigh = calculate(
    A.GasScatterer(
        A.ConstantGasScattering(cross_section), A.RayleighGasScattering(0.0)
    )
)

assert isotropic[0] > 0.0
assert np.allclose(rayleigh[0] / isotropic[0], 1.5, rtol=1e-10)
assert np.allclose(rayleigh[1:], 0.0, atol=1e-40)

# A horizontal view is a 90-degree scattering geometry for the overhead sun.
# Rayleigh scattering is then fully linearly polarized, while the tiny U term
# only reflects changing local frames along the curved limb path.
ws.ray_pathGeometric(pos=[1.0, 0.0, 0.0], los=[90.0, 0.0], max_stepsize=1e4)
ws.ray_path_suns_pathFromPathObserver(just_hit=1)
right_angle = calculate(
    A.GasScatterer(
        A.ConstantGasScattering(1e-34), A.RayleighGasScattering(0.0)
    )
)
assert np.allclose(right_angle[1] / right_angle[0], -1.0, rtol=1e-6)
assert np.allclose(right_angle[2:], 0.0, atol=1e-30)
