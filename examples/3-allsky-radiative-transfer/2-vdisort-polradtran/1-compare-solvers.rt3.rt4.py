"""
VDISORT compared with RT4 and RT3

This script runs the same plane-parallel problems through VDISORT and
Evans' doubling-adding solvers RT4 (thermal, azimuthally symmetric, [I, Q])
and RT3 (solar beam and all Fourier azimuth modes, [I, Q, U, V]) and plots
the radiances and their differences.  All solvers get the same streams and
phase matrix, so they differ only by RT3's and RT4's doubling error, which is
below max_delta_tau / mu0 of the largest radiance.

It needs an ARTS built with -DENABLE_RT4=ON and -DENABLE_RT3=ON.
"""

import os

import matplotlib.pyplot as plt
import pyarts3 as pyarts

prt = pyarts.polradtran

# %% The preset problems, which CI also runs

for problem in prt.cases.all():
    comparison = prt.compare(problem)
    print(comparison)
    comparison.check()

# %% Thermal emission over a Fresnel surface: all three solvers, at the surface

comparison = prt.compare(prt.cases.thermal_rayleigh(nstokes=2))
fig, ax = prt.plot(comparison, level=comparison.problem.nlay)

# %% A solar beam: U is largest away from the principal plane

comparison = prt.compare(prt.cases.solar_rayleigh())
fig, ax = prt.plot(comparison, level=0, azimuth=90.0)

# %% A problem of your own: a polarizing Henyey-Greenstein cloud (g = 0.6)
# under a Rayleigh layer.  legendre_series projects any scattering matrix
# onto RT3's Legendre columns; this one is regular at 0 and 180 degrees
# (F22 = F33 forward, F22 = -F33 backward)


def polarized_hg(c, g=0.6, p=0.4):
    hg = (1 - g**2) / (1 + g**2 - 2 * g * c) ** 1.5
    rayleigh_shape = 1 + c**2
    return (
        hg,
        -p * hg * (1 - c**2) / rayleigh_shape,
        hg * 2 * c / rayleigh_shape,
        0 * c,
        hg,
        hg * 2 * c / rayleigh_shape,
    )


legendre = prt.legendre_series(
    polarized_hg, degree=13
)  # RT3 keeps up to 2 nmu - 3 = 13

problem = prt.Problem(
    name="solar polarizing HG cloud (g = 0.6) under a Rayleigh layer",
    nstokes=4,
    nmu=8,
    nfourier=8,
    frequency=pyarts.arts.convert.wavelen2freq(0.6e-6),
    dz=[1.0, 2.0],
    gas=[0.0, 0.01],
    sets=[
        prt.ScatteringSet(0.2, 0.2, prt.RAYLEIGH),
        prt.ScatteringSet(1.0, 0.95, legendre),
    ],
    set_index=[0, 1],
    temperature=[1.0, 1.0, 1.0],
    sky=0.0,
    surface_temperature=0.0,
    surface=("lambertian", 0.2),
    thermal=False,
    beam=(1.0, 0.5),
)
comparison = prt.compare(problem)
print(comparison)
comparison.check()
fig, ax = prt.plot(comparison, level=0, azimuth=60.0)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()
