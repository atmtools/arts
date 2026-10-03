"""VDISORT, RT4 and RT3 run the same way on shared problems.

Compares the polarized discrete-ordinate solver VDISORT with Evans'
doubling-adding solvers RT4 (thermal, azimuthally symmetric, [I, Q]) and RT3
(solar beam and all Fourier azimuth modes, [I, Q, U, V]) of polradtran,
which ARTS builds when a Fortran compiler is found.

Every preset problem of pyarts3.polradtran, and a problem defined here with
pyarts3.polradtran.legendre_series, runs through every solver that can
represent it (VDISORT always; RT4 for thermal problems with nstokes <= 2; RT3
unless a beam falls on a non-Lambertian surface).  All solvers get the same
streams and phase matrix, so they differ only by RT3's and RT4's
first-order doubling error, and every pair must agree to
Comparison.tolerance(), 10 times the bound on that error.  The problem
set-up, the solver mappings and the tolerance live in pyarts3.polradtran, so
this test checks exactly what pyarts3.polradtran.plot() shows; run it without
ARTS_HEADLESS to see the plots.

Collected only when ENABLE_RT3 and ENABLE_RT4 are both on (the ".rt3." and
".rt4." in the file name).
"""

import os

import pyarts3 as pyarts

prt = pyarts.polradtran
assert prt.available_solvers() == [
    "VDISORT",
    "RT4",
    "RT3",
], "collected only with ENABLE_RT3=ON and ENABLE_RT4=ON"


def polarized_hg(c, g=0.6, p=0.4):
    """Henyey-Greenstein with Rayleigh-like polarization.

    Regular at 0 and 180 degrees: F22 = F33 forward and F22 = -F33 backward.
    """
    hg = (1 - g**2) / (1 + g**2 - 2 * g * c) ** 1.5
    rayleigh_shape = 1 + c**2
    f33 = hg * 2 * c / rayleigh_shape
    return hg, -p * hg * (1 - c**2) / rayleigh_shape, f33, 0 * c, hg, f33


hg_cloud = prt.Problem(
    name="solar polarizing HG cloud (g = 0.6) under a Rayleigh layer",
    nstokes=4,
    nmu=8,
    nfourier=8,
    frequency=pyarts.arts.convert.wavelen2freq(0.6e-6),
    dz=[1.0, 2.0],
    gas=[0.0, 0.01],
    sets=[
        prt.ScatteringSet(0.2, 0.2, prt.RAYLEIGH),
        # RT3 keeps Legendre degrees up to 2 nmu - 3 = 13 with double-Gauss streams
        prt.ScatteringSet(1.0, 0.95, prt.legendre_series(polarized_hg, degree=13)),
    ],
    set_index=[0, 1],
    temperature=[1.0, 1.0, 1.0],
    sky=0.0,
    surface_temperature=0.0,
    surface=("lambertian", 0.2),
    thermal=False,
    beam=(1.0, 0.5),
)

# Each problem, the solvers it must reach, the minimum max |Q|, |U|, |V| /
# max |I| that shows it exercises polarization, and where to plot it
CASES = [
    (prt.cases.thermal_rayleigh(2), ["VDISORT", "RT4", "RT3"], (1e-3,), (4, 0.0)),
    (prt.cases.thermal_rayleigh(1), ["VDISORT", "RT4", "RT3"], (), (4, 0.0)),
    (prt.cases.thermal_mie(), ["VDISORT", "RT4", "RT3"], (1e-4,), (2, 0.0)),
    (prt.cases.solar_rayleigh(), ["VDISORT", "RT3"], (0.1, 0.1), (0, 90.0)),
    (prt.cases.solar_thermal_multilayer(), ["VDISORT", "RT3"], (0.05, 0.05, 1e-5), (0, 90.0)),
    (hg_cloud, ["VDISORT", "RT3"], (0.1, 0.1), (0, 60.0)),
]
assert len(prt.cases.all()) == 5, "every preset must be listed in CASES"

for problem, solvers, polarization, (level, azimuth) in CASES:
    comparison = prt.compare(problem)
    print(comparison)
    assert (
        list(comparison.results) == solvers
    ), f"{problem.name}: ran {list(comparison.results)}, expected {solvers}"
    pol = comparison.polarization()
    for s, minimum in enumerate(polarization):
        assert (
            pol[s] > minimum
        ), f"{problem.name}: Stokes {s + 1} must be exercised (> {minimum} of I), got {pol[s]:.2e}"
    comparison.check()

    if "ARTS_HEADLESS" not in os.environ:
        prt.plot(comparison, level=level, azimuth=azimuth)

if "ARTS_HEADLESS" not in os.environ:
    import matplotlib.pyplot as plt

    plt.show()
