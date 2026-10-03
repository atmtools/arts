"""VDISORT, RT4 and RT3 run the same way on shared problems.

Runs every preset problem of pyarts3.polradtran through every solver that
can represent it (VDISORT always; RT4 for thermal problems with
nstokes <= 2; RT3 unless a beam falls on a non-Lambertian surface) and
requires every pair to agree to Comparison.tolerance(), 10 times the bound
on RT3's and RT4's first-order doubling error.  The problem set-up, the
solver mappings and the tolerance live in pyarts3.polradtran, so this test
checks exactly what pyarts3.polradtran.plot() shows.

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

# Which solvers each preset must reach, and that it exercises polarization
EXPECTED = [
    (["VDISORT", "RT4", "RT3"], (1e-3,)),  # thermal Rayleigh, nstokes 2: Q
    (["VDISORT", "RT4", "RT3"], ()),  # thermal Rayleigh, nstokes 1
    (["VDISORT", "RT4", "RT3"], (1e-4,)),  # thermal Mie: Q
    (["VDISORT", "RT3"], (0.1, 0.1)),  # solar Rayleigh: Q, U
    (["VDISORT", "RT3"], (0.05, 0.05, 1e-5)),  # solar + thermal multilayer: Q, U, V
]

cases = prt.cases.all()
assert len(cases) == len(EXPECTED)
for problem, (solvers, polarization) in zip(cases, EXPECTED):
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
        import matplotlib.pyplot as plt

        prt.plot(comparison, level=0, azimuth=90.0)
        plt.show()
