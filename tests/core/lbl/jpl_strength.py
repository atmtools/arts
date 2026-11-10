"""The line strengths read from the JPL catalog must reproduce the JPL intensities

NH3-4111 has a TIPS partition function, which includes a nuclear spin factor that
the JPL partition function leaves out, so its JPL intensity is only reproduced at
the reference temperature.  HNC-142 has the JPL partition function table as its
partition function, so its JPL intensity is reproduced at all the temperatures of
that table.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

c = pyarts.arts.constants.c
h = pyarts.arts.constants.h
k = pyarts.arts.constants.k

# JPL line records: NH3 1(0)-0(0) and HNC 1-0
NH3 = "  572498.1599  0.0020 -1.5128 3    0.0000 12 -170021303 1 0 0       0 0 1       "
HNC = "   90663.5930  0.0090 -2.4372 2    0.0000  3 -27002 101 1           0           "

# JPL partition functions of the entries (17002 and 27002)
QJPL = {
    "NH3-4111": {300.0: 588.7816},
    "HNC-142": {300.0: 138.311, 225.0: 103.801, 150.0: 69.301, 75.0: 34.814},
}


def jpl_intensity(record, T, Q):
    """The JPL intensity [m2 Hz] of the record at T for the JPL partition functions Q"""
    f = float(record[0:13]) * 1e6
    s300 = 10 ** float(record[21:29]) * 1e-12
    E = float(record[31:41]) * 100 * h * c

    return (
        s300
        * Q[300.0]
        / Q[T]
        * np.exp(-E / k * (1 / T - 1 / 300))
        * np.expm1(-h * f / (k * T))
        / np.expm1(-h * f / (k * 300))
    )


def arts_intensity(line, isot, T):
    """The line strength of ARTS [m2 Hz] in the convention of the catalogs"""
    Q = pyarts.arts.SpeciesIsotope(isot).Q(T)
    f = line.f0
    return line.s(T, Q) * f * -np.expm1(-h * f / (k * T)) * c**2 / (8 * np.pi)


ws = pyarts.Workspace()

for record, isot in [(NH3, "NH3-4111"), (HNC, "HNC-142")]:
    with tempfile.TemporaryDirectory() as tmp:
        fn = os.path.join(tmp, "lines.cat")
        with open(fn, "w") as fil:
            fil.write(record + "\n")
        ws.abs_bandsReadJPL(file=fn)

    assert len(ws.abs_bands) == 1
    key, band = list(ws.abs_bands.items())[0]
    assert str(key.isot) == isot, f"{key.isot} != {isot}"
    line = band.lines[0]

    for T in QJPL[isot]:
        x = arts_intensity(line, isot, T)
        y = jpl_intensity(record, T, QJPL[isot])
        assert np.isclose(x, y, rtol=1e-12, atol=0), (
            f"{isot} at {T} K: ARTS {x} != JPL {y}"
        )
