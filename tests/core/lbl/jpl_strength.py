"""The line strengths read from the JPL catalog must reproduce the JPL intensities

NH3-4111 has a TIPS partition function, which includes a nuclear spin factor that
the JPL partition function leaves out, so its JPL intensity is only reproduced at
the reference temperature.  HNC-142 has the JPL partition function table as its
partition function, so its JPL intensity is reproduced at all the temperatures of
that table.

The expected values are the ARTS line strengths, line.s(T, Q(T)), that reproduce
the JPL intensities of the records, scaled with the JPL partition functions of the
entries, to 1e-12.  They were confirmed with independent code.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

# JPL line records: NH3 1(0)-0(0) and HNC 1-0
NH3 = "  572498.1599  0.0020 -1.5128 3    0.0000 12 -170021303 1 0 0       0 0 1       "
HNC = "   90663.5930  0.0090 -2.4372 2    0.0000  3 -27002 101 1           0           "

# Line strengths at the reference temperature, and at the temperatures of the JPL
# partition function table of 27002 HNC
EXPECTED = {
    "NH3-4111": {300.0: 1.7136987718879266e-40},
    "HNC-142": {
        300.0: 7.827566789366163e-40,
        225.0: 1.0429943740465154e-39,
        150.0: 1.562226504962444e-39,
        75.0: 3.109779370954281e-39,
    },
}

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

    for T, y in EXPECTED[isot].items():
        x = line.s(T, pyarts.arts.SpeciesIsotope(isot).Q(T))
        assert np.isclose(x, y, rtol=1e-12, atol=0), f"{isot} at {T} K: ARTS {x} != JPL {y}"
