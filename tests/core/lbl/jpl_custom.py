"""Reading JPL line records with custom JPL species information and coded degeneracies

A line record whose species tag is not in the builtin JPL species information can
be read by giving that information explicitly, also with intensities at another
reference temperature.  Upper state degeneracies with a letter for the hundreds
are only read when asked for, and a degeneracy with an asterisk is unknown.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

# JPL 27002 HNC 1-0, and the JPL partition function of HNC at 300 K
HNC = "   90663.5930  0.0090 -2.4372 2    0.0000  3 -27002 101 1           0           "
Q300 = 138.311
ISOT = "HNC-142"
Q = pyarts.arts.SpeciesIsotope(ISOT).Q

# The ARTS line strengths, line.s(T, Q(T)), that reproduce the JPL intensity of
# the HNC record at 300 K and its JPL scaling to 150 K, confirmed with independent
# code.  The ARTS partition function of HNC is the JPL table.
S = {300.0: 7.827566789366163e-40, 150.0: 1.562226504962444e-39}

# The JPL log intensity of the HNC record at 296 K, to the 4 decimals of the file,
# and the ARTS line strength at 296 K that reproduces it
LGINT296 = -2.4256
S296 = 7.933040934457255e-40


def record(lgint=None, gup=None, tag=None):
    """The HNC record with another log intensity, degeneracy field or species tag"""
    out = HNC
    if lgint is not None:
        out = out[:21] + f"{lgint:8.4f}" + out[29:]
    if gup is not None:
        out = out[:41] + gup + out[44:]
    if tag is not None:
        out = out[:44] + f"{tag:7d}" + out[51:]
    return out


def strength(line, T):
    return line.s(T, Q(T))


def read(rec, **kwargs):
    ws = pyarts.Workspace()
    with tempfile.TemporaryDirectory() as tmp:
        fn = os.path.join(tmp, "lines.cat")
        with open(fn, "w") as fil:
            fil.write(rec + "\n")
        ws.abs_bandsReadJPL(file=fn, **kwargs)
    assert len(ws.abs_bands) == 1
    key, band = list(ws.abs_bands.items())[0]
    assert str(key.isot) == ISOT, f"{key.isot} != {ISOT}"
    return band.lines[0]


def fails(rec, **kwargs):
    try:
        read(rec, **kwargs)
    except Exception:
        return True
    return False


def info(tag, QT0=Q300, T0=300.0):
    return pyarts.arts.JplSpeciesInfo(id=tag, spec=ISOT, QT0=QT0, T0=T0)


# %% A species tag that JPL does not have is only read with custom information
wrong = record(tag=123456)
assert fails(wrong)
line = read(wrong, jpl_species=pyarts.arts.ArrayOfJplSpeciesInfo([info(123456)]))
for T in S:
    assert np.isclose(strength(line, T), S[T], rtol=1e-12, atol=0)

# %% Only the custom information is used when it is given
assert fails(HNC, jpl_species=pyarts.arts.ArrayOfJplSpeciesInfo([info(123456)]))

# %% Intensities at 296 K reproduce the 300 K intensities
rec296 = record(lgint=LGINT296)
line296 = read(
    rec296,
    jpl_species=pyarts.arts.ArrayOfJplSpeciesInfo([info(27002, Q(296.0), 296.0)]),
)
assert np.isclose(strength(line296, 296.0), S296, rtol=1e-12, atol=0)
# Only limited by the 4 decimals of the log intensity in the file
assert np.isclose(strength(line296, 300.0), S[300.0], rtol=3e-4, atol=0)

# %% Degeneracies with a letter for the hundreds are only read when asked for
lettered = record(gup="A32")
assert fails(lettered)
line = read(lettered, parse_letters=1)
assert np.isclose(line.gu, 1032, rtol=1e-12)
assert np.isclose(strength(line, 300.0), S[300.0], rtol=1e-12, atol=0)

# %% A degeneracy with an asterisk is unknown, which keeps the line strength
starred = record(gup="*32")
line = read(starred)
assert line.gu < 0 and line.a < 0
assert np.isclose(strength(line, 300.0), S[300.0], rtol=1e-12, atol=0)

# %% Custom information must be unique and positive
assert fails(HNC, jpl_species=pyarts.arts.ArrayOfJplSpeciesInfo([info(27002), info(27002)]))
assert fails(HNC, jpl_species=pyarts.arts.ArrayOfJplSpeciesInfo([info(27002, QT0=0.0)]))
