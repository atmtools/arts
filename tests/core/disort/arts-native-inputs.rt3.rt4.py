"""RT4, RT3 and VDISORT scattering inputs from ARTS-native data through pyarts3.

Exercises rt4/rt3/vdisort.scattering_optics for ARTS's Rayleigh GasScatterer
against tabulated values of Rayleigh scattering (Chandrasekhar 1960).  The path
builders and the solvers are compared in vdisort-polradtran.rt3.rt4.py.  The
full C++ tests are cpp.fast.rt4-arts-test, cpp.fast.rt3-arts-test,
cpp.fast.vdisort-arts-test and cpp.fast.vdisort-arts-comparison.
"""

import numpy as np
import pyarts3 as pa

A = pa.arts
CROSS_SECTION = 1e-30  # m2
FREQ = 89e9


def rayleigh_species():
    species = A.ArrayOfScatteringSpecies()
    species.add(
        A.GasScatterer(
            A.ConstantGasScattering(CROSS_SECTION), A.RayleighGasScattering(0.0)
        )
    )
    return species


def air(p, t):
    a = A.AtmPoint()
    a.pressure = p
    a.temperature = t
    return a


# The m = 0 Rayleigh [I, Q] block (normalised to 1 over 4 pi, Q = I_v - I_h) for outgoing and incident |mu| in MU;
# it depends on mu^2 only
MU = np.array([0.2, 0.6, 1.0])
RAYLEIGH_M0 = np.array([
    [[[1.0968, 0.3168], [0.3168, 1.0368]], [[0.9912, 0.2112], [-0.0288, 0.6912]], [[0.78, 0.0], [-0.72, 0.0]]],
    [[[0.9912, -0.0288], [0.2112, 0.6912]], [[1.0008, -0.0192], [-0.0192, 0.4608]], [[1.02, 0.0], [-0.48, 0.0]]],
    [[[0.78, -0.72], [0.0, 0.0]], [[1.02, -0.48], [0.0, 0.0]], [[1.5, 0.0], [0.0, 0.0]]],
])  # fmt: skip


species = rayleigh_species()
atm = air(8e4, 260.0)
sigma = CROSS_SECTION * A.physics.number_density(atm.pressure, atm.temperature)

# RT3: the Legendre series of Rayleigh scattering, (F11, F12, F22, F33, F34, F44)
s = A.rt3.scattering_optics(species, atm, A.AtmSpectralPoint(), FREQ, 2, 1e-12)
ref = np.array(
    [
        [1.0, -0.5, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.5, 0.0, 1.5],
        [0.5, 0.5, 0.5, 0.0, 0.0, 0.0],
    ]
)
assert np.abs(np.asarray(s.legendre) - ref).max() < 1e-14
assert (
    abs(s.extinction - sigma) < 1e-14 * sigma
    and abs(s.scattering - sigma) < 1e-14 * sigma
)

# RT4: the azimuthal mean of the laboratory-frame phase matrix on the streams, the same in all hemispheres
o = A.rt4.scattering_optics(species, atm, A.AtmSpectralPoint(), FREQ, MU, 2)
dev = np.abs(np.asarray(o.phase) * 4 * np.pi / sigma - RAYLEIGH_M0).max()
assert dev < 1e-13, dev

# VDISORT: the m = 0 cosine coefficients on signed streams
signed = np.concatenate((MU, -MU))
f = A.vdisort.scattering_optics(species, atm, A.AtmSpectralPoint(), FREQ, signed, signed, 2, 1e-12)
c0 = np.asarray(f.cosine)[0]
dev = np.abs(c0[..., :2, :2] - np.tile(RAYLEIGH_M0, (2, 2, 1, 1))).max()
assert dev < 1e-13, dev
assert abs(c0[0, 0, 3, 3] - 0.06) < 1e-13  # 3/2 mu_out mu_in
