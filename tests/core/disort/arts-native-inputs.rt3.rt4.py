"""RT4, RT3 and VDISORT scattering inputs from ARTS-native data through pyarts3.

Exercises rt4/rt3/vdisort.scattering_optics for ARTS's Rayleigh GasScatterer
against closed forms of Rayleigh scattering (Chandrasekhar 1960).  The path
builders and the solvers are compared in vdisort-polradtran.rt3.rt4.py.  The
full C++ tests are cpp.fast.rt4-arts-test, cpp.fast.rt3-arts-test,
cpp.fast.vdisort-arts-test and cpp.fast.vdisort-arts-comparison.
"""

import numpy as np
import pyarts3 as pa

A = pa.arts
CROSS_SECTION = 1e-30  # m2
FREQ = 89e9
K_BOLTZMANN = 1.380649e-23


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


def rayleigh_m0(mo, mi):
    """m = 0 Rayleigh [I, Q] block (normalised to 1 over 4 pi), Q = I_v - I_h."""
    a, b = mo * mo, mi * mi
    return np.array(
        [
            [3 / 8 * (3 - a - b + 3 * a * b), 3 / 8 * (1 - 3 * a) * (1 - b)],
            [3 / 8 * (1 - a) * (1 - 3 * b), 9 / 8 * (1 - a) * (1 - b)],
        ]
    )


species = rayleigh_species()
atm = air(8e4, 260.0)
sigma = CROSS_SECTION * atm.pressure / (K_BOLTZMANN * atm.temperature)

# RT3: the Legendre series of Rayleigh scattering, RT3's column order
s = A.rt3.scattering_optics(species, atm, FREQ, 2, 8, 1e-12)
ref = np.array(
    [
        [1.0, -0.5, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.5, 0.0, 0.0, 1.5],
        [0.5, 0.5, 0.0, 0.0, 0.5, 0.0],
    ]
)
assert np.abs(np.asarray(s.legendre) - ref).max() < 1e-14
assert (
    abs(s.extinction - sigma) < 1e-14 * sigma
    and abs(s.scattering - sigma) < 1e-14 * sigma
)

# RT4: the azimuthal mean of the laboratory-frame phase matrix on the streams
mu = np.asarray(A.rt4.get_quadrature(8).mu)
o = A.rt4.scattering_optics(species, atm, FREQ, mu, 2, 16)
phase = np.asarray(o.phase)
dev = 0.0
for ho in (A.rt4.down, A.rt4.up):
    for hi in (A.rt4.down, A.rt4.up):
        for io, x in enumerate(mu):
            for ii, y in enumerate(mu):
                mo = x if ho == A.rt4.up else -x
                mi = y if hi == A.rt4.up else -y
                dev = max(
                    dev,
                    np.abs(
                        phase[ho, hi, io, ii] * 4 * np.pi / sigma - rayleigh_m0(mo, mi)
                    ).max(),
                )
assert dev < 1e-13, dev

# VDISORT: the m = 0 cosine coefficients on signed streams
signed = np.concatenate((mu, -mu))
f = A.vdisort.scattering_optics(species, atm, FREQ, signed, signed, 2, 16, 16, 1e-12)
c0 = np.asarray(f.cosine)[0]
dev = max(
    np.abs(c0[o_, i_, :2, :2] - rayleigh_m0(signed[o_], signed[i_])).max()
    for o_ in range(len(signed))
    for i_ in range(len(signed))
)
assert dev < 1e-13, dev
assert abs(np.asarray(f.cosine)[0, 0, 0, 3, 3] - 1.5 * signed[0] ** 2) < 1e-13
