"""Azimuthally random orientation through the Python API.

ARTS's grids against the T-matrix's angles, the normalization of the tilt
weights and the rules of the arguments; tmatrix-aro-test checks the physics.
"""

import numpy as np
import pyarts3 as pa

TM = pa.arts.tmatrix
EPS = np.finfo(float).eps
ZA = np.array([10.0, 90.0, 150.0])
DAA = np.array([-120.0, 0.0, 45.0])
PARTICLE = dict(
    radius=1.5,
    wavelength=2 * np.pi,
    shape=TM.Shape.spheroid(1.6),
    refractive_real=1.5,
    refractive_imag=0.02,
)


def aro(tilts, weights, **kw):
    return TM.azimuthally_random(
        **PARTICLE,
        za_inc=ZA,
        delta_aa=DAA,
        za_scat=ZA,
        tilt_angles=tilts,
        tilt_weights=weights,
        **kw,
    )


# A vertical axis is invariant to the azimuth, so its average is the fixed
# orientation, with the incident azimuth 90 - aa = 90 degrees and the
# scattered one 90 - delta_aa
vertical = aro([0.0], [1.0])
phase = np.array(vertical.phase)
extinction = np.array(vertical.extinction)
assert phase.shape == (ZA.size, DAA.size, ZA.size, 4, 4)
assert extinction.shape == (ZA.size, 3)
assert np.array(vertical.absorption).shape == (ZA.size, 2)

geometries = [
    TM.Geometry(
        theta_incident=zi,
        theta_scattered=zs,
        phi_incident=90.0,
        phi_scattered=(90.0 - a) % 360.0,
        alpha=0.0,
        beta=0.0,
    )
    for zi in ZA
    for a in DAA
    for zs in ZA
]
fixed = TM.fixed_batch(**PARTICLE, geometries=geometries)
reference = np.array([np.array(f.phase) for f in fixed]).reshape(phase.shape)
# Rounding: a few ulp of the largest element, grown with the nmax^2 terms of
# each amplitude
bound = 16 * EPS * vertical.order * np.max(np.abs(reference))
assert np.max(np.abs(phase - reference)) <= bound, np.max(np.abs(phase - reference))

# The extinction of the forward amplitude, lambda (Im(S11 + S22),
# Im(S11 - S22), Re(S11 - S22)) in ARTS's Stokes convention
forward = TM.fixed_batch(
    **PARTICLE,
    geometries=[
        TM.Geometry(
            theta_incident=z,
            theta_scattered=z,
            phi_incident=90.0,
            phi_scattered=90.0,
            alpha=0.0,
            beta=0.0,
        )
        for z in ZA
    ],
)
for k, f in zip(extinction, forward):
    S = np.array(f.amplitude)
    K = PARTICLE["wavelength"] * np.array(
        [(S[0, 0] + S[1, 1]).imag, (S[0, 0] - S[1, 1]).imag, (S[0, 0] - S[1, 1]).real]
    )
    assert np.max(np.abs(k - K)) <= 16 * EPS * vertical.order * K[0], (k, K)

# The weights are normalized: 2 and 6 are 1/4 and 3/4, exactly
for a, b in [
    (aro([30.0, 70.0], [2.0, 6.0]), aro([30.0, 70.0], [0.25, 0.75])),
]:
    for field in ["phase", "extinction", "absorption"]:
        assert np.array_equal(np.array(getattr(a, field)), np.array(getattr(b, field))), field

# The rules of the arguments
for tilts, weights, grids in [
    ([30.0, 70.0], [1.0], {}),
    ([], [], {}),
    ([30.0], [-1.0], {}),
    ([30.0, 70.0], [0.0, 0.0], {}),
    ([181.0], [1.0], {}),
    ([30.0], [np.nan], {}),
    ([30.0], [1.0], {"za_inc": [181.0]}),
    ([30.0], [1.0], {"za_scat": []}),
    ([30.0], [1.0], {"delta_aa": [np.inf]}),
]:
    args = dict(
        **PARTICLE,
        za_inc=ZA,
        delta_aa=DAA,
        za_scat=ZA,
        tilt_angles=tilts,
        tilt_weights=weights,
    )
    args.update(grids)
    try:
        TM.azimuthally_random(**args)
    except Exception as e:
        assert "tilt" in str(e) or "angle" in str(e) or "delta_aa" in str(e), e
    else:
        raise AssertionError(f"no error for {tilts}, {weights}, {grids}")

try:
    aro([30.0, 70.0], [1.0])
except Exception as e:
    assert "tilt_angles has 2 elements and tilt_weights 1" in str(e), e
