"""Original ARTS2 reference values, exercised through ARTS3's Python API."""

from pathlib import Path
import re
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import numpy as np
import pyarts3 as pa

TM = pa.arts.tmatrix
# Random orientation computes Clebsch-Gordan coefficients with the Wigner tables
pa.Workspace().WignerInit()
REFERENCES = Path(__file__).resolve().parents[3] / "3rdparty" / "tmatrix"
# The original programs print Mishchenko's Stokes convention, whose U is the
# negative of ARTS's: D Z D, D = diag(1, 1, -1, 1), is the same matrix in ARTS's.
D = np.diag([1.0, 1.0, -1.0, 1.0])


def number(text):
    return float(text.replace("D", "e"))


# The original program moves its polar angles 1e-7 rad towards 90 degrees and
# its azimuths towards 180, and prints the matrices there; ARTS evaluates the
# angles it is given, so it is given the moved ones. Beta, 52 degrees, is not
# within 1e-7 rad of 90 and stays.
MOVED = np.rad2deg(1e-7)


def fixed(radius=10.0, wavelength=2 * np.pi, **kw):
    kw.setdefault("equivalent_radius", "volume")
    args = dict(
        radius=radius,
        wavelength=wavelength,
        aspect_ratio=0.5,
        refractive_real=1.5,
        refractive_imag=0.02,
        theta_incident=56.0 + MOVED,
        theta_scattered=65.0 + MOVED,
        phi_incident=114.0 + MOVED,
        phi_scattered=128.0 + MOVED,
        alpha=145.0,
        beta=52.0,
    )
    args.update(kw)
    return TM.fixed(**args)


def random(radius=1.0, quadrature=5):
    # The 19 scattering angles the original program prints
    return TM.random(
        radius=radius,
        wavelength=0.5,
        angles=19,
        aspect_ratio=2.0,
        refractive_real=1.53,
        refractive_imag=0.008,
        equivalent_radius="surface_area",
        distribution="power_law",
        b=0.1,
        gamma=0.5,
        size_quadrature=quadrature,
        lower_radius_ratio=0.89031,
        upper_radius_ratio=1.56538,
    )


def small_random():
    # Ownership and thread checks do not need the expensive reference size
    # distribution. Keep a nonspherical, absorbing particle in these checks.
    return TM.random(radius=0.1, wavelength=1.0, aspect_ratio=2.0,
                     refractive_real=1.5, refractive_imag=0.01)


def test_fixed():
    text = (
        REFERENCES / "tmatrix_ampld.ref"
    ).read_text()
    result = fixed()
    assert isinstance(result.phase, pa.arts.Muelmat)
    assert result.order == int(re.search(r"NMAX=\s*(\d+)", text)[1])
    values = re.findall(r"S\d\d=\s*([^ ]+) \+ i\*\s*([^\n]+)", text)
    actual = np.asarray(result.amplitude).ravel()
    for value, (real, imag) in zip(actual, values, strict=True):
        # One unit in the last printed decimal place, independently for each
        # real/imaginary component. Never recompute or replace reference values.
        for got, reference in [(value.real, real), (value.imag, imag)]:
            exponent = int(reference.split("D")[1])
            decimals = len(reference.split("D")[0].split(".")[1])
            np.testing.assert_allclose(
                got, number(reference), rtol=0, atol=10.0 ** (exponent - decimals)
            )
    phase = D @ np.loadtxt(text.split("PHASE MATRIX\n")[1].splitlines()[:4]) @ D
    np.testing.assert_allclose(result.phase, phase, rtol=0, atol=1e-4)
    # The returned matrices own their data; a new COMMON-block calculation
    # must not alter a previously returned result.
    saved = np.array(result.amplitude)
    scaled = fixed(radius=10e-6, wavelength=2e-6 * np.pi)
    np.testing.assert_allclose(
        np.asarray(scaled.amplitude) / 1e-6, saved, rtol=2e-6, atol=1e-8
    )
    np.testing.assert_allclose(
        np.asarray(scaled.phase) / 1e-12, result.phase, rtol=2e-6, atol=1e-6
    )
    small_random()
    np.testing.assert_array_equal(result.amplitude, saved)
    # The shape is used: the cylinder of the same aspect ratio scatters differently
    cylinder = fixed(shape="cylinder")
    assert not np.allclose(cylinder.amplitude, saved)
    np.testing.assert_array_equal(result.amplitude, saved)
    return saved


def test_random():
    text = (REFERENCES / "tmatrix_tmd.ref").read_text()
    summary = re.findall(
        r"CEXT=\s*(\S+)\s+CSCA=\s*(\S+)\s+W=\s*(\S+)\s+<COS>=\s*(\S+)", text
    )
    tables = re.findall(r"      <\s+F11[^\n]+\n((?:[^\n]+\n){19})", text)
    assert len(tables) == 2
    # The original single call scans radii 1 and .5, with 7 and 4 size
    # quadrature points. The ARTS3 interface computes one distribution/call.
    for i, (radius, quadrature) in enumerate([(1.0, 5), (0.5, 2)]):
        result = random(radius, quadrature)
        # The printed F11, F22, F33, F44, F12, F34; F34 in ARTS's convention, as D above
        expected = np.loadtxt(tables[i].splitlines())[:, 1:] * [1, 1, 1, 1, 1, -1]
        assert isinstance(result.phase, pa.arts.CompactPlanarMuelmatVector)
        # The compact F11, F12, F22, F33, F34, F44 in the printed order
        actual = np.asarray(result.phase)[:, [0, 2, 3, 5, 1, 4]]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-4)
        for name, reference in zip(
            ["extinction", "scattering", "albedo", "asymmetry"], summary[i], strict=True
        ):
            np.testing.assert_allclose(
                getattr(result, name),
                number(reference),
                rtol=0,
                atol=1e-5 if name in ["extinction", "scattering"] else 1e-6,
            )
        np.testing.assert_allclose(
            result.scattering / result.extinction, result.albedo, rtol=1e-14
        )
        np.testing.assert_allclose(result.effective_radius, radius, rtol=0, atol=5e-5)
        np.testing.assert_allclose(result.effective_variance, 0.1, rtol=0, atol=5e-5)
    # The same particle specified in metres.
    # Reuse the second reference result, including its four size nodes.
    base = result
    si = TM.random(
        radius=0.5e-6,
        wavelength=0.5e-6,
        angles=19,
        aspect_ratio=2.0,
        refractive_real=1.53,
        refractive_imag=0.008,
        equivalent_radius="surface_area",
        distribution="power_law",
        b=0.1,
        gamma=0.5,
        size_quadrature=2,
    )
    np.testing.assert_allclose(si.phase, base.phase, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(si.scattering / 1e-12, base.scattering, rtol=1e-12)
    # The third CEXT/CSCA footer in the ARTS2 reference is the old wrapper bug:
    # last-particle cross sections paired with distribution-averaged albedo.
    # Compare the correct, unchanged distribution summaries above instead.
    assert len(summary) == 3


def test_errors_and_threads():
    for kwargs in [
        {"radius": -1},
        {"accuracy": float("nan")},
        {"shape": "chebyshev"},
        {"shape": "chebyshev", "degree": 2, "aspect_ratio": 1.5},
        {"degree": 2},
        {"beta": 181},
        {"phi_incident": -1},
    ]:
        try:
            fixed(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid input accepted: {kwargs}")
    # A strongly absorbing particle overflows the Bessel functions inside it at once
    try:
        fixed(radius=1.0, refractive_imag=1e3)
    except RuntimeError as e:
        assert "converge" in str(e) and "not finite" in str(e)
    else:
        raise AssertionError("a T-matrix that is not finite did not produce a Python exception")
    # (r2 - r1) / ln(r2 / r1) = a has a root r1 in (0, a] for every b > 0: here r1 is about 3e-9 a.
    # One size node is the middle of [r1, r2 = 2 a (1 + b) - r1], a (1 + b).
    wide = TM.random(radius=0.05, wavelength=1.0, aspect_ratio=2.0, refractive_real=1.5,
                     refractive_imag=0.01, distribution="power_law", b=10.0)
    np.testing.assert_allclose(wide.effective_radius, 0.05 * 11, rtol=1e-14)
    assert np.isfinite(wide.extinction) and np.isfinite(wide.scattering)
    # but not one that double precision resolves once 1 + 2 b rounds to 1
    try:
        TM.random(radius=0.05, wavelength=1.0, aspect_ratio=2.0, refractive_real=1.5,
                  refractive_imag=0.01, distribution="power_law", b=1e-17)
    except RuntimeError as e:
        assert "power law" in str(e)
    else:
        raise AssertionError("power law without radii did not produce a Python exception")
    # r1 is about 2 a (1 + b) exp(-2 (1 + b)): 1e-260 a at b = 300, where y / x overflows,
    # and below the subnormal range at b = 400
    TM.random(radius=1e-4, wavelength=1.0, aspect_ratio=2.0, refractive_real=1.5,
              refractive_imag=0.01, distribution="power_law", b=300.0)
    try:
        TM.random(radius=1e-4, wavelength=1.0, aspect_ratio=2.0, refractive_real=1.5,
                  refractive_imag=0.01, distribution="power_law", b=400.0)
    except RuntimeError as e:
        assert "not resolvable" in str(e)
    else:
        raise AssertionError("a subnormal power-law radius did not produce a Python exception")

    # Exercise real concurrent entry with cheap particles
    expected = np.array(fixed(radius=0.1, wavelength=1.0).amplitude)
    expected_random = np.array(small_random().legendre)
    barrier = Barrier(4)

    def evaluate(i):
        barrier.wait(timeout=30)
        random_result = small_random() if i % 2 else None
        result = fixed(radius=0.1, wavelength=1.0)
        if random_result is not None:
            np.testing.assert_array_equal(random_result.legendre, expected_random)
        return np.array(result.amplitude)

    with ThreadPoolExecutor(max_workers=4) as pool:
        for actual in pool.map(evaluate, range(8)):
            np.testing.assert_array_equal(actual, expected)


def test_fixed_batch():
    args = (100e-6, 299792458 / 230e9, 1.5, 1.78, .003)
    geometries = [TM.Geometry(theta_incident=0., theta_scattered=90., phi_scattered=30.),
                  TM.Geometry(theta_incident=56., theta_scattered=65., phi_incident=114.,
                              phi_scattered=128., alpha=145., beta=52.),
                  TM.Geometry(theta_incident=90., theta_scattered=90.)]
    assert geometries[1].alpha == 145. and geometries[2].phi_incident == 0.
    batch = TM.fixed_batch(*args, geometries)
    assert len(batch) == len(geometries)
    saved = [np.array(r.phase) for r in batch]
    # New solver calls must not overwrite earlier geometries or retained batches.
    for i, g in enumerate(geometries):
        scalar = TM.fixed(*args, theta_incident=g.theta_incident, theta_scattered=g.theta_scattered,
                          phi_incident=g.phi_incident, phi_scattered=g.phi_scattered,
                          alpha=g.alpha, beta=g.beta)
        np.testing.assert_array_equal(batch[i].amplitude, scalar.amplitude)
        np.testing.assert_array_equal(batch[i].phase, scalar.phase)
        np.testing.assert_array_equal(batch[i].phase, saved[i])
    assert TM.fixed_batch(*args, []) == []
    for invalid in [[TM.Geometry(beta=181.)], [TM.Geometry(phi_incident=float('nan'))]]:
        try:
            TM.fixed_batch(*args, invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid batch geometry accepted')
    # Rows of numbers are no geometries
    try:
        TM.fixed_batch(*args, [[0., 90., 0., 30., 0., 0.]])
    except TypeError:
        pass
    else:
        raise AssertionError('a row of numbers was accepted as a geometry')


assert TM.available()
test_fixed()
test_random()
test_errors_and_threads()

test_fixed_batch()
