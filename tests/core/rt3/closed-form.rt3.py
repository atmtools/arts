"""RT3 through pyarts3.arts.rt3 against references from outside RT3.

* The quadrature: Gauss-Legendre from numpy, Lobatto from a table checked
  with mpmath.
* Single scattering by an optically thin Rayleigh layer: a table of the
  exact single-scattering Stokes vectors on the two extra streams.
* A gas-only atmosphere over a Fresnel surface: ARTS's linear-in-tau
  transmission (rtepack.tran), Planck (physics.planck) and Fresnel
  reflection (vdisort.fresnel_reflection).

No RT3 output enters a reference.
"""

from pathlib import Path

import numpy as np
from pyarts3 import arts

rt3 = arts.rt3
polradtran = arts.polradtran

C = arts.constants.c

RUNMIETEST = Path(__file__).resolve().parents[3] / "3rdparty/polradtran/runmietest"


def heredoc(text, name):
    """The body of 'cat >name <<EOF ... EOF' in a csh script."""
    return text.split(f"cat >{name} <<EOF\n")[1].split("\nEOF")[0]


def from_rt3_file(columns):
    """The [nleg + 1, 6] series of an RT3 scattering file, columns (F11, F12,
    F33, F34, F22, F44), in rtepack's order (F11, F12, F22, F33, F34, F44)."""
    return np.asarray(columns)[:, [0, 1, 4, 2, 3, 5]]


# Positive half of the 2 nmu-point Gauss-Lobatto rule on [-1, 1], (mu, weights)
LOBATTO = {
    1: ([1.0], [1.0]),
    2: ([0.44721359549995794, 1.0], [0.83333333333333333, 0.16666666666666667]),
    5: (
        [0.16527895766638702, 0.4779249498104445, 0.73877386510550508, 0.91953390816645881, 1.0],
        [0.32753976118389746, 0.29204268367968376, 0.22488934206312645, 0.13330599085107011, 0.022222222222222222],
    ),
    8: (
        [0.10132627352194945, 0.29983046890076321, 0.48605942188713761, 0.65238870288249309,
         0.79200829186181506, 0.89920053309347209, 0.96956804627021793, 1.0],
        [0.20195830817822987, 0.19369002382520358, 0.17749191339170413, 0.15402698080716428,
         0.1242553821325141, 0.089393697325930801, 0.050850361005919905, 0.0083333333333333333],
    ),
}  # fmt: skip


def expected_quadrature(nmu, quad):
    """The rules as documented: numpy's Gauss-Legendre, the Lobatto table."""
    if quad == polradtran.QuadratureType.double_gauss:
        x, w = np.polynomial.legendre.leggauss(nmu)
        return 0.5 * (x + 1.0), 0.5 * w
    if quad == polradtran.QuadratureType.gauss:
        x, w = np.polynomial.legendre.leggauss(2 * nmu)
        return x[nmu:], w[nmu:]
    return np.array(LOBATTO[nmu])


def test_quadrature():
    for quad in polradtran.QuadratureType:
        for nmu in (1, 2, 5, 8):
            q = polradtran.get_quadrature(nmu, quad)
            mu, w = expected_quadrature(nmu, quad)
            np.testing.assert_allclose(np.asarray(q.mu), mu, rtol=0, atol=1e-14)
            np.testing.assert_allclose(np.asarray(q.weights), w, rtol=0, atol=1e-14)
        print(f"quadrature {quad.name:13s} == reference for nmu = 1, 2, 5, 8")
    assert rt3.max_legendre_degree(8, polradtran.QuadratureType.gauss) == 29
    assert rt3.max_legendre_degree(8, polradtran.QuadratureType.double_gauss) == 13
    assert rt3.max_legendre_degree(8, polradtran.QuadratureType.lobatto) == 27


def test_single_scattering():
    """A conservative Rayleigh layer of tau = 1e-6 over a black surface, lit
    by a beam of flux F on the horizontal at mu0 = 0.6 propagating toward
    phi = 0: exact single scattering, up at the top and down at the bottom,
    on the extra streams mu = 0.5 and 0.77.  The neglected multiple
    scattering is of relative order tau."""
    tau, mu0, flux = 1e-6, 0.6, 2.5
    rayleigh = np.array(  # F11, F12, F22, F33, F34, F44
        [
            [1.0, -0.5, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.5, 0.0, 1.5],
            [0.5, 0.5, 0.5, 0.0, 0.0, 0.0],
        ]
    )
    p = rt3.Problem(
        nstokes=4,
        nmu=8,
        extra_mu=[0.5, 0.77],
        aziorder=4,
        max_delta_tau=1e-9,
        direct_flux=flux,
        direct_mu=mu0,
        thermal=False,
        frequency=C / 0.5e-6,
        height=np.array([1.0, 0.0]),
        temperature=np.zeros(2),
        gas_extinction=np.zeros(1),
        scattering_sets=[
            rt3.ScatteringSet(extinction=tau, scattering=tau, legendre=rayleigh)
        ],
        layer_scattering_index=[0],
        ground=polradtran.LambertianSurface(0.0),
    )
    r = rt3.solve(p)
    phi = np.radians([0.0, 40.0, 90.0, 160.0, 200.0, 300.0])
    up = np.asarray(rt3.azimuth_radiance(r.up, phi))[0]
    down = np.asarray(rt3.azimuth_radiance(r.down, phi))[1]

    # [phi, mu, (I, Q, U)] in units of 1e-7; V is 0
    ref_up = 1e-7 * np.array([
        [[5.7410455, -4.2061202, 0.0000000], [3.2371769, -3.2220239, 0.0000000]],
        [[5.2383608, -2.0784502, 4.2252678], [3.2458735, -1.5053042, 2.8389314]],
        [[5.4212053, 1.8402257, 4.1349591], [3.9189392, 1.5936269, 1.9782017]],
        [[9.4720570, 0.2695936, 0.3912129], [6.0933151, 0.1176874, -0.3464419]],
        [[9.4720570, 0.2695936, -0.3912129], [6.0933151, 0.1176874, 0.3464419]],
        [[4.9842955, -0.1882307, -4.9592994], [3.3676948, 0.0089104, -3.0914931]],
    ])
    ref_down = 1e-7 * np.array([
        [[9.8760046, -0.0711611, 0.0000000], [6.2836075, -0.1755932, 0.0000000]],
        [[8.4059233, 1.0891123, 1.0905332], [5.5795747, 0.8283971, -0.2958043]],
        [[5.4212053, 1.8402257, 4.1349591], [3.9189392, 1.5936269, 1.9782017]],
        [[5.5864664, -3.6159970, 2.4372657], [3.2306067, -2.7450210, 1.6996116]],
        [[5.5864664, -3.6159970, -2.4372657], [3.2306067, -2.7450210, -1.6996116]],
        [[7.0517750, 1.8792489, -2.2026599], [4.8909101, 1.5321257, -0.3348527]],
    ])  # fmt: skip
    extra = slice(p.nmu, None)
    assert np.array_equal(np.asarray(r.mu)[extra], [0.5, 0.77])
    up, down = up[:, extra], down[:, extra]
    ref_up, ref_down = (np.pad(x, ((0, 0), (0, 0), (0, 1))) for x in (ref_up, ref_down))
    scale = max(np.abs(ref_up).max(), np.abs(ref_down).max())
    dev = max(np.abs(up - ref_up).max(), np.abs(down - ref_down).max()) / scale
    print(
        f"thin Rayleigh layer vs table: max |dev| / max |I| = {dev:.2e} (tau = {tau:.0e})"
    )
    assert dev < 10 * tau, dev


def test_gas_fresnel():
    """Gas-only, Planck linear in tau, polarizing Fresnel surface, nstokes 4;
    every level and stream, m = 0, and every m > 0 mode exactly 0.

    Each layer is ARTS's linear-in-tau step, I_end = T (I_start - J_start) +
    L (J_start - J_end) + J_end with T and L of rtepack.tran and J = [B, 0, 0, 0].
    The surface emits (1 - R) [B_s, 0, 0, 0] and reflects R I_down, R of
    vdisort.fresnel_reflection."""
    f, n = 89e9, 3.0 + 0.2j
    height = np.array([3000.0, 2000.0, 1000.0, 0.0])
    temperature = np.array([220.0, 240.0, 265.0, 285.0])
    gas = np.array([1e-4, 3e-4, 5e-4])
    sky, tsurf = 2.725, 290.0
    p = rt3.Problem(
        nstokes=4,
        nmu=8,
        extra_mu=[0.45, 1.0],
        aziorder=2,
        frequency=f,
        height=height,
        temperature=temperature,
        gas_extinction=gas,
        layer_scattering_index=[-1, -1, -1],
        sky_temperature=sky,
        surface_temperature=tsurf,
        ground=polradtran.FresnelSurface(n),
    )
    r = rt3.solve(p)
    mu = np.asarray(r.mu)
    j = np.zeros((4, 4))
    j[:, 0] = arts.physics.planck(f, temperature)
    dz = np.abs(np.diff(height))

    def step(i0, l, m, start, end):
        t = arts.rtepack.tran(arts.Propmat(gas[l]), arts.Propmat(gas[l]), dz[l] / m)
        return np.asarray(t()) @ (i0 - j[start]) + np.asarray(t.linsrc()) @ (j[start] - j[end]) + j[end]

    down = np.zeros((4, len(mu), 4))
    up = np.zeros((4, len(mu), 4))
    for i, m in enumerate(mu):
        down[0, i, 0] = arts.physics.planck(f, sky)
        for l in range(3):
            down[l + 1, i] = step(down[l, i], l, m, l, l + 1)
        refl = np.asarray(arts.vdisort.fresnel_reflection(m, n))
        up[3, i] = (np.eye(4)[0] - refl[:, 0]) * arts.physics.planck(f, tsurf) + refl @ down[3, i]
        for l in reversed(range(3)):
            up[l, i] = step(up[l + 1, i], l, m, l + 1, l)

    got_up, got_down = np.asarray(r.up), np.asarray(r.down)
    assert np.all(got_up[:, 1:] == 0.0) and np.all(got_down[:, 1:] == 0.0)
    dev = max(
        np.max(np.abs(got_up[:, 0] - up) / up[..., :1]),
        np.max(np.abs(got_down[:, 0] - down) / down[..., :1]),
    )
    assert np.all(
        got_up[-1, 0, mu < 0.9, 1] > 0.0
    )  # warm dielectric: vertically polarized
    print(f"gas-only + Fresnel {n}, nstokes 4: max rel dev {dev:.2e}")
    assert dev < 1e-12, dev


def expect_error(what, p):
    """solve() must raise."""
    try:
        rt3.solve(p)
    except RuntimeError:
        print(f"error path: {what:44s} raises")
        return
    raise AssertionError(f"{what} did not raise")


def test_errors():
    base = dict(
        nstokes=4,
        nmu=4,
        aziorder=2,
        direct_flux=1.0,
        direct_mu=0.5,
        frequency=C / 3e-6,
        height=[1.0, 0.0],
        temperature=[250.0, 260.0],
        gas_extinction=[0.0],
        layer_scattering_index=[0],
        ground=polradtran.LambertianSurface(0.1),
    )
    mie = from_rt3_file(
        [
            [float(x) for x in line.split()[1:]]
            for line in heredoc(RUNMIETEST.read_text(), "mietest.sca").splitlines()[4:]
        ]
    )
    ok = rt3.ScatteringSet(1.0, 0.9, mie)
    rt3.solve(rt3.Problem(**base, scattering_sets=[ok]))

    expect_error(
        "nstokes = 5",
        rt3.Problem(**{**base, "nstokes": 5}, scattering_sets=[ok]),
    )
    expect_error(
        "beam over a Fresnel surface",
        rt3.Problem(
            **{**base, "ground": polradtran.FresnelSurface(1.5)}, scattering_sets=[ok]
        ),
    )
    expect_error(
        "Mie series of degree 11, nmu = 2",
        rt3.Problem(**{**base, "nmu": 2}, scattering_sets=[ok]),
    )
    expect_error(
        "unnormalised phase function",
        rt3.Problem(**base, scattering_sets=[rt3.ScatteringSet(1.0, 0.9, 0.5 * mie)]),
    )
    # Extra angles with every quadrature (RT3's 'E' type had them with gauss only)
    r = rt3.solve(
        rt3.Problem(**base, quad=polradtran.QuadratureType.lobatto, extra_mu=[0.5], scattering_sets=[ok])
    )
    assert np.asarray(r.mu)[-1] == 0.5 and np.all(np.isfinite(np.asarray(r.up)))
    print("extra_mu with lobatto runs")


test_quadrature()
test_single_scattering()
test_gas_fresnel()
test_errors()
