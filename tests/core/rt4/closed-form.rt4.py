"""RT4 through pyarts3.arts.rt4 against references from outside RT4.

* The quadrature: Gauss-Legendre from numpy, Lobatto from a table checked
  with mpmath.
* Gas-only layers over a Fresnel surface: ARTS's linear-in-tau transmission
  (rtepack.tran), Planck (physics.planck) and Fresnel reflection
  (vdisort.fresnel_reflection).
* A particle layer with lower-triangular K: a table of the exact solution.
* Isothermal Kirchhoff with ARTS's Rayleigh phase matrix: I = B, Q = 0.

No RT4 output enters a reference.
"""

import numpy as np
from pyarts3 import arts

rt4 = arts.rt4
polradtran = arts.polradtran

FREQUENCY = 50e9
DOWN, UP = rt4.down, rt4.up
assert (DOWN, UP) == (0, 1)

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


HEIGHT = np.array([3000.0, 2000.0, 1000.0, 0.0])
TEMPERATURE = np.array([220.0, 240.0, 265.0, 285.0])
GAS = np.array([1e-4, 3e-4, 5e-4])  # vertical tau 0.1, 0.3, 0.5
SKY = 2.725
SURFACE = 290.0
INDEX = 3.0 + 0.2j


def planck(t):
    return arts.physics.planck(FREQUENCY, t)


def gas_references(mu):
    """Down and up [level, mu, 2] for the 3 gas-only layers over the Fresnel surface.

    Each layer is ARTS's linear-in-tau step, I_end = T (I_start - J_start) +
    L (J_start - J_end) + J_end with T and L of rtepack.tran and J = [B, 0].
    The surface emits (1 - R) [B_s, 0] and reflects R [I_d, Q_d], R the
    [I, Q] block of vdisort.fresnel_reflection.
    """
    nlay = len(GAS)
    dz = np.abs(np.diff(HEIGHT))
    j = np.column_stack([planck(TEMPERATURE), np.zeros(nlay + 1)])

    def step(i0, l, m, start, end):
        t = arts.rtepack.tran(arts.Propmat(GAS[l]), arts.Propmat(GAS[l]), dz[l] / m)
        T, L = np.asarray(t())[:2, :2], np.asarray(t.linsrc())[:2, :2]
        return T @ (i0 - j[start]) + L @ (j[start] - j[end]) + j[end]

    down = np.zeros((nlay + 1, len(mu), 2))
    up = np.zeros((nlay + 1, len(mu), 2))
    for i, m in enumerate(mu):
        down[0, i] = [planck(SKY), 0.0]
        for l in range(nlay):
            down[l + 1, i] = step(down[l, i], l, m, l, l + 1)
        refl = np.asarray(arts.vdisort.fresnel_reflection(m, INDEX))[:2, :2]
        up[nlay, i] = (np.array([1.0, 0.0]) - refl[:, 0]) * planck(SURFACE) + refl @ down[nlay, i]
        for l in reversed(range(nlay)):
            up[l, i] = step(up[l + 1, i], l, m, l + 1, l)
    return down, up


def max_rel_dev(r, down, up):
    """max over levels and streams of |dI| / I_ref and |dQ| / I_ref."""
    dev = []
    for got, ref in ((np.asarray(r.down), down), (np.asarray(r.up), up)):
        dev.append(np.max(np.abs(got - ref) / ref[..., :1], axis=(0, 1)))
    return np.max(dev, axis=0)


def problem(**kw):
    """The 3-layer atmosphere over the Fresnel surface; numpy and list inputs."""
    args = {
        "nstokes": 2,
        "nmu": 8,
        "extra_mu": np.array([0.45, 1.0]),
        "frequency": FREQUENCY,
        "height": HEIGHT,
        "temperature": TEMPERATURE,
        "gas_extinction": GAS,
        "layer_optics_index": [-1, -1, -1],
        "sky_temperature": SKY,
        "surface_temperature": SURFACE,
        "ground": polradtran.FresnelSurface(INDEX),
    }
    args.update(kw)
    return rt4.Problem(**args)


def test_quadrature():
    """RT4's rules equal numpy's Gauss-Legendre and the Lobatto table, ascending."""
    for quad in polradtran.QuadratureType:
        for nmu in (1, 2, 5, 8):
            q = polradtran.get_quadrature(nmu, quad)
            mu, w = expected_quadrature(nmu, quad)
            np.testing.assert_allclose(
                np.asarray(q.mu), mu, rtol=0, atol=1e-14, err_msg=str(quad)
            )
            np.testing.assert_allclose(
                np.asarray(q.weights), w, rtol=0, atol=1e-14, err_msg=str(quad)
            )
        print(f"quadrature {quad.name:13s} == reference for nmu = 1, 2, 5, 8")


def test_gas_fresnel():
    """Gas-only, Planck linear in tau, polarizing Fresnel surface; every level and stream.

    RT4 solves gas-only layers analytically, so only round-off remains.
    """
    for quad in polradtran.QuadratureType:
        p = problem(quad=quad)
        r = rt4.solve(p)
        mu = np.asarray(r.mu)
        mu_q, w_q = expected_quadrature(p.nmu, quad)
        np.testing.assert_allclose(
            mu, np.concatenate([mu_q, [0.45, 1.0]]), rtol=0, atol=1e-14
        )
        np.testing.assert_allclose(
            np.asarray(r.weights), np.concatenate([w_q, [0.0, 0.0]]), rtol=0, atol=1e-14
        )

        down, up = gas_references(mu)
        # Oblique emission from a warm dielectric is vertically polarized
        assert np.all(np.asarray(r.up)[-1, mu < 0.9, 1] > 0.0) and np.all(
            up[-1, mu < 0.9, 1] > 0.0
        )
        assert np.all(np.asarray(r.down)[..., 1] == 0.0)
        dev = max_rel_dev(r, down, up)
        print(
            f"gas-only + Fresnel {INDEX}, {quad.name:13s} max rel dev I {dev[0]:.2e}, Q {dev[1]:.2e}"
        )
        assert np.all(dev < 1e-12), dev


def test_layout():
    """A zero-phase particle layer with stream-dependent, lower-triangular K and
    Stokes-asymmetric absorption, over the Fresnel surface.

    The phase matrix is a non-zero forward-only matrix that scatters nothing
    in net (it is offset by K), so it must be put in the forward (same
    hemisphere, same stream) quadrant: phase[h, h, i, i] = s / (2 pi w_i),
    K11 = k + s.  Any hemisphere, stream or Stokes transpose of the numpy
    arrays changes the answer.  The doubling is first order in the initial
    sublayer, whose largest slant thickness is max_delta_tau / mu_min.

    The reference is the exact solution of the transfer equation, tabulated
    in units of B(SURFACE) on the 8 double-Gauss streams and the two extra ones.
    """
    p = problem(max_delta_tau=1e-7, layer_optics_index=[-1, 0, -1])
    q = polradtran.get_quadrature(p.nmu, p.quad)
    mu = np.concatenate([np.asarray(q.mu), np.asarray(p.extra_mu)])
    w = np.concatenate([np.asarray(q.weights), np.zeros(len(p.extra_mu))])
    nmu = len(mu)

    kp = 2e-4 * (1.0 + 0.5 * mu)
    kappa = 0.6e-4 * (1.0 + mu)
    a_i = kp.copy()
    a_q = 0.3e-4 * (1.0 + 2.0 * mu)
    s = 1e-4 * (1.0 + mu)  # forward "scattering" on the quadrature streams only
    s[w == 0.0] = 0.0

    ext = np.zeros((2, nmu, 2, 2))
    absorption = np.zeros((2, nmu, 2))
    phase = np.zeros((2, 2, nmu, nmu, 2, 2))
    for h in (DOWN, UP):
        ext[h, :, 0, 0] = ext[h, :, 1, 1] = kp + s
        ext[h, :, 1, 0] = kappa
        absorption[h, :, 0] = a_i
        absorption[h, :, 1] = a_q
        for i in range(nmu):
            if w[i] > 0:
                phase[h, h, i, i, 0, 0] = phase[h, h, i, i, 1, 1] = s[i] / (
                    2.0 * np.pi * w[i]
                )
    p.optics = [rt4.LayerOptics(extinction=ext, absorption=absorption, phase=phase)]

    r = rt4.solve(p)

    unit = planck(SURFACE)
    down = np.zeros((4, nmu, 2))
    up = np.zeros((4, nmu, 2))
    down[0, :, 0] = planck(SKY)
    down[1:, :, 0] = unit * np.array([
        [0.808327951, 0.501675448, 0.277190572, 0.177018956, 0.128306931,
         0.102597682, 0.088790320, 0.082212476, 0.162831251, 0.080770376],
        [0.910012128, 0.894147565, 0.818092644, 0.692668456, 0.585919243,
         0.512440550, 0.467187994, 0.444078159, 0.665455856, 0.438872318],
        [0.979937217, 0.968568018, 0.942236070, 0.877881060, 0.795249551,
         0.723702821, 0.673806953, 0.646631964, 0.859245100, 0.640353327],
    ])
    down[2:, :, 1] = unit * np.array([
        [-0.053967400, -0.048839283, -0.022604802, 0.003848349, 0.018881263,
         0.026697383, 0.030699487, 0.032532238, 0.008174764, 0.032926507],
        [0.000000000, -0.000357188, -0.002747075, 0.001130887, 0.008110601,
         0.013860689, 0.017595789, 0.019532952, 0.002691075, 0.019970936],
    ])
    up[..., 0] = unit * np.array([
        [0.771303854, 0.808140511, 0.839079722, 0.863775394, 0.878606859,
         0.885117421, 0.887125664, 0.887424366, 0.868068560, 0.887415141],
        [0.830295973, 0.844101013, 0.864905206, 0.884018649, 0.894686965,
         0.898233072, 0.898373250, 0.897710828, 0.887255637, 0.897484786],
        [0.916185971, 0.927411786, 0.941932550, 0.947937140, 0.943837948,
         0.936019337, 0.928931322, 0.924554177, 0.947675112, 0.923491020],
        [0.982455267, 0.982370324, 0.979660243, 0.964210611, 0.944629047,
         0.928345324, 0.917064639, 0.910860612, 0.959692617, 0.909415676],
    ])
    up[..., 1] = unit * np.array([
        [-0.000325071, -0.019252087, -0.031446592, -0.029912729, -0.023418950,
         -0.019286350, -0.018759093, -0.019613032, -0.028530227, -0.019941128],
        [-0.050038078, -0.051481634, -0.047933551, -0.038214402, -0.027730840,
         -0.021988055, -0.020967967, -0.021719706, -0.035629941, -0.022038355],
        [0.000000000, 0.000068069, 0.002319821, 0.009922881, 0.016241719,
         0.014984427, 0.009307062, 0.004379691, 0.011878847, 0.003050901],
        [0.001965771, 0.009307195, 0.019089065, 0.033767034, 0.037810288,
         0.028861840, 0.016238091, 0.007294398, 0.036084752, 0.005030086],
    ])  # fmt: skip
    dev = max_rel_dev(r, down, up)
    tol = p.max_delta_tau / mu[0]
    print(
        f"layout, lower-triangular K(mu) + forward phase   max rel dev I {dev[0]:.2e}, Q {dev[1]:.2e} (tol {tol:.1e})"
    )
    assert np.all(dev < tol), dev


def test_kirchhoff(reciprocal):
    """Isothermal sky, layers and surface: I = B, Q = 0 exactly on the streams
    when a(h, mu_j) = K [1, 0] - 2 pi sum_{i, h'} w_i Z(h, j <- h', i) [1, 0].

    The non-reciprocal variant multiplies the Q row of Z (Z_QI, Z_QQ) by
    (1 + 0.3 mu_out), which makes a_Q non-zero and catches a full transpose
    of Z that the reciprocal Rayleigh matrix is blind to.  Z_II stays
    reciprocal, so the medium conserves energy (rt4.solve checks it); a
    non-reciprocal Z_II could not obey Kirchhoff's law and energy
    conservation with one absorption.  Z is rt4.scattering_optics's phase
    matrix of ARTS's Rayleigh GasScatterer, scaled to the scattering
    coefficient sigma.  Only round-off remains, which grows with the
    number of sublayers 2^n.
    """
    t = 260.0
    sigma, kabs, gas = 6e-4, 4e-4, np.array([5e-5, 2e-4])
    p = problem(
        height=np.array([2000.0, 1000.0, 0.0]),
        temperature=np.full(3, t),
        gas_extinction=gas,
        layer_optics_index=[0, -1],
        sky_temperature=t,
        surface_temperature=t,
        extra_mu=[1.0],
        max_delta_tau=1e-7,
    )
    q = polradtran.get_quadrature(p.nmu, p.quad)
    mu = np.concatenate([np.asarray(q.mu), [1.0]])
    w = np.concatenate([np.asarray(q.weights), [0.0]])
    nmu = len(mu)

    species = arts.ArrayOfScatteringSpecies()
    species.add(arts.GasScatterer(arts.ConstantGasScattering(1e-30), arts.RayleighGasScattering(0.0)))
    atm = arts.AtmPoint()
    atm.pressure, atm.temperature = 1e5, t
    rayleigh = rt4.scattering_optics(species, atm, FREQUENCY, mu, 2)
    phase = sigma / np.asarray(rayleigh.extinction)[0, 0, 0, 0] * np.asarray(rayleigh.phase)
    if not reciprocal:
        phase[..., 1, :] *= (1.0 + 0.3 * mu)[:, None, None]  # [h, h', i_out, i_in, s, s']

    ext = np.zeros((2, nmu, 2, 2))
    ext[:, :, 0, 0] = ext[:, :, 1, 1] = sigma + kabs
    absorption = np.zeros((2, nmu, 2))
    absorption[:, :, 0] = sigma + kabs
    absorption -= 2.0 * np.pi * np.einsum("i,hgjis->hjs", w, phase[..., 0])
    if reciprocal:
        assert np.max(np.abs(absorption[..., 1])) < 1e-12 * sigma
    p.optics = [rt4.LayerOptics(ext, absorption, phase)]

    r = rt4.solve(p)
    b = planck(t)
    out = np.stack([np.asarray(r.down), np.asarray(r.up)])
    dev_i = np.max(np.abs(out[..., 0] - b)) / b
    dev_q = np.max(np.abs(out[..., 1])) / b
    tau = (sigma + kabs + gas[0]) * 1000.0
    tol = 2.0 ** (np.floor(np.log2(tau / p.max_delta_tau)) + 1) * np.finfo(float).eps
    name = "Rayleigh" if reciprocal else "non-reciprocal Rayleigh"
    print(
        f"isothermal Kirchhoff, {name:23s} + gas + Fresnel max rel dev I {dev_i:.2e}, Q {dev_q:.2e} (tol {tol:.1e})"
    )
    assert dev_i < tol and dev_q < tol, (dev_i, dev_q)


def expect_error(what, p):
    """solve() must raise."""
    try:
        rt4.solve(p)
    except RuntimeError:
        print(f"error path: {what:40s} raises")
        return
    raise AssertionError(f"{what} did not raise")


def test_errors():
    expect_error("nstokes = 3", problem(nstokes=3))
    # Evans' fixed array sizes are not limits of the port: 2 * (31 + 2) > 64
    r = rt4.solve(problem(nmu=31))
    assert np.asarray(r.up).shape == (4, 33, 2) and np.all(np.isfinite(np.asarray(r.up)))
    print("nstokes * nmu_total = 66 > 64 runs")

    p = problem(layer_optics_index=[-1, 0, -1])
    nmu = p.nmu + len(p.extra_mu)
    ext = np.zeros((2, nmu, 2, 2))
    ext[DOWN, :, 0, 0] = ext[DOWN, :, 1, 1] = 1e-4
    ext[UP, :, 0, 0] = ext[UP, :, 1, 1] = 2e-4
    p.optics = [
        rt4.LayerOptics(ext, np.zeros((2, nmu, 2)), np.zeros((2, 2, nmu, nmu, 2, 2)))
    ]
    expect_error("up/down asymmetric extinction", p)

    p.optics = [
        rt4.LayerOptics(
            np.zeros((2, nmu, 2, 2)),
            np.zeros((2, nmu, 2)),
            np.zeros((2, 2, nmu, nmu, 1, 1)),
        )
    ]
    expect_error("wrong phase shape", p)


test_quadrature()
test_gas_fresnel()
test_layout()
test_kirchhoff(reciprocal=True)
test_kirchhoff(reciprocal=False)
test_errors()
