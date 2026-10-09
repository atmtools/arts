"""RT4 through pyarts3.arts.rt4 against closed forms derived here.

Every reference is computed below with numpy only (Planck from the exact SI
constants, Fresnel from Snell's law, the transfer equation integrated
analytically along each stream).  No RT4 output and no ARTS physics function
enters a reference.
"""

import numpy as np
from pyarts3 import arts

rt4 = arts.rt4
polradtran = arts.polradtran
assert rt4.available(), "this test is only collected with ENABLE_RT4=ON"

H = 6.62607015e-34
C = 299792458.0
K = 1.380649e-23

FREQUENCY = 50e9
DOWN, UP = rt4.down, rt4.up
assert (DOWN, UP) == (0, 1)


def planck(f, t):
    """B_nu(T) in W m-2 Hz-1 sr-1."""
    return 2.0 * H * f**3 / C**2 / np.expm1(H * f / (K * t))


def fresnel_vh(n, mu):
    """|r_v|^2 and |r_h|^2 for a medium of index 1 above index n."""
    cos_t = np.sqrt(1.0 - (1.0 - mu**2) / n**2 + 0j)
    rv = (n * mu - cos_t) / (n * mu + cos_t)
    rh = (mu - n * cos_t) / (mu + n * cos_t)
    return np.abs(rv) ** 2, np.abs(rh) ** 2


def own_quadrature(nmu, quad):
    """The rules as documented, from numpy's Gauss-Legendre and Legendre roots."""
    if quad == polradtran.QuadratureType.double_gauss:
        x, w = np.polynomial.legendre.leggauss(nmu)
        return 0.5 * (x + 1.0), 0.5 * w
    if quad == polradtran.QuadratureType.gauss:
        x, w = np.polynomial.legendre.leggauss(2 * nmu)
        return x[nmu:], w[nmu:]
    n = 2 * nmu  # Lobatto: +-1 and the roots of P'_{n-1}
    p = np.polynomial.legendre.Legendre.basis(n - 1)
    x = np.sort(np.concatenate([[-1.0, 1.0], p.deriv().roots().real]))
    w = 2.0 / (n * (n - 1) * p(x) ** 2)
    return x[nmu:], w[nmu:]


def propagate(v0, length, k, kappa, a_i, a_q, b_start, b_end):
    """Exact solution of dV/ds = -K V + a B(s) along a path of the given length.

    K = [[k, 0], [kappa, k]], a = [a_i, a_q], B linear in s from b_start to
    b_end.  With exp(-K u) = exp(-k u) (1 - kappa u N), N = [[0, 0], [1, 0]],
    and E_n = int_0^L u^n exp(-k u) du:
      I = I0 e^{-kL} + a_i (B_e E0 - g E1)
      Q = e^{-kL} (Q0 - kappa L I0) + a_q (B_e E0 - g E1) - kappa a_i (B_e E1 - g E2)
    with g = (B_e - B_s) / L.
    """
    x = k * length
    ex = np.exp(-x)
    e0 = -np.expm1(-x) / k
    e1 = (1.0 - (1.0 + x) * ex) / k**2
    e2 = (2.0 - (2.0 + 2.0 * x + x * x) * ex) / k**3
    g = (b_end - b_start) / length
    i = v0[0] * ex + a_i * (b_end * e0 - g * e1)
    q = (
        ex * (v0[1] - kappa * length * v0[0])
        + a_q * (b_end * e0 - g * e1)
        - kappa * a_i * (b_end * e1 - g * e2)
    )
    return np.array([i, q])


HEIGHT = np.array([3000.0, 2000.0, 1000.0, 0.0])
TEMPERATURE = np.array([220.0, 240.0, 265.0, 285.0])
GAS = np.array([1e-4, 3e-4, 5e-4])  # vertical tau 0.1, 0.3, 0.5
SKY = 2.725
SURFACE = 290.0
INDEX = 3.0 + 0.2j


def references(mu, stream_optics, surface_n):
    """Closed-form down and up [level, mu, 2] for the 3-layer atmosphere.

    stream_optics(l, i) gives (k, kappa, a_i, a_q) of layer l at stream i,
    identical in both hemispheres.  The surface is Fresnel with index
    surface_n: I_v and I_h each reflect specularly, so in (I, Q) the surface
    returns (1 - R) [B_s, 0] + R [I_d, Q_d] with R = [[R1, R2], [R2, R1]].
    """
    nlay = len(GAS)
    dz = np.abs(np.diff(HEIGHT))
    b = planck(FREQUENCY, TEMPERATURE)
    bs = planck(FREQUENCY, SURFACE)
    down = np.zeros((nlay + 1, len(mu), 2))
    up = np.zeros((nlay + 1, len(mu), 2))
    for i, m in enumerate(mu):
        down[0, i] = [planck(FREQUENCY, SKY), 0.0]
        for l in range(nlay):
            down[l + 1, i] = propagate(
                down[l, i], dz[l] / m, *stream_optics(l, i), b[l], b[l + 1]
            )
        rv, rh = fresnel_vh(surface_n, m)
        r1, r2 = 0.5 * (rv + rh), 0.5 * (rv - rh)
        refl = np.array([[r1, r2], [r2, r1]])
        up[nlay, i] = (np.array([1.0, 0.0]) - refl[:, 0]) * bs + refl @ down[nlay, i]
        for l in reversed(range(nlay)):
            up[l, i] = propagate(
                up[l + 1, i], dz[l] / m, *stream_optics(l, i), b[l + 1], b[l]
            )
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
    """RT4's rules equal numpy's, ascending, for each type."""
    for quad in polradtran.QuadratureType:
        for nmu in (1, 2, 5, 8):
            q = polradtran.get_quadrature(nmu, quad)
            mu, w = own_quadrature(nmu, quad)
            np.testing.assert_allclose(
                np.asarray(q.mu), mu, rtol=0, atol=1e-14, err_msg=str(quad)
            )
            np.testing.assert_allclose(
                np.asarray(q.weights), w, rtol=0, atol=1e-14, err_msg=str(quad)
            )
        print(f"quadrature {quad.name:13s} == numpy for nmu = 1, 2, 5, 8")


def test_gas_fresnel():
    """Gas-only, Planck linear in tau, polarizing Fresnel surface; every level and stream.

    RT4 solves gas-only layers analytically, so only round-off remains.
    """
    gas = lambda l, i: (GAS[l], 0.0, GAS[l], 0.0)
    for quad in polradtran.QuadratureType:
        p = problem(quad=quad)
        r = rt4.solve(p)
        mu = np.asarray(r.mu)
        mu_q, w_q = own_quadrature(p.nmu, quad)
        np.testing.assert_allclose(
            mu, np.concatenate([mu_q, [0.45, 1.0]]), rtol=0, atol=1e-14
        )
        np.testing.assert_allclose(
            np.asarray(r.weights), np.concatenate([w_q, [0.0, 0.0]]), rtol=0, atol=1e-14
        )

        down, up = references(mu, gas, INDEX)
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

    def optics(l, i):
        if l == 1:
            return GAS[l] + kp[i], kappa[i], GAS[l] + a_i[i], a_q[i]
        return GAS[l], 0.0, GAS[l], 0.0

    down, up = references(mu, optics, INDEX)
    dev = max_rel_dev(r, down, up)
    tol = p.max_delta_tau / mu[0]
    print(
        f"layout, lower-triangular K(mu) + forward phase   max rel dev I {dev[0]:.2e}, Q {dev[1]:.2e} (tol {tol:.1e})"
    )
    assert np.all(dev < tol), dev


def rayleigh_m0(mo, mi):
    """Azimuthal mean of the Rayleigh phase matrix, 4 pi normalised, (I, Q),
    mo outgoing and mi incident mu (meridional basis, Q = I_v - I_h)."""
    a, b = mo**2, mi**2
    return np.array(
        [
            [3 / 8 * (3 - a - b + 3 * a * b), 3 / 8 * (1 - 3 * a) * (1 - b)],
            [3 / 8 * (1 - a) * (1 - 3 * b), 9 / 8 * (1 - a) * (1 - b)],
        ]
    )


def test_kirchhoff(reciprocal):
    """Isothermal sky, layers and surface: I = B, Q = 0 exactly on the streams
    when a(h, mu_j) = K [1, 0] - 2 pi sum_{i, h'} w_i Z(h, j <- h', i) [1, 0].

    The non-reciprocal variant multiplies the Q row of Z (Z_QI, Z_QQ) by
    (1 + 0.3 mu_out), which makes a_Q non-zero and catches a full transpose
    of Z that the reciprocal Rayleigh matrix is blind to.  Z_II stays
    reciprocal, so the medium conserves energy (rt4.solve checks it); a
    non-reciprocal Z_II could not obey Kirchhoff's law and energy
    conservation with one absorption.  Only round-off remains, which grows with the
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

    z = np.zeros((nmu, nmu, 2, 2))
    for i in range(nmu):
        for j in range(nmu):
            z[i, j] = sigma * rayleigh_m0(mu[i], mu[j]) / (4.0 * np.pi)
            if not reciprocal:
                z[i, j, 1, :] *= 1.0 + 0.3 * mu[i]
    phase = np.broadcast_to(z, (2, 2) + z.shape).copy()

    ext = np.zeros((2, nmu, 2, 2))
    ext[:, :, 0, 0] = ext[:, :, 1, 1] = sigma + kabs
    absorption = np.zeros((2, nmu, 2))
    absorption[:, :, 0] = sigma + kabs
    # both incoming hemispheres scatter into (h, j) identically
    absorption -= 2.0 * np.pi * 2.0 * np.einsum("i,jis->js", w, z[..., 0])[None, :, :]
    if reciprocal:
        assert np.max(np.abs(absorption[..., 1])) < 1e-12 * sigma
    p.optics = [rt4.LayerOptics(ext, absorption, phase)]

    r = rt4.solve(p)
    b = planck(FREQUENCY, t)
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


def expect_error(what, p, match):
    """solve() must raise before the Fortran code (a STOP there would end this process)."""
    try:
        rt4.solve(p)
    except RuntimeError as e:
        assert match in str(e), str(e)
        print(f"error path: {what:40s} raises")
        return
    raise AssertionError(f"{what} did not raise")


def test_errors():
    expect_error("nstokes = 3", problem(nstokes=3), "nstokes 1 ([I]) or 2 ([I, Q])")
    expect_error("nstokes * nmu_total > 64", problem(nmu=31), "<= 64")

    p = problem(layer_optics_index=[-1, 0, -1])
    nmu = p.nmu + len(p.extra_mu)
    ext = np.zeros((2, nmu, 2, 2))
    ext[DOWN, :, 0, 0] = ext[DOWN, :, 1, 1] = 1e-4
    ext[UP, :, 0, 0] = ext[UP, :, 1, 1] = 2e-4
    p.optics = [
        rt4.LayerOptics(ext, np.zeros((2, nmu, 2)), np.zeros((2, 2, nmu, nmu, 2, 2)))
    ]
    expect_error("up/down asymmetric extinction", p, "mirror symmetric")

    p.optics = [
        rt4.LayerOptics(
            np.zeros((2, nmu, 2, 2)),
            np.zeros((2, nmu, 2)),
            np.zeros((2, 2, nmu, nmu, 1, 1)),
        )
    ]
    expect_error("wrong phase shape", p, "must have extinction")


test_quadrature()
test_gas_fresnel()
test_layout()
test_kirchhoff(reciprocal=True)
test_kirchhoff(reciprocal=False)
test_errors()
