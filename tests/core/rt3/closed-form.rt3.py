"""RT3 through pyarts3.arts.rt3 against references from outside RT3.

* Evans' Mie benchmark: the scattering file and the expected output are read
  from his original script 3rdparty/polradtran/runmietest (the output of
  the original RT3 program, per micrometre).
* Single scattering by an optically thin Rayleigh layer, with the Stokes
  vector of the scattered field built from vector geometry in numpy.
* A gas-only atmosphere over a Fresnel surface, integrated analytically.

No RT3 output enters a reference.
"""

from pathlib import Path

import numpy as np
from pyarts3 import arts

rt3 = arts.rt3
assert rt3.available(), "this test is only collected with ENABLE_RT3=ON"

H = 6.62607015e-34
C = 299792458.0
K = 1.380649e-23

RUNMIETEST = Path(__file__).resolve().parents[3] / "3rdparty/polradtran/runmietest"


def heredoc(text, name):
    """The body of 'cat >name <<EOF ... EOF' in a csh script."""
    return text.split(f"cat >{name} <<EOF\n")[1].split("\nEOF")[0]


def planck(f, t):
    """B_nu(T) in W m-2 Hz-1 sr-1."""
    return 2.0 * H * f**3 / C**2 / np.expm1(H * f / (K * t))


def own_quadrature(nmu, quad):
    """The rules as documented, from numpy's Gauss-Legendre and Legendre roots."""
    if quad == rt3.QuadratureType.double_gauss:
        x, w = np.polynomial.legendre.leggauss(nmu)
        return 0.5 * (x + 1.0), 0.5 * w
    if quad == rt3.QuadratureType.gauss:
        x, w = np.polynomial.legendre.leggauss(2 * nmu)
        return x[nmu:], w[nmu:]
    n = 2 * nmu  # Lobatto: +-1 and the roots of P'_{n-1}
    p = np.polynomial.legendre.Legendre.basis(n - 1)
    x = np.sort(np.concatenate([[-1.0, 1.0], p.deriv().roots().real]))
    w = 2.0 / (n * (n - 1) * p(x) ** 2)
    return x[nmu:], w[nmu:]


def test_quadrature():
    for quad in rt3.QuadratureType:
        for nmu in (1, 2, 5, 8):
            q = rt3.get_quadrature(nmu, quad)
            mu, w = own_quadrature(nmu, quad)
            np.testing.assert_allclose(np.asarray(q.mu), mu, rtol=0, atol=1e-14)
            np.testing.assert_allclose(np.asarray(q.weights), w, rtol=0, atol=1e-14)
        print(f"quadrature {quad.name:13s} == numpy for nmu = 1, 2, 5, 8")
    assert rt3.max_legendre_degree(8, rt3.QuadratureType.gauss) == 29
    assert rt3.max_legendre_degree(8, rt3.QuadratureType.double_gauss) == 13
    assert rt3.max_legendre_degree(8, rt3.QuadratureType.lobatto) == 27


def print_half_unit(v):
    """Half a unit in the 6th significant digit of an E13.6 value; 0 for a printed 0."""
    v = np.abs(np.asarray(v))
    with np.errstate(divide="ignore"):
        e = np.floor(np.log10(np.where(v > 0, v, 1.0))) + 1.0
    return np.where(v > 0, 0.5 * 10.0 ** (e - 6), 0.0)


def real4_bound(c, phi_deg):
    """Bound on the REAL*4 error of rt3.f's OUTPUT_FILE for sum_m c_m t(m phi).

    PHI is REAL*4, M*PHI and the cosine are REAL*4, and so is the running
    sum: per term |c_m| (m |phi_f - phi| + 2^-24 m phi + 2^-23), plus one
    rounding of 2^-24 sum_m |c_m| per term.
    """
    eps = 2.0**-24
    phi = np.radians(phi_deg)
    phif = float(np.float32(np.pi) * np.float32(phi_deg) / np.float32(180.0))
    m = np.arange(c.shape[0])
    a = np.abs(c)
    return np.sum(a * (m * abs(phif - phi) + eps * m * phi + 2 * eps)) + len(
        m
    ) * eps * np.sum(a)


def test_mietest():
    """Evans and Stephens (1991) Mie case: tau = 1, omega = 0.99, Lambertian
    0.1, solar flux 0.2 pi on the horizontal at the zenith angle 78.46 deg,
    G quadrature with 8 nodes, aziorder 8, nstokes 4."""
    text = RUNMIETEST.read_text()
    sca = heredoc(text, "mietest.sca").splitlines()
    extinction, scattering = float(sca[0].split()[0]), float(sca[1].split()[0])
    nleg = int(sca[3].split()[0])
    legendre = np.array(
        [[float(x) for x in line.split()[1:]] for line in sca[4 : 4 + nleg + 1]]
    )
    table = np.array(
        [
            [float(x) for x in line.split()]
            for line in heredoc(text, "mietest.out.check").splitlines()
            if not line.startswith("C")
        ]
    )  # Z, PHI, MU, I, Q, U, V

    lam = 0.951
    f = C / (lam * 1e-6)
    per_um_to_per_hz = lam / f
    p = rt3.Problem(
        nstokes=4,
        nmu=8,
        quad=rt3.QuadratureType.gauss,
        aziorder=8,
        direct_flux=0.628318531 * per_um_to_per_hz,
        direct_mu=abs(
            np.cos(0.017453292 * 78.46304097)
        ),  # rt3.f USER_INPUT, truncated pi / 180
        thermal=False,
        frequency=f,
        height=[1.0, 0.0],
        temperature=[0.0, 0.0],
        gas_extinction=[0.0],
        scattering_sets=[rt3.ScatteringSet(extinction, scattering, legendre)],
        layer_scattering_index=[0],
        ground=rt3.LambertianSurface(0.1),
    )
    r = rt3.solve(p)
    mu = np.asarray(r.mu)
    up = np.asarray(r.up) / per_um_to_per_hz
    down = np.asarray(r.down) / per_um_to_per_hz
    phis = np.unique(table[:, 1])
    rad_up = np.asarray(rt3.azimuth_radiance(r.up, np.radians(phis))) / per_um_to_per_hz
    rad_down = (
        np.asarray(rt3.azimuth_radiance(r.down, np.radians(phis))) / per_um_to_per_hz
    )

    worst = 0.0
    for z, phi, m, *ref in table:
        level = 0 if z > 0.5 else 1
        if abs(m) > 1.5:
            flux = (
                np.asarray(r.up_flux if m < 0 else r.down_flux)[level]
                / per_um_to_per_hz
            )
            got, bound = flux, 2.0**-24 * np.abs(flux)
        else:
            i = int(np.argmin(np.abs(mu - abs(m))))
            assert abs(mu[i] - abs(m)) <= 5.000001e-6
            k = int(np.argmin(np.abs(phis - phi)))
            got = (rad_up if m < 0 else rad_down)[level, k, i]
            c = (up if m < 0 else down)[level, :, i]
            bound = np.array([real4_bound(c[:, s], phi) for s in range(4)])
        tol = print_half_unit(ref) + bound
        dev = np.abs(got - ref)
        # A printed 0 with a zero bound (e.g. U at phi = 0) must be exactly 0
        ratio = np.divide(dev, tol, out=np.where(dev > 0, np.inf, 0.0), where=tol > 0)
        worst = max(worst, float(np.max(ratio)))
    print(
        f"runmietest: {len(table)} rows, max |RT3 - table| / (print half-unit + REAL*4 bound) = {worst:.2f}"
    )
    assert worst <= 1.0, worst


def rayleigh_column(k_out, k_in):
    """Dipole scattering of unpolarized light from k_in into k_out.

    The Stokes vector [I, Q, U, V] of the scattered field in the meridional
    basis of k_out (z up, h = k x z / |k x z|, v = h x k, Q = I_v - I_h,
    U = 2 E_v E_h), scaled so that I = 3/4 (1 + cos^2 Theta).
    """

    def basis(k):
        h = np.cross(k, [0.0, 0.0, 1.0])
        h /= np.linalg.norm(h)
        return np.cross(h, k), h

    vi, hi = basis(k_in)
    vo, ho = basis(k_out)
    z = np.zeros(4)
    for e in (vi, hi):
        ev, eh = vo @ e, ho @ e
        z += 0.75 * np.array([ev * ev + eh * eh, ev * ev - eh * eh, 2 * ev * eh, 0.0])
    return z


def direction(mu_z, phi):
    s = np.sqrt(1.0 - mu_z**2)
    return np.array([s * np.cos(phi), s * np.sin(phi), mu_z])


def test_single_scattering():
    """A conservative Rayleigh layer of tau = 1e-6 over a black surface, lit
    by a beam of flux F on the horizontal at mu0 = 0.6 propagating toward
    phi = 0: exact single scattering, up at the top and down at the bottom.
    The neglected multiple scattering is of relative order tau."""
    tau, mu0, flux = 1e-6, 0.6, 2.5
    rayleigh = np.array(
        [
            [1.0, -0.5, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.5, 0.0, 0.0, 1.5],
            [0.5, 0.5, 0.0, 0.0, 0.5, 0.0],
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
        ground=rt3.LambertianSurface(0.0),
    )
    r = rt3.solve(p)
    phi = np.radians([0.0, 40.0, 90.0, 160.0, 200.0, 300.0])
    up = np.asarray(rt3.azimuth_radiance(r.up, phi))[0]
    down = np.asarray(rt3.azimuth_radiance(r.down, phi))[1]

    k0 = direction(-mu0, 0.0)
    ref_up, ref_down = np.zeros_like(up), np.zeros_like(down)
    for k, ph in enumerate(phi):
        for i, mu in enumerate(np.asarray(r.mu)):
            ref_up[k, i] = (
                flux
                / (4 * np.pi)
                * rayleigh_column(direction(mu, ph), k0)
                * -np.expm1(-tau * (1 / mu0 + 1 / mu))
                / (mu0 + mu)
            )
            ref_down[k, i] = (
                flux
                / (4 * np.pi)
                * rayleigh_column(direction(-mu, ph), k0)
                * np.exp(-tau / mu0)
                * np.expm1(tau / mu0 - tau / mu)
                / (mu - mu0)
            )
    scale = max(np.abs(ref_up).max(), np.abs(ref_down).max())
    assert np.abs(ref_up[..., 2]).max() > 0.1 * scale
    dev = max(np.abs(up - ref_up).max(), np.abs(down - ref_down).max()) / scale
    print(
        f"thin Rayleigh layer vs vector geometry: max |dev| / max |I| = {dev:.2e} (tau = {tau:.0e})"
    )
    assert dev < 10 * tau, dev


def test_gas_fresnel():
    """Gas-only, Planck linear in tau, polarizing Fresnel surface, nstokes 4;
    every level and stream, m = 0, and every m > 0 mode exactly 0."""
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
        ground=rt3.FresnelSurface(n),
    )
    r = rt3.solve(p)
    mu = np.asarray(r.mu)
    b = planck(f, temperature)
    bs = planck(f, tsurf)
    dtau = gas * np.abs(np.diff(height))

    def path(i0, x, b0, b1):
        ex = np.exp(-x)
        return i0 * ex + b1 * -np.expm1(-x) - (b1 - b0) * (1 - (1 + x) * ex) / x

    down = np.zeros((4, len(mu), 4))
    up = np.zeros((4, len(mu), 4))
    down[0, :, 0] = planck(f, sky)
    for l in range(3):
        down[l + 1, :, 0] = path(down[l, :, 0], dtau[l] / mu, b[l], b[l + 1])
    cos_t = np.sqrt(1.0 - (1.0 - mu**2) / n**2 + 0j)
    rv2 = np.abs((n * mu - cos_t) / (n * mu + cos_t)) ** 2
    rh2 = np.abs((mu - n * cos_t) / (mu + n * cos_t)) ** 2
    r1, r2 = 0.5 * (rv2 + rh2), 0.5 * (rv2 - rh2)
    up[3, :, 0] = (1 - r1) * bs + r1 * down[3, :, 0]
    up[3, :, 1] = r2 * (down[3, :, 0] - bs)
    for l in reversed(range(3)):
        up[l, :, 0] = path(up[l + 1, :, 0], dtau[l] / mu, b[l + 1], b[l])
        up[l, :, 1] = up[l + 1, :, 1] * np.exp(-dtau[l] / mu)

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


def expect_error(what, p, match):
    """solve() must raise before the Fortran code (a STOP there would end this process)."""
    try:
        rt3.solve(p)
    except RuntimeError as e:
        assert match in str(e), str(e)
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
        ground=rt3.LambertianSurface(0.1),
    )
    mie = np.array(
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
        "nstokes 1 to 4",
    )
    expect_error(
        "beam over a Fresnel surface",
        rt3.Problem(
            **{**base, "ground": rt3.FresnelSurface(1.5)}, scattering_sets=[ok]
        ),
        "Lambertian",
    )
    expect_error(
        "Mie series of degree 11, nmu = 2",
        rt3.Problem(**{**base, "nmu": 2}, scattering_sets=[ok]),
        "silently truncate",
    )
    expect_error(
        "unnormalised phase function",
        rt3.Problem(**base, scattering_sets=[rt3.ScatteringSet(1.0, 0.9, 0.5 * mie)]),
        "normalised",
    )
    expect_error(
        "extra_mu with lobatto",
        rt3.Problem(
            **base,
            quad=rt3.QuadratureType.lobatto,
            extra_mu=[0.5],
            scattering_sets=[ok],
        ),
        "only to the gauss",
    )


test_quadrature()
test_mietest()
test_single_scattering()
test_gas_fresnel()
test_errors()
