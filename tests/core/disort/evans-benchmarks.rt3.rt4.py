"""Evans' PolRadTran benchmarks with ARTS's RT3, RT4 and VDISORT, plotted.

Evans' four scripts in 3rdparty/polradtran (runmietest, runtesta, runtestr,
runtestc) define the problems and hold his expected outputs, the tables.  The
problems are read from the scripts and solved by

- ARTS's RT3 and RT4 (pyarts3.arts.rt3, rt4) at Evans' settings;
- VDISORT (pyarts3.arts.cppvdisort) on its own double-Gauss streams;
- RT3 with Gauss quadrature of increasing order, and with VDISORT's
  double-Gauss streams.

What the tables are.  Evans' tables are solutions with Gauss quadrature of
few streams (runmietest 8, runtesta 4, runtestr 8 per hemisphere), not the
exact radiance.  Gauss quadrature on [-1, 1] handles the jump of the radiance
at the horizon poorly, so it converges slowly (about like 1 / nmu).  So:

- ARTS's RT3 and RT4 at Evans' settings must reproduce his tables (to print
  precision), as his own programs did;
- RT3 and VDISORT on the same double-Gauss streams solve the same discrete
  problem and must agree to RT3's doubling error (about 1e-7);
- RT3 with Gauss quadrature must approach VDISORT as its number of streams
  grows.  At Evans' number it is his table, so VDISORT differs from the
  tables by their quadrature error, not by its own.

runtestc (horizontally oriented ice columns) has direction-dependent
extinction, which only RT4 represents; it is solved by RT4 alone.

VDISORT is evaluated on its streams (64 per hemisphere); its values at
Evans' angles are cubic Lagrange interpolations in mu, checked against 128
streams.  Without ARTS_HEADLESS the solutions and differences are plotted.
"""

import os
from pathlib import Path

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts
POLRADTRAN = Path(__file__).resolve().parents[3] / "3rdparty" / "polradtran"
HEADLESS = "ARTS_HEADLESS" in os.environ

# Exact SI radiation constants in Evans' units: 2 h c^2 [W m-2 sr-1 um^4], h c / k [um K]
H, C, K = 6.62607015e-34, 299792458.0, 1.380649e-23
C1, C2 = 2 * H * C * C * 1e24, H * C / K * 1e6


###############################################################################
# Evans' scripts
###############################################################################


def read_script(name):
    """The here-documents of a script: written files and program runs (program, answers)."""
    files, runs = {}, []
    lines = (POLRADTRAN / name).read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if " <<EOF" not in line:
            continue
        head = line.split(" <<EOF")[0]
        body = []
        while lines[i] != "EOF":
            body.append(lines[i])
            i += 1
        i += 1
        if head.startswith("cat >"):
            files[head[5:]] = "\n".join(body) + "\n"
        else:
            runs.append((head, body))
    check = next(f for f in files if f.endswith(".check"))
    return {"name": name, "files": files, "runs": runs, "check": check}


def read_output(text):
    """Rows (z, phi, mu, values) of an rt3.f or rt4.f output; phi is 0 for rt4.f."""
    rows, has_phi = [], False
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("C"):
            if " MU " in line:
                has_phi = " PHI " in line
            continue
        x = [float(t) for t in line.split()]
        if has_phi:
            rows.append((x[0], x[1], x[2], np.array(x[3:])))
        else:
            rows.append((x[0], 0.0, x[1], np.array(x[2:])))
    return rows


def read_layers(text):
    levels = []
    for line in text.splitlines():
        parts = line.split("'")
        numbers = parts[0].split()
        if len(numbers) < 3:
            continue
        levels.append((float(numbers[0]), float(numbers[1]), float(numbers[2]), parts[1].strip() if len(parts) > 1 else ""))
    return levels


def read_legendre(text):
    """An RT3 scattering file: extinction, scattering and the [L + 1, 6] series (F11, F12, F33, F34, F22, F44)."""
    lines = [l for l in text.splitlines() if not l.startswith("C")]
    ext, sca, degree = float(lines[0].split()[0]), float(lines[1].split()[0]), int(lines[3].split()[0])
    coef = np.array([[float(t) for t in lines[4 + l].split()[1:7]] for l in range(degree + 1)])
    return ext, sca, coef


def read_rt4_scattering(text, ns):
    """An RT4 scattering file (GET_SCAT_FILE): extinction [2, n, ns, ns], absorption [2, n, ns],
    phase [2 out, 2 in, n out, n in, ns, ns]; hemisphere 0 is mu > 0 in the file, propagating down."""
    lines = [l for l in text.splitlines() if l.strip() and l.strip()[0] != "C"]
    n, naz = (int(t) for t in lines[0].split()[:2])
    assert naz == 0
    x = np.array([float(t) for l in lines[1:] for t in l.split()])
    pos = 0

    def take(k):
        nonlocal pos
        pos += k
        return x[pos - k : pos]

    phase = np.zeros((2, 2, n, n, ns, ns))
    for l1 in range(2):
        for j1 in range(n):
            for l2 in range(2):
                for j2 in range(n):
                    take(3)
                    phase[l2, l1, j2, j1] = take(16).reshape(4, 4)[:ns, :ns]
    extinction = np.zeros((2, n, ns, ns))
    for l in range(2):
        for j in range(n):
            take(1)
            extinction[l, j] = take(16).reshape(4, 4)[:ns, :ns]
    absorption = np.zeros((2, n, ns))
    for l in range(2):
        for j in range(n):
            take(1)
            absorption[l, j] = take(4)[:ns]
    return extinction, absorption, phase


def answers(run):
    return [a.split()[0] if a.split() else "" for a in run[1]]


def planck5(lam, t):
    """Evans' 5-digit Planck function [W m-2 sr-1 um-1]."""
    t = np.asarray(t, dtype=float)
    with np.errstate(divide="ignore", over="ignore"):
        b = 1.1911e8 / lam**5 / np.expm1(1.4388e4 / (lam * t))
    return np.where(t > 0, b, 0.0)


def exact_t(lam, t):
    """The temperature at which the exact Planck function is Evans' 5-digit one at t."""
    if t == 0:
        return 0.0
    return C2 / (lam * np.log1p(C1 / (lam**5 * planck5(lam, t))))


def brightness_vh(i, q, lam, flux=False):
    """Evans' CONVERT_OUTPUT (units T, polarization VH): [I, Q] per um to the EBB temperatures of V and H."""
    out = []
    for rad in (i + q, i - q):
        rad = rad / np.pi if flux else rad
        with np.errstate(divide="ignore"):  # rad = 0 gives T = 0
            out.append(np.sign(rad) * 1.4388e4 / (lam * np.log1p(1.1911e8 / (np.abs(rad) * lam**5))))
    return np.array(out)


def read_case(name):
    """A script's problem in Evans' units (heights in km, per um), with his table."""
    s = read_script(name)
    solver, a = s["runs"][-1][0], answers(s["runs"][-1])
    case = {"name": name, "script": s, "program": solver, "table": read_output(s["files"][s["check"]])}
    if solver == "rt3":
        case.update(nstokes=int(a[0]), nmu=int(a[1]), quad=a[2][0], aziorder=int(a[3]), layer_file=a[4])
        src = int(a[6])
        k = 7
        case["flux"], case["mu0"] = 0.0, 1.0
        if src in (1, 3):
            case["flux"], case["mu0"] = float(a[7]), abs(np.cos(0.017453292 * float(a[8])))
            k = 9
        case["thermal"] = src >= 2
        case["ground_t"], gtype = float(a[k]), a[k + 1][0]
        rest = s["runs"][-1][1][k + 2 :]
        case["brightness"] = False
    else:
        case.update(nstokes=int(a[0]), nmu=int(a[1]), quad=a[2][0], aziorder=0, layer_file=a[3])
        case["flux"], case["mu0"], case["thermal"] = 0.0, 1.0, True
        case["ground_t"], gtype = float(a[4]), a[5][0]
        rest = s["runs"][-1][1][6:]
        case["brightness"] = True
    if gtype == "F":
        case["fresnel"] = complex(rest[0].split()[0].replace("(", "").replace(")", "").replace(",", "+").replace("+-", "-") + "j")
        case["albedo"] = 0.0
    else:
        case["fresnel"], case["albedo"] = None, float(rest[0].split()[0])
    case["sky_t"], case["wavelength"] = float(rest[1].split()[0]), float(rest[2].split()[0])

    # The layers; RT4 files made by scatcnv come from its Legendre input
    legendre_of = {}
    for prog, ans in s["runs"]:
        if prog == "scatcnv":
            legendre_of[ans[1].split()[0]] = ans[0].split()[0]
    levels = read_layers(s["files"][case["layer_file"]])
    case["height"] = np.array([l[0] for l in levels])
    case["temperature"] = np.array([l[1] for l in levels])
    case["gas"] = np.array([l[2] for l in levels[:-1]])
    case["files"] = [l[3] for l in levels[:-1]]
    case["legendre"] = {}
    for f in set(case["files"]) - {""}:
        source = legendre_of.get(f, f)
        if source in s["files"]:
            case["legendre"][f] = read_legendre(s["files"][source])
    case["randomly_oriented"] = all(f == "" or f in case["legendre"] for f in case["files"])
    return case


###############################################################################
# Phase matrices from a Legendre series, by vector geometry
###############################################################################


def _direction(mu, phi):
    s = np.sqrt(np.maximum(0.0, 1.0 - mu * mu))
    return np.stack(np.broadcast_arrays(s * np.cos(phi), s * np.sin(phi), mu), -1)


def _rotation(normal, mu, phi):
    """The Stokes rotation from the meridional basis of (mu, phi) to the scattering plane of normal."""
    mu, phi = np.broadcast_arrays(mu, phi)
    s = np.sqrt(np.maximum(0.0, 1.0 - mu * mu))
    ev = np.stack([mu * np.cos(phi), mu * np.sin(phi), -s], -1)
    eh = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], -1)
    par = np.cross(normal, _direction(mu, phi))
    c, d = np.sum(par * ev, -1), np.sum(par * eh, -1)
    L = np.zeros(c.shape + (4, 4))
    L[..., 0, 0] = L[..., 3, 3] = 1.0
    L[..., 1, 1] = L[..., 2, 2] = c * c - d * d
    L[..., 1, 2] = 2 * c * d
    L[..., 2, 1] = -2 * c * d
    return L


def lab_frame(coef, mu_in, phi_in, mu_out):
    """Z(mu_out, 0 <- mu_in, phi_in) of the series coef (F11, F12, F33, F34, F22, F44), normalised to 1 over 4 pi."""
    mu_in, phi_in, mu_out = np.broadcast_arrays(mu_in, phi_in, mu_out)
    ki, ko = _direction(mu_in, phi_in), _direction(mu_out, 0.0 * mu_out)
    normal = np.cross(ki, ko)
    norm = np.linalg.norm(normal, axis=-1)
    fallback = np.stack([-np.sin(phi_in), np.cos(phi_in), np.zeros_like(phi_in)], -1)
    normal = np.where((norm < 1e-12)[..., None], fallback, normal / np.maximum(norm, 1e-300)[..., None])
    x = np.clip(np.sum(ki * ko, -1), -1.0, 1.0)
    f = [np.polynomial.legendre.legval(x, coef[:, k]) for k in range(6)]
    F = np.zeros(x.shape + (4, 4))
    F[..., 0, 0], F[..., 0, 1], F[..., 1, 0], F[..., 1, 1] = f[0], f[1], f[1], f[4]
    F[..., 2, 2], F[..., 2, 3], F[..., 3, 2], F[..., 3, 3] = f[2], f[3], -f[3], f[5]
    Lo = _rotation(normal, mu_out, 0.0 * mu_out)
    Li = _rotation(normal, mu_in, phi_in)
    return np.swapaxes(Lo, -1, -2) @ F @ Li


def fourier_modes(coef, mu_out, mu_in, nfourier, nstokes, nphi=32):
    """Ordinary C^m, S^m [m, out, in, 4, 4] = (1 / 2 pi) int Z(mu_out, 0; mu_in, phi) {cos, sin}(m phi) dphi.

    The trapezoid rule on nphi samples is exact for a series of degree L and m < nfourier if L + m < nphi
    (here at most 13 + 8)."""
    degree = coef.shape[0] - 1
    assert degree + nfourier - 1 < nphi, "too few azimuth samples for an exact transform"
    phi = 2 * np.pi * np.arange(nphi) / nphi
    Z = lab_frame(coef, mu_in[None, :, None], phi[None, None, :], mu_out[:, None, None])  # [out, in, phi, 4, 4]
    Z[..., nstokes:, :] = 0.0
    Z[..., :, nstokes:] = 0.0
    m = np.arange(nfourier)[:, None]
    cos, sin = np.cos(m * phi) / nphi, np.sin(m * phi) / nphi
    return np.einsum("mk,oikab->moiab", cos, Z), np.einsum("mk,oikab->moiab", sin, Z)


###############################################################################
# The solvers
###############################################################################


def layers(case):
    """Per layer: total extinction per km, single-scattering albedo, and the series (or None)."""
    out = []
    for gas, f in zip(case["gas"], case["files"]):
        ext, sca, coef = case["legendre"][f] if f else (0.0, 0.0, None)
        out.append((gas + ext, sca / (gas + ext), coef))
    return out


def vdisort(case, nmu):
    """VDISORT on nmu double-Gauss streams per hemisphere, in Evans' units (5-digit Planck function, per um).
    Returns the signed streams (> 0 upward) and the radiance [level, psi = 0, 90, 180 deg, stream, 4]."""
    nodes, _ = np.polynomial.legendre.leggauss(nmu)
    pos = 0.5 * (nodes + 1.0)
    mu = np.concatenate((pos, -pos))
    nq, nf, lam = 2 * nmu, case["aziorder"] + 1, case["wavelength"]
    lay = layers(case)
    nl = len(lay)
    dz = -np.diff(case["height"])
    tau = np.cumsum([k * z for (k, _, _), z in zip(lay, dz)])

    C = np.zeros((nf, nl, nq, nq, 4, 4))
    S = np.zeros_like(C)
    Cb = np.zeros((nf, nl, nq, 4, 4))
    Sb = np.zeros_like(Cb)
    for l, (_, _, coef) in enumerate(lay):
        if coef is None:
            continue
        C[:, l], S[:, l] = fourier_modes(coef, mu, mu, nf, case["nstokes"])
        if case["flux"] > 0:
            c, s = fourier_modes(coef, mu, np.array([-case["mu0"]]), nf, case["nstokes"])
            Cb[:, l], Sb[:, l] = c[:, :, 0], s[:, :, 0]

    B = planck5(lam, case["temperature"]) if case["thermal"] else np.zeros(nl + 1)
    source = np.zeros((nl, 2, 4))
    top = 0.0
    for l in range(nl):
        slope = (B[l + 1] - B[l]) / (tau[l] - top)
        source[l, 0, 0], source[l, 1, 0] = B[l] - slope * top, slope
        top = tau[l]
    b_top = np.zeros((2, nf, nmu, 4))
    b_top[0, 0, :, 0] = planck5(lam, case["sky_t"])
    b_bottom = np.zeros((2, nf, nmu, 4))
    Bs = planck5(lam, case["ground_t"])
    if case["fresnel"] is not None:
        for i in range(nmu):
            R = np.asarray(A.vdisort.fresnel_reflection(pos[i], case["fresnel"]))
            b_bottom[0, 0, i] = Bs * (np.array([1.0, 0, 0, 0]) - R[:, 0])
        brdf = A.vdisort.fresnel_fourier_modes(case["fresnel"], nf)
    else:
        b_bottom[0, 0, :, 0] = (1 - case["albedo"]) * Bs if case["thermal"] else 0.0
        brdf = A.vdisort.lambertian_fourier_modes(case["albedo"], nf)

    kwargs = {}
    if case["flux"] > 0:
        kwargs["beam_phase_matrix"] = np.asarray(A.vdisort.combine_beam_phase_matrices(Cb, Sb))
    model = A.cppvdisort(
        tau_arr=tau,
        omega_arr=np.array([w for _, w, _ in lay]),
        NQuad=nq,
        NFourier=nf,
        phase_matrix=np.asarray(A.vdisort.combine_phase_matrices(C, S)),
        mu0=case["mu0"],
        beam_stokes=np.array([case["flux"] / case["mu0"] if case["flux"] > 0 else 0.0, 0, 0, 0]),
        phi0=0.0,
        b_pos=b_bottom,
        b_neg=b_top,
        BDRF_Fourier_modes=brdf,
        s_poly_coeffs=source,
        **kwargs,
    )
    assert np.allclose(np.asarray(model.mu), mu, rtol=0, atol=1e-14), "VDISORT's streams must be double-Gauss"
    u = np.asarray(model.u(tau=np.concatenate(([0.0], tau)), phi=np.array([0.0, np.pi / 2, np.pi])))
    return mu, u


def vdisort_at(solution, mus, level, k):
    """VDISORT at signed Evans cosines (< 0 upwelling, Evans' sign), cubic Lagrange in mu per hemisphere."""
    mu, u = solution
    n = len(mu) // 2
    pos = mu[:n]
    out = []
    for m in np.atleast_1d(mus):
        values = u[level, k, :n] if m < 0 else u[level, k, n:]  # Evans' mu < 0 is upwelling, VDISORT's mu > 0
        x = abs(m)
        j = np.clip(np.searchsorted(pos, x) - 2, 0, n - 4)
        xs = pos[j : j + 4]
        w = [np.prod([(x - xs[b]) / (xs[a] - xs[b]) for b in range(4) if b != a]) for a in range(4)]
        out.append(sum(w[a] * values[j + a] for a in range(4)))
    return np.array(out)


def rt3(case, nmu, quad, max_delta_tau=1e-7):
    """ARTS's RT3 (exact Planck function: the temperatures are those of Evans' 5-digit values, results per um).
    Returns the streams and a function (level, k, mu) -> [values], with Evans' signed mu."""
    lam = case["wavelength"]
    f = C / (lam * 1e-6)
    per_um = lam / f
    sets, index = [], []
    for name in case["files"]:
        if not name:
            index.append(-1)
            continue
        ext, sca, coef = case["legendre"][name]
        index.append(len(sets))
        sets.append(A.rt3.ScatteringSet(ext, sca, coef))
    ground = A.polradtran.FresnelSurface(case["fresnel"]) if case["fresnel"] is not None else A.polradtran.LambertianSurface(case["albedo"])
    p = A.rt3.Problem(
        nstokes=case["nstokes"],
        nmu=nmu,
        quad=getattr(A.polradtran.QuadratureType, quad),
        aziorder=case["aziorder"],
        max_delta_tau=max_delta_tau,
        direct_flux=case["flux"] * per_um,
        direct_mu=case["mu0"],
        thermal=case["thermal"],
        frequency=f,
        height=case["height"],
        temperature=np.array([exact_t(lam, t) for t in case["temperature"]]),
        gas_extinction=case["gas"],
        scattering_sets=sets,
        layer_scattering_index=index,
        sky_temperature=exact_t(lam, case["sky_t"]),
        surface_temperature=exact_t(lam, case["ground_t"]),
        ground=ground,
    )
    r = A.rt3.solve(p)
    psi = np.array([0.0, np.pi / 2, np.pi])
    up = np.asarray(A.rt3.azimuth_radiance(r.up, psi)) / per_um
    down = np.asarray(A.rt3.azimuth_radiance(r.down, psi)) / per_um
    mus = np.asarray(r.mu)
    fluxes = (np.asarray(r.up_flux) / per_um, np.asarray(r.down_flux) / per_um)

    def value(level, k, m):
        if abs(m) == 2:
            return fluxes[0 if m < 0 else 1][level]
        j = int(np.argmin(np.abs(mus - abs(m))))
        return (up if m < 0 else down)[level, k, j]

    return mus, value


def rt4_optics_from_legendre(coef, ext, sca, mu, ns):
    """RT4 optics on the streams mu of a randomly oriented series: the azimuthal mean of k_s Z / (4 pi)."""
    signed = np.concatenate((-mu, mu))  # RT4 hemisphere 0 propagates down
    c0, _ = fourier_modes(coef, signed, signed, 1, ns)
    n = len(mu)
    phase = np.zeros((2, 2, n, n, ns, ns))
    for ho in range(2):
        for hi in range(2):
            phase[ho, hi] = sca / (4 * np.pi) * c0[0, ho * n : (ho + 1) * n, hi * n : (hi + 1) * n, :ns, :ns]
    extinction = np.zeros((2, n, ns, ns))
    for s in range(ns):
        extinction[:, :, s, s] = ext
    absorption = np.zeros((2, n, ns))
    absorption[:, :, 0] = ext - sca
    return extinction, absorption, phase


def rt4(case):
    """ARTS's RT4 at Evans' settings: (level, mu) -> [I, Q] per um, with Evans' optics."""
    lam, ns = case["wavelength"], case["nstokes"]
    f = C / (lam * 1e-6)
    per_um = lam / f
    quad = {"G": A.polradtran.QuadratureType.gauss, "L": A.polradtran.QuadratureType.lobatto, "D": A.polradtran.QuadratureType.double_gauss}[case["quad"]]
    mu = np.asarray(A.polradtran.get_quadrature(case["nmu"], quad).mu)
    optics, index, made = [], [], {}
    for name in case["files"]:
        if not name:
            index.append(-1)
            continue
        if name not in made:
            if name in case["legendre"]:
                ext, sca, coef = case["legendre"][name]
                o = rt4_optics_from_legendre(coef, ext, sca, mu, ns)
            else:
                o = read_rt4_scattering((POLRADTRAN / name).read_text(), ns)
            made[name] = len(optics)
            optics.append(A.rt4.LayerOptics(*o))
        index.append(made[name])
    ground = A.polradtran.FresnelSurface(case["fresnel"]) if case["fresnel"] is not None else A.polradtran.LambertianSurface(case["albedo"])
    p = A.rt4.Problem(
        nstokes=ns,
        nmu=case["nmu"],
        quad=quad,
        max_delta_tau=1e-6,
        frequency=f,
        height=case["height"],
        temperature=np.array([exact_t(lam, t) for t in case["temperature"]]),
        gas_extinction=case["gas"],
        optics=optics,
        layer_optics_index=index,
        sky_temperature=exact_t(lam, case["sky_t"]),
        surface_temperature=exact_t(lam, case["ground_t"]),
        ground=ground,
    )
    r = A.rt4.solve(p)
    up, down = np.asarray(r.up) / per_um, np.asarray(r.down) / per_um
    w, mus = np.asarray(r.weights), np.asarray(r.mu)

    def value(level, k, m):
        if abs(m) == 2:
            rad = up if m < 0 else down
            return 2 * np.pi * np.sum(w[:, None] * mus[:, None] * rad[level], axis=0)
        j = int(np.argmin(np.abs(mus - abs(m))))
        return (up if m < 0 else down)[level, j]

    return mus, value


###############################################################################
# Comparison and plots
###############################################################################


def level_of(case, z):
    return int(np.argmin(np.abs(case["height"] - z)))


def as_table(case, row, values):
    """A solution's Stokes values in the table's quantity: radiance, or V, H brightness temperature."""
    if case["brightness"]:
        return brightness_vh(values[0], values[1], case["wavelength"], abs(row[2]) == 2)
    return np.asarray(values)[: len(row[3])]


def table_deviation(case, solution, rows=None):
    """max |solution - table| over the rows: radiances relative to max I and fluxes to max F, or in K."""
    rows = case["table"] if rows is None else rows
    scale = {False: 1.0, True: 1.0}
    if not case["brightness"]:
        for flux in (False, True):
            scale[flux] = max([abs(r[3][0]) for r in case["table"] if (abs(r[2]) == 2) == flux] or [1.0])
    d = 0.0
    for r in rows:
        x = as_table(case, r, solution(level_of(case, r[0]), int(round(r[1] / 90)), r[2]))
        d = max(d, np.max(np.abs(x - r[3][: len(x)])) / scale[abs(r[2]) == 2])
    return d


def radiance_rows(case):
    return [r for r in case["table"] if abs(r[2]) != 2]


results = {}
for name in ("runmietest", "runtesta", "runtestr", "runtestc"):
    case = read_case(name)
    unit = "K" if case["brightness"] else "of max I"
    print(f"{name}: {case['program']}, nstokes {case['nstokes']}, {case['quad']} nmu {case['nmu']}, "
          f"{len(case['height']) - 1} layer(s), {'Fresnel' if case['fresnel'] is not None else 'Lambertian'} ground")
    res = {"case": case}

    if case["program"] == "rt3":
        evans = rt3(case, case["nmu"], {"G": "gauss", "D": "double_gauss", "L": "lobatto"}[case["quad"]], 1e-6)
    else:
        evans = rt4(case)
    d = table_deviation(case, evans[1])
    print(f"    ARTS {case['program'].upper()} at Evans' settings vs his table: {d:.2e} {unit}")
    assert d <= (0.01 + 1e-9 if case["brightness"] else 2e-6)
    res["arts"] = evans

    if case["randomly_oriented"]:
        if case["program"] == "rt4":
            rt3_evans = rt3(case, case["nmu"], "gauss", 1e-6)
            d = table_deviation(case, rt3_evans[1])
            print(f"    ARTS RT3 at Evans' settings vs his RT4 table: {d:.2e} {unit}")
            assert d <= 0.01 + 1e-9
        v64, v32 = vdisort(case, 64), vdisort(case, 32)
        res["vdisort"] = v64
        vfun = lambda sol: (lambda level, k, m: vdisort_at(sol, [m], level, k)[0])  # noqa: E731
        rows = radiance_rows(case)
        # VDISORT at Evans' angles: 32 against 64 streams bounds the error of the 64-stream interpolation
        conv = 0.0
        for r in rows:
            a = as_table(case, r, vfun(v64)(level_of(case, r[0]), int(round(r[1] / 90)), r[2]))
            b = as_table(case, r, vfun(v32)(level_of(case, r[0]), int(round(r[1] / 90)), r[2]))
            conv = max(conv, np.max(np.abs(a - b)))
        scale = 1.0 if case["brightness"] else max(abs(r[3][0]) for r in rows)
        conv /= scale
        print(f"    VDISORT at Evans' angles, 32 vs 64 streams: {conv:.2e} {unit}")
        d_interp = table_deviation(case, vfun(v64), rows)
        print(f"    VDISORT (64 streams) at Evans' angles vs his table: {d_interp:.2e} {unit}")
        res["vdisort_vs_table"] = d_interp

        # RT3 and VDISORT on the same double-Gauss streams
        same = 0.0
        for nmu in (8, 16):
            mus, f3 = rt3(case, nmu, "double_gauss")
            mu_v, u_v = vdisort(case, nmu)
            scale = max(abs(u_v[..., 0]).max(), 1e-300)
            for level in range(len(case["height"])):
                for k in range(3 if case["aziorder"] > 0 else 1):
                    for j in range(nmu):
                        for sign, idx in ((-1, j), (1, nmu + j)):  # Evans' -mu is VDISORT's up stream j
                            a = f3(level, k, sign * mus[j])[: case["nstokes"]]
                            b = u_v[level, k, idx, : case["nstokes"]]
                            same = max(same, np.max(np.abs(a - b)) / scale)
        print(f"    ARTS RT3 and VDISORT on the same double-Gauss streams (8, 16): {same:.2e} of max I")
        assert same <= 2e-6
        res["same_streams"] = same

        # RT3 with Gauss quadrature approaches VDISORT
        sweep = []
        for nmu in (4, 6, 8, 12, 16):
            if case["nstokes"] * nmu > 64 or 4 * nmu - 3 < max(c.shape[0] - 1 for _, _, c in case["legendre"].values()):
                continue
            mus, f3 = rt3(case, nmu, "gauss")
            gauss_rows = [(z, 90.0 * k, s * m, np.zeros(case["nstokes"] if not case["brightness"] else 2))
                          for z in case["height"] for k in range(3 if case["aziorder"] > 0 else 1) for s in (-1, 1) for m in mus]
            d = 0.0
            sc = 1.0 if case["brightness"] else max(abs(r[3][0]) for r in rows)  # max I of the table, as above
            for r in gauss_rows:
                level, k = level_of(case, r[0]), int(round(r[1] / 90))
                a = as_table(case, r, f3(level, k, r[2])[: case["nstokes"]])
                b = as_table(case, r, vdisort_at(v64, [r[2]], level, k)[0][: case["nstokes"]])
                d = max(d, np.max(np.abs(a - b)) / sc)
            sweep.append((nmu, d))
            print(f"    ARTS RT3 Gauss nmu {nmu:2d} vs VDISORT at RT3's nodes: {d:.2e} {unit}")
        assert all(b[1] < a[1] for a, b in zip(sweep, sweep[1:])), "RT3 with Gauss quadrature must approach VDISORT"
        assert conv <= 0.1 * sweep[-1][1], "VDISORT at Evans' angles must be 10 times more accurate than the smallest difference shown"
        res["sweep"] = sweep
    else:
        print("    oriented particles: VDISORT and RT3 do not apply (direction-dependent extinction)")
    results[name] = res


if not HEADLESS:
    import matplotlib.pyplot as plt

    for name, res in results.items():
        case = res["case"]
        labels = ["V [K]", "H [K]"] if case["brightness"] else ["I", "Q", "U", "V"][: case["nstokes"]]
        panels = sorted({(r[0], r[1]) for r in radiance_rows(case)}, key=lambda x: (-x[0], x[1]))
        fig, axes = plt.subplots(len(labels), len(panels), figsize=(3.2 * len(panels), 2.6 * len(labels)),
                                 squeeze=False, sharex=True, sharey="row", constrained_layout=True)
        for col, (z, phi) in enumerate(panels):
            rows = [r for r in radiance_rows(case) if r[0] == z and r[1] == phi]
            mu_t = np.array([r[2] for r in rows])
            level, k = level_of(case, z), int(round(phi / 90))
            if "vdisort" in res:
                mu_v, u_v = res["vdisort"]
                n = len(mu_v) // 2
                x = np.concatenate((-mu_v[:n][::-1], mu_v[:n]))  # Evans' sign: upwelling negative
                vals = np.concatenate((u_v[level, k, :n][::-1], u_v[level, k, n:]))
                y = np.array([as_table(case, (z, phi, xx, np.zeros(4)), vv) for xx, vv in zip(x, vals)])
            for s, label in enumerate(labels):
                ax = axes[s, col]
                if "vdisort" in res:
                    ax.plot(x, y[:, s], "-", color="C0", lw=1.5, label="VDISORT (64 double-Gauss streams)")
                ax.plot(mu_t, [r[3][s] for r in rows], "o", mfc="none", color="k", ms=7, label="Evans' table")
                arts_values = [as_table(case, r, res["arts"][1](level, k, r[2]))[s] for r in rows]
                ax.plot(mu_t, arts_values, ".", color="C3", ms=5, label=f"ARTS {case['program'].upper()} at Evans' settings")
                ax.grid(True, alpha=0.3)
                if col == 0:
                    ax.set_ylabel(label)
                if s == 0:
                    ax.set_title(f"Z = {z:g} km" + ("" if case["brightness"] else f", phi = {phi:g} deg"))
                if s == len(labels) - 1:
                    ax.set_xlabel("mu (< 0 upwelling, as Evans prints it)")
        axes[0, 0].legend(fontsize=7)
        fig.suptitle(f"{name}: solutions")

        if "sweep" in res:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
            nmu, d = np.array(res["sweep"]).T
            ax1.loglog(nmu, d, "o-", label="ARTS RT3, Gauss quadrature, vs VDISORT")
            ax1.loglog(nmu, d[0] * nmu[0] / nmu, "k:", label="1 / nmu")
            ax1.axvline(case["nmu"], color="gray", ls="--", label=f"Evans' table (nmu {case['nmu']})")
            if case["brightness"]:
                ax1.text(0.03, 0.05, f"RT3 and VDISORT on the same double-Gauss streams: {res['same_streams']:.1e} of max I",
                         transform=ax1.transAxes, color="C2", fontsize=8)
            else:
                ax1.axhline(res["same_streams"], color="C2", label="RT3 and VDISORT on the same double-Gauss streams")
            ax1.set(xlabel="RT3 streams per hemisphere", ylabel=f"max difference ({'K' if case['brightness'] else 'of max I'})",
                    title=f"{name}: the table is a Gauss-quadrature solution")
            ax1.legend(fontsize=7)
            ax1.grid(True, which="both", alpha=0.3)
            rows = radiance_rows(case)
            diffs = []
            for r in rows:
                level, k = level_of(case, r[0]), int(round(r[1] / 90))
                v = as_table(case, r, vdisort_at(res["vdisort"], [r[2]], level, k)[0])
                diffs.append((r[2], (r[3][: len(v)] - v)[0]))
            m, dd = np.array(diffs).T
            ax2.plot(m, dd, "o", mfc="none")
            ax2.set(xlabel="mu (< 0 upwelling)", ylabel=f"table - VDISORT, {labels[0]}",
                    title="Evans' table minus VDISORT, every level and azimuth")
            ax2.grid(True, alpha=0.3)
    plt.show()
