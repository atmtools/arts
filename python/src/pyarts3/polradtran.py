"""VDISORT and Evans' polradtran (RT3, RT4) on the same problem.

This module defines a plane-parallel problem once and runs it through every
solver that can represent it, the same way:

* VDISORT (:mod:`pyarts3.arts.vdisort`), always;
* RT4 (:mod:`pyarts3.arts.rt4`, Evans' polarized doubling-adding for
  azimuthally symmetric media) when the problem is thermal only, without a
  beam, with ``nstokes <= 2``;
* RT3 (:mod:`pyarts3.arts.rt3`, Evans' doubling-adding with a solar beam and
  all Fourier azimuth modes) unless a beam falls on a non-Lambertian
  surface.

RT3 and RT4 are optional (``ENABLE_RT3`` and ``ENABLE_RT4``); a solver whose
backend is not built is skipped.  All solvers get the same double-Gauss
streams and the same phase matrix.  RT3 is given the Legendre series of the
scattering matrix directly.  VDISORT and RT4 are given its stream-pair
Fourier modes, from the lab-frame Mueller matrix built from vector geometry
and sampled at RT3's own azimuths (:func:`rt3_azimuth_samples`).
The solvers therefore solve the same semi-discrete problem, and RT3 and RT4
differ from VDISORT only by their first-order doubling error.  That error is
below ``max_delta_tau / mu0`` of the largest radiance (``mu0 = 1`` without a
beam), and :meth:`Comparison.tolerance` is 10 times it.

Conventions (those of :mod:`pyarts3.arts.rt3`):

* Layers and levels are ordered top-down; level 0 is the top.
* Directions are ``k = (sin cos phi, sin sin phi, mu)`` with ``z`` up.  The
  beam propagates downward toward azimuth 0, so ``phi = 0`` is the forward
  half for upwelling radiation.
* The Stokes basis is meridional: ``h = k x z / |k x z|``, ``v = h x k``,
  ``Q = I_v - I_h``, ``U = 2 Re(E_v E_h*)``.
* Radiances are in W m-2 Hz-1 sr-1, fluxes in W m-2 Hz-1.
* ``nstokes < 4`` means the problem truncated to the first ``nstokes``
  Stokes components (for ``nstokes = 1`` the scalar problem).  VDISORT,
  which always carries four components, is given that truncated problem.

Example
-------

.. code-block:: python

    import pyarts3 as pyarts

    comparison = pyarts.polradtran.compare(pyarts.polradtran.cases.solar_rayleigh())
    print(comparison)
    fig, ax = pyarts.polradtran.plot(comparison, level=0, azimuth=90.0)

The C++ comparisons ``src/core/disort-cpp/test/vdisort/vdisort-rt{3,4}-comparison.cpp``
measure the doubling error, its convergence and Richardson extrapolation, and
show that deliberate input mistakes are detected.
"""

from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np
from pyarts3 import arts

__all__ = [
    "EVANS_MIE",
    "RAYLEIGH",
    "Comparison",
    "Problem",
    "Radiance",
    "ScatteringSet",
    "applicable_solvers",
    "available_solvers",
    "cases",
    "compare",
    "fourier_modes",
    "legendre_series",
    "planck",
    "plot",
    "rt3_azimuth_samples",
    "solve",
]

_H, _C, _K = 6.62607015e-34, 299792458.0, 1.380649e-23

#: Rayleigh scattering as RT3 Legendre columns (F11, F12, F33, F34, F22, F44)
RAYLEIGH = np.array(
    [
        [1.0, -0.5, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.5, 0.0, 0.0, 1.5],
        [0.5, 0.5, 0.0, 0.0, 0.5, 0.0],
    ]
)

#: Evans' Mie scattering matrix of polradtran's ``runmietest`` (RT3 Legendre columns)
EVANS_MIE = np.array(
    [
        [1.00000000, -0.32071711, 0.71206342, -0.01882245, 1.00000000, 0.71206342],
        [1.45529318, -0.20350675, 1.76014119, -0.04725108, 1.45529318, 1.76014119],
        [1.05402631, 0.24638948, 1.06682431, 0.00894436, 1.05402631, 1.06682431],
        [0.39758994, 0.18605748, 0.39651104, 0.04505815, 0.39758994, 0.39651104],
        [0.11659302, 0.07124848, 0.09576412, 0.00958275, 0.11659302, 0.09576412],
        [0.02387477, 0.01700757, 0.01765088, 0.00215761, 0.02387477, 0.01765088],
        [0.00395010, 0.00302534, 0.00261549, 0.00029195, 0.00395010, 0.00261549],
        [0.00053888, 0.00043592, 0.00032713, 0.00003502, 0.00053888, 0.00032713],
        [0.00006372, 0.00005326, 0.00003583, 0.00000337, 0.00006372, 0.00003583],
        [0.00000667, 0.00000572, 0.00000351, 0.00000029, 0.00000667, 0.00000351],
        [0.00000063, 0.00000055, 0.00000031, 0.00000002, 0.00000063, 0.00000031],
        [0.00000006, 0.00000005, 0.00000003, 0.00000000, 0.00000006, 0.00000003],
    ]
)


def planck(frequency, temperature):
    """Planck radiance B_nu(T) [W m-2 Hz-1 sr-1]; 0 for T <= 0."""
    if temperature <= 0:
        return 0.0
    return (
        2.0 * _H * frequency**3 / _C**2 / np.expm1(_H * frequency / (_K * temperature))
    )


@dataclass
class ScatteringSet:
    """Particle optics of a homogeneous layer.

    Attributes
    ----------
    extinction : float
        Particle extinction coefficient per unit length.
    scattering : float
        Particle scattering coefficient per unit length.
    legendre : numpy.ndarray
        ``[nleg + 1, 6]`` scattering matrix as plain Legendre series in
        ``cos(Theta)`` with RT3's column order (F11, F12, F33, F34, F22, F44),
        including the factor ``2 l + 1`` and with ``legendre[0, 0] == 1``.
    """

    extinction: float
    scattering: float
    legendre: np.ndarray


@dataclass
class Problem:
    """A plane-parallel problem that every applicable solver runs the same way.

    Attributes
    ----------
    name : str
        A label for printing and plot titles.
    nstokes : int
        1 to 4.
    nmu : int
        Double-Gauss streams per hemisphere.
    nfourier : int
        Fourier azimuth modes (1 for an azimuthally symmetric problem).
    frequency : float
        [Hz].
    dz : list of float
        Layer thicknesses, top-down, in the inverse unit of the extinctions.
    gas : list of float
        Unpolarized gas extinction per layer.
    sets : list of ScatteringSet
        Scattering properties.
    set_index : list of int
        Per layer, the index into ``sets``, or -1 for a gas-only layer.
    temperature : list of float
        [K] at the ``len(dz) + 1`` interfaces, top-down.  Planck is linear in
        optical depth within each layer.
    sky : float
        [K] of the isotropic unpolarized blackbody incident at the top.
    surface_temperature : float
        [K].
    surface : tuple
        ``("lambertian", albedo)`` or ``("fresnel", refractive_index)``.
    thermal : bool
        Thermal emission of the layers and of a Lambertian surface.  The sky
        and a Fresnel surface always emit (as in RT3); set their temperatures
        to 0 K to remove them.
    beam : tuple or None
        ``(flux on the horizontal plane [W m-2 Hz-1], mu0)`` of a solar beam
        propagating toward azimuth 0, or None.
    """

    name: str
    nstokes: int
    nmu: int
    nfourier: int
    frequency: float
    dz: list
    gas: list
    sets: list
    set_index: list
    temperature: list
    sky: float
    surface_temperature: float
    surface: tuple
    thermal: bool = True
    beam: tuple | None = None
    height: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.height = np.concatenate((np.cumsum(self.dz[::-1])[::-1], [0.0]))

    @property
    def nlay(self):
        return len(self.dz)

    @property
    def mu0(self):
        """Cosine of the beam zenith angle, 1 without a beam."""
        return self.beam[1] if self.beam else 1.0

    def extinction(self, layer):
        i = self.set_index[layer]
        return self.gas[layer] + (self.sets[i].extinction if i >= 0 else 0.0)

    def scattering(self, layer):
        i = self.set_index[layer]
        return self.sets[i].scattering if i >= 0 else 0.0


@dataclass
class Radiance:
    """One solver's solution on the shared streams.

    Attributes
    ----------
    solver : str
    mu : numpy.ndarray
        ``[nmu]`` ascending stream cosines, the same in both hemispheres.
    phi : numpy.ndarray
        ``[nphi]`` azimuths [rad].
    up, down : numpy.ndarray
        ``[nlay + 1, nphi, nmu, nstokes]`` radiance propagating upward and
        downward at each level (0 = top).
    """

    solver: str
    mu: np.ndarray
    phi: np.ndarray
    up: np.ndarray
    down: np.ndarray


# --------------------------------------------------------------------------
# Phase matrix on the streams
# --------------------------------------------------------------------------
def _direction(mu, phi):
    s = np.sqrt(np.maximum(0.0, 1.0 - mu**2))
    return np.stack(np.broadcast_arrays(s * np.cos(phi), s * np.sin(phi), mu), -1)


def _scattering_matrix(legendre, cos_theta):
    f = [np.polynomial.legendre.legval(cos_theta, legendre[:, c]) for c in range(6)]
    F = np.zeros(cos_theta.shape + (4, 4))
    (
        F[..., 0, 0],
        F[..., 0, 1],
        F[..., 2, 2],
        F[..., 2, 3],
        F[..., 1, 1],
        F[..., 3, 3],
    ) = f
    F[..., 1, 0] = F[..., 0, 1]
    F[..., 3, 2] = -F[..., 2, 3]
    return F


def _stokes_rotation(eta):
    c, s = np.cos(2 * eta), np.sin(2 * eta)
    M = np.zeros(eta.shape + (4, 4))
    M[..., 0, 0] = M[..., 3, 3] = 1.0
    M[..., 1, 1] = M[..., 2, 2] = c
    M[..., 1, 2], M[..., 2, 1] = s, -s
    return M


def _lab_phase(k_out, k_in, legendre):
    def meridional(k):
        h = np.cross(k, [0.0, 0.0, 1.0])
        h /= np.linalg.norm(h, axis=-1)[..., None]
        return np.cross(h, k), h

    v_in, h_in = meridional(k_in)
    v_out, h_out = meridional(k_out)
    perp = np.cross(k_in, k_out)
    norm = np.linalg.norm(perp, axis=-1)
    parallel_rays = norm < 1e-12
    perp = np.where(
        parallel_rays[..., None],
        h_in,
        perp / np.where(parallel_rays, 1.0, norm)[..., None],
    )
    par_in, par_out = np.cross(perp, k_in), np.cross(perp, k_out)
    eta_in = np.arctan2(np.sum(par_in * h_in, -1), np.sum(par_in * v_in, -1))
    eta_out = np.arctan2(np.sum(par_out * h_out, -1), np.sum(par_out * v_out, -1))
    F = _scattering_matrix(legendre, np.clip(np.sum(k_out * k_in, -1), -1.0, 1.0))
    return _stokes_rotation(-eta_out) @ F @ _stokes_rotation(eta_in)


def fourier_modes(legendre, mu_out, mu_in, nfourier, nphi=1024):
    """Ordinary Fourier coefficients of the lab-frame phase matrix.

    ``C^m, S^m = (1 / 2 pi) int Z(mu_out, 0; mu_in, phi') {cos, sin}(m phi') dphi'``
    on signed stream cosines (> 0 upward), without the factor ``2 - delta_m0``.
    The phase matrix is normalized to 1 over 4 pi.

    Returns
    -------
    C, S : numpy.ndarray
        ``[nfourier, len(mu_out), len(mu_in), 4, 4]``.
    """
    mu_out, mu_in = np.asarray(mu_out, dtype=float), np.asarray(mu_in, dtype=float)
    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    k_out, k_in = np.broadcast_arrays(
        _direction(mu_out[:, None, None], 0.0), _direction(mu_in[None, :, None], phi)
    )
    Z = _lab_phase(k_out, k_in, legendre)
    C = np.stack(
        [np.einsum("oipab,p->oiab", Z, np.cos(m * phi)) / nphi for m in range(nfourier)]
    )
    S = np.stack(
        [np.einsum("oipab,p->oiab", Z, np.sin(m * phi)) / nphi for m in range(nfourier)]
    )
    return C, S


def rt3_azimuth_samples(legendre, nfourier):
    """The number of equidistant azimuths at which RT3 samples the phase matrix.

    RT3 (``SCATTERING`` in radscat3.f) takes the Fourier modes of the
    lab-frame matrix by the trapezoid rule at ``2 pi k / n``.  Using the same
    samples for VDISORT and RT4 gives every solver the identical discrete
    phase matrix, also for a scattering matrix that is not regular at
    ``Theta = 0`` or ``pi`` (then the modes depend on the samples).
    """
    nonzero = np.flatnonzero(np.any(np.asarray(legendre) != 0.0, axis=1))
    degree = int(nonzero.max()) if nonzero.size else 0
    if nfourier <= 1:
        return 2 * ((degree + 1) // 2) + 4
    return 2 * 2 ** int(np.log(degree + 4.0) / np.log(2.0) + 1.0)


def legendre_series(scattering_matrix, degree, nquad=512):
    """RT3 Legendre columns of a scattering matrix given as a function.

    Parameters
    ----------
    scattering_matrix : callable
        ``f(cos_theta) -> (F11, F12, F33, F34, F22, F44)`` in the
        scattering-plane basis (``Q = I_par - I_perp``), arrays like
        ``cos_theta``.  For a physical matrix, F12 = F34 = 0 and F22 = F33 at
        ``Theta = 0``, and F22 = -F33 at ``Theta = pi``.
    degree : int
        Highest Legendre degree kept.  RT3 with double-Gauss streams keeps at
        most ``2 nmu - 3``.
    nquad : int, optional
        Gauss-Legendre nodes for the projection.

    Returns
    -------
    numpy.ndarray
        ``[degree + 1, 6]``, normalized so that ``legendre[0, 0] == 1``.
    """
    x, w = np.polynomial.legendre.leggauss(nquad)
    F = np.stack(scattering_matrix(x), axis=-1)
    P = np.stack(
        [np.polynomial.legendre.Legendre.basis(l)(x) for l in range(degree + 1)]
    )
    coef = (2 * np.arange(degree + 1)[:, None] + 1) / 2 * (P * w) @ F
    return coef / coef[0, 0]


def _streams(nmu):
    x, w = np.polynomial.legendre.leggauss(nmu)
    mu = 0.5 * (x + 1.0)
    return mu, 0.5 * w, np.concatenate((mu, -mu))  # VDISORT: upward first


def _fresnel(mu, n):
    cos_t = np.sqrt(1.0 - (1.0 - mu**2) / n**2 + 0j)
    rv = (n * mu - cos_t) / (n * mu + cos_t)
    rh = (mu - n * cos_t) / (mu + n * cos_t)
    return 0.5 * (abs(rv) ** 2 + abs(rh) ** 2), 0.5 * (abs(rv) ** 2 - abs(rh) ** 2)


# --------------------------------------------------------------------------
# The solvers
# --------------------------------------------------------------------------
def _vdisort(p, phi, max_delta_tau):
    vd = arts.vdisort
    mu, _, signed = _streams(p.nmu)
    N, NF, NL, ns = p.nmu, p.nfourier, p.nlay, p.nstokes
    tau = np.cumsum([p.extinction(l) * p.dz[l] for l in range(NL)])
    omega = np.array([p.scattering(l) / p.extinction(l) for l in range(NL)])

    phase = np.zeros((2, NF, NL, 2 * N, 2 * N, 4, 4))
    beam_phase = np.zeros((2, NF, NL, 2 * N, 4, 4))
    for i, s in enumerate(p.sets):
        layers = [l for l in range(NL) if p.set_index[l] == i]
        if not layers:
            continue
        nphi = rt3_azimuth_samples(s.legendre, NF)
        C, S = fourier_modes(s.legendre, signed, signed, NF, nphi)
        phase[:, :, layers] = np.asarray(
            vd.combine_phase_matrices(C[:, None], S[:, None])
        )[:, :, 0][:, :, None]
        if p.beam:
            Cb, Sb = fourier_modes(s.legendre, signed, np.array([-p.mu0]), NF, nphi)
            combined = vd.combine_beam_phase_matrices(
                Cb[:, None, :, 0], Sb[:, None, :, 0]
            )
            beam_phase[:, :, layers] = np.asarray(combined)[:, :, 0][:, :, None]
    # The problem truncated to nstokes components, as RT3 and RT4 solve it
    phase[..., ns:, :] = phase[..., :, ns:] = 0.0
    beam_phase[..., ns:, :] = beam_phase[..., :, ns:] = 0.0

    top, bottom = np.zeros((2, NF, N, 4)), np.zeros((2, NF, N, 4))
    top[vd.cosine_mode, 0, :, 0] = planck(p.frequency, p.sky)
    Bs = planck(p.frequency, p.surface_temperature)
    if p.surface[0] == "fresnel":
        R1, R2 = _fresnel(mu, p.surface[1])
        bottom[vd.cosine_mode, 0, :, 0], bottom[vd.cosine_mode, 0, :, 1] = (
            1 - R1
        ) * Bs, -R2 * Bs
        brdf = vd.fresnel_fourier_modes(p.surface[1], NF)
    else:
        bottom[vd.cosine_mode, 0, :, 0] = (1 - p.surface[1]) * Bs if p.thermal else 0.0
        brdf = vd.lambertian_fourier_modes(p.surface[1], NF)
    bottom[..., ns:] = 0.0

    source = np.zeros((NL, 2, 4))
    if p.thermal:
        tau_top = np.concatenate(([0.0], tau[:-1]))
        B = np.array([planck(p.frequency, t) for t in p.temperature])
        slope = (B[1:] - B[:-1]) / (tau - tau_top)
        source[:, 0, 0], source[:, 1, 0] = B[:-1] - slope * tau_top, slope

    flux, mu0 = p.beam if p.beam else (0.0, 0.5)
    model = arts.cppvdisort(
        tau_arr=tau,
        omega_arr=omega,
        NQuad=2 * N,
        NFourier=NF,
        phase_matrix=phase,
        mu0=mu0,
        beam_stokes=np.array([flux / mu0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        b_pos=bottom,
        b_neg=top,
        BDRF_Fourier_modes=brdf,
        s_poly_coeffs=source,
        beam_phase_matrix=beam_phase,
    )
    u = np.asarray(model.u(tau=np.concatenate(([0.0], tau)), phi=phi))
    return Radiance("VDISORT", mu, phi, u[:, :, :N, :ns], u[:, :, N:, :ns])


def _rt4(p, phi, max_delta_tau):
    rt4 = arts.rt4
    mu, _, signed = _streams(p.nmu)
    N, ns = p.nmu, p.nstokes

    def stream(h, i):
        return N + i if h == rt4.down else i

    optics = []
    for s in p.sets:
        C0 = fourier_modes(
            s.legendre, signed, signed, 1, rt3_azimuth_samples(s.legendre, 1)
        )[0][0]
        Z = np.zeros((2, 2, N, N, ns, ns))
        for ho in (rt4.down, rt4.up):
            for hi in (rt4.down, rt4.up):
                for io in range(N):
                    for ii in range(N):
                        Z[ho, hi, io, ii] = (
                            s.scattering
                            / (4 * np.pi)
                            * C0[stream(ho, io), stream(hi, ii)][:ns, :ns]
                        )
        ext = np.zeros((2, N, ns, ns))
        ext[..., range(ns), range(ns)] = s.extinction
        absorption = np.zeros((2, N, ns))
        absorption[..., 0] = s.extinction - s.scattering
        optics.append(rt4.LayerOptics(ext, absorption, Z))
    ground = (
        rt4.FresnelSurface(p.surface[1])
        if p.surface[0] == "fresnel"
        else rt4.LambertianSurface(p.surface[1])
    )
    r = rt4.solve(
        rt4.Problem(
            nstokes=ns,
            nmu=N,
            quad=rt4.QuadratureType.double_gauss,
            max_delta_tau=max_delta_tau,
            frequency=p.frequency,
            height=p.height,
            temperature=np.asarray(p.temperature, dtype=float),
            gas_extinction=np.asarray(p.gas, dtype=float),
            optics=optics,
            layer_optics_index=p.set_index,
            sky_temperature=p.sky,
            surface_temperature=p.surface_temperature,
            ground=ground,
        )
    )
    # Fourier mode 0 only: the same at every azimuth
    up, down = (
        np.repeat(np.asarray(x)[:, None], len(phi), axis=1) for x in (r.up, r.down)
    )
    return Radiance("RT4", mu, phi, up, down)


def _rt3(p, phi, max_delta_tau):
    rt3 = arts.rt3
    mu, _, _ = _streams(p.nmu)
    ground = (
        rt3.FresnelSurface(p.surface[1])
        if p.surface[0] == "fresnel"
        else rt3.LambertianSurface(p.surface[1])
    )
    flux, mu0 = p.beam if p.beam else (0.0, 1.0)
    r = rt3.solve(
        rt3.Problem(
            nstokes=p.nstokes,
            nmu=p.nmu,
            quad=rt3.QuadratureType.double_gauss,
            aziorder=p.nfourier - 1,
            max_delta_tau=max_delta_tau,
            direct_flux=flux,
            direct_mu=mu0,
            thermal=p.thermal,
            frequency=p.frequency,
            height=p.height,
            temperature=np.asarray(p.temperature, dtype=float),
            gas_extinction=np.asarray(p.gas, dtype=float),
            scattering_sets=[
                rt3.ScatteringSet(s.extinction, s.scattering, s.legendre)
                for s in p.sets
            ],
            layer_scattering_index=p.set_index,
            sky_temperature=p.sky,
            surface_temperature=p.surface_temperature,
            ground=ground,
        )
    )
    up, down = (np.asarray(rt3.azimuth_radiance(x, phi)) for x in (r.up, r.down))
    return Radiance("RT3", mu, phi, up, down)


_SOLVERS = {"VDISORT": _vdisort, "RT4": _rt4, "RT3": _rt3}


def available_solvers():
    """The solvers built into this pyarts: VDISORT, plus RT4 and RT3 when enabled."""
    return [
        name
        for name, ok in (
            ("VDISORT", True),
            ("RT4", arts.rt4.available()),
            ("RT3", arts.rt3.available()),
        )
        if ok
    ]


def applicable_solvers(problem):
    """The available solvers that can represent ``problem``, VDISORT first."""
    can = {
        "VDISORT": True,
        "RT4": problem.beam is None and problem.nstokes <= 2,
        "RT3": problem.beam is None or problem.surface[0] == "lambertian",
    }
    return [name for name in available_solvers() if can[name]]


def solve(problem, solver, *, phi=None, max_delta_tau=1e-7):
    """Run one solver.

    Parameters
    ----------
    problem : Problem
    solver : str
        "VDISORT", "RT4" or "RT3".
    phi : array_like, optional
        Azimuths [rad].  Defaults to 0 to 360 degrees in steps of 15.
    max_delta_tau : float, optional
        Initial doubling-layer thickness of RT3 and RT4 (unused by VDISORT).

    Returns
    -------
    Radiance
    """
    if solver not in applicable_solvers(problem):
        raise ValueError(
            f"{solver} cannot run {problem.name!r} in this build; applicable: {applicable_solvers(problem)}"
        )
    phi = (
        np.radians(np.arange(0.0, 360.0, 15.0))
        if phi is None
        else np.asarray(phi, dtype=float)
    )
    return _SOLVERS[solver](problem, phi, max_delta_tau)


@dataclass
class Comparison:
    """Solutions of one problem by several solvers.

    Attributes
    ----------
    problem : Problem
    results : dict[str, Radiance]
        By solver name, VDISORT first.
    max_delta_tau : float
        The initial doubling-layer thickness RT3 and RT4 were run with.
    """

    problem: Problem
    results: dict
    max_delta_tau: float

    @property
    def scale(self):
        """Largest |I| of the reference (VDISORT) solution, the unit of :meth:`deviation`."""
        r = self.results["VDISORT"]
        return max(np.abs(r.up[..., 0]).max(), np.abs(r.down[..., 0]).max())

    def tolerance(self):
        """10 times the bound on RT3's and RT4's doubling error, relative to :attr:`scale`."""
        return 10.0 * self.max_delta_tau / self.problem.mu0 + 1e-12

    def deviation(self, a, b):
        """Max |a - b| over levels, azimuths, streams and both directions, per Stokes, / :attr:`scale`."""
        ra, rb = self.results[a], self.results[b]
        return (
            np.array(
                [
                    max(
                        np.abs(ra.up[..., s] - rb.up[..., s]).max(),
                        np.abs(ra.down[..., s] - rb.down[..., s]).max(),
                    )
                    for s in range(self.problem.nstokes)
                ]
            )
            / self.scale
        )

    def pairs(self):
        names = list(self.results)
        return [
            (names[i], names[j])
            for i in range(len(names))
            for j in range(i + 1, len(names))
        ]

    def polarization(self):
        """max |Q|, |U|, |V| over max |I| of the reference solution."""
        r = self.results["VDISORT"]
        return (
            np.array(
                [
                    max(np.abs(r.up[..., s]).max(), np.abs(r.down[..., s]).max())
                    for s in range(1, self.problem.nstokes)
                ]
            )
            / self.scale
        )

    def check(self):
        """Raise AssertionError unless every pair agrees to :meth:`tolerance`."""
        for a, b in self.pairs():
            dev = self.deviation(a, b).max()
            if not dev <= self.tolerance():
                raise AssertionError(
                    f"{self.problem.name}: {a} and {b} must agree to {self.tolerance():.1e} of max I "
                    f"(max_delta_tau = {self.max_delta_tau}), got {dev:.3e}"
                )

    def __str__(self):
        lines = [
            (
                f"{self.problem.name}: {', '.join(self.results)}; "
                f"max |Q|, |U|, |V| / max |I| = {np.round(self.polarization(), 4)}"
            )
        ]
        for a, b in self.pairs():
            lines.append(
                f"    {a:7s} vs {b:7s}  max |diff| / max |I| (I, Q, U, V) = {self.deviation(a, b)}"
                f"  tolerance {self.tolerance():.1e}"
            )
        return "\n".join(lines)


def compare(problem, *, solvers=None, phi=None, max_delta_tau=1e-7):
    """Run ``problem`` through ``solvers`` (default: every applicable one).

    Returns
    -------
    Comparison
    """
    solvers = applicable_solvers(problem) if solvers is None else list(solvers)
    if "VDISORT" not in solvers:
        solvers = ["VDISORT"] + solvers
    return Comparison(
        problem,
        {s: solve(problem, s, phi=phi, max_delta_tau=max_delta_tau) for s in solvers},
        max_delta_tau,
    )


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
_STOKES = ["I", "Q", "U", "V"]
_STYLE = {
    "VDISORT": {"color": "C0", "ls": "-", "marker": "."},
    "RT4": {"color": "C1", "ls": "", "marker": "o", "mfc": "none"},
    "RT3": {"color": "C2", "ls": "", "marker": "x"},
}


def plot(comparison, *, level=0, azimuth=0.0, fig=None, ax=None, stokes=None, **kwargs):
    """Plot the solvers' radiances and their differences from VDISORT.

    The top row shows each Stokes component against the line-of-sight zenith
    angle (ARTS convention: below 90 degrees looks up and sees the downward
    radiance, above 90 looks down and sees the upward radiance) at one level
    and azimuth.  The bottom row shows ``|solver - VDISORT| / max |I|`` on a
    logarithmic scale, with the tolerance of :meth:`Comparison.tolerance`
    (the comparison passes when every point is below it).

    .. rubric:: Example

    .. code-block:: python

        import pyarts3 as pyarts

        comparison = pyarts.polradtran.compare(pyarts.polradtran.cases.thermal_rayleigh())
        fig, ax = pyarts.polradtran.plot(comparison, level=0)

    Parameters
    ----------
    comparison : Comparison
    level : int, optional
        Level index (0 = top of the atmosphere, ``nlay`` = just above the surface).
    azimuth : float, optional
        Azimuth [degree] of the line of sight's propagation direction; the
        nearest computed azimuth is used.
    fig : ~matplotlib.figure.Figure, optional
    ax : numpy.ndarray of ~matplotlib.axes.Axes, optional
        ``[2, len(stokes)]`` axes to draw on.
    stokes : list of int, optional
        Stokes components to show; defaults to all ``nstokes``.
    **kwargs
        Passed to :func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    fig : ~matplotlib.figure.Figure
    ax : numpy.ndarray of ~matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    stokes = list(range(comparison.problem.nstokes)) if stokes is None else list(stokes)
    if fig is None:
        fig = plt.figure(figsize=(3.6 * len(stokes), 6.0), constrained_layout=True)
    if ax is None:
        ax = np.asarray(fig.subplots(2, len(stokes), sharex=True, squeeze=False))

    ref = comparison.results["VDISORT"]
    k = int(np.argmin(np.abs((np.degrees(ref.phi) - azimuth + 180.0) % 360.0 - 180.0)))
    # A NaN at 90 degrees separates the two hemispheres, which are not continuous there
    za = np.concatenate(
        (
            np.degrees(np.arccos(ref.mu[::-1])),
            [90.0],
            180.0 - np.degrees(np.arccos(ref.mu)),
        )
    )

    def los(r, s):
        # looking up (za < 90) sees down[mu descending]; looking down sees up[mu ascending]
        return np.concatenate(
            (r.down[level, k, ::-1, s], [np.nan], r.up[level, k, :, s])
        )

    tol = comparison.tolerance()
    for col, s in enumerate(stokes):
        top, bottom = ax[0, col], ax[1, col]
        for name, r in comparison.results.items():
            top.plot(za, los(r, s), label=name, **{**_STYLE.get(name, {}), **kwargs})
            if name != "VDISORT":
                diff = np.abs(los(r, s) - los(ref, s)) / comparison.scale
                bottom.semilogy(
                    za,
                    np.maximum(diff, 1e-18),
                    label=f"{name} - VDISORT",
                    **{**_STYLE.get(name, {}), **kwargs},
                )
        if len(comparison.results) > 1:
            bottom.axhline(tol, color="k", ls="--", lw=0.8, label="tolerance")
        top.set_title(f"Stokes {_STOKES[s]}")
        top.axvline(90.0, color="0.7", lw=0.5)
        bottom.axvline(90.0, color="0.7", lw=0.5)
        bottom.set_xlabel("Line-of-sight zenith angle [deg]")
        top.grid(True, alpha=0.3)
        bottom.grid(True, alpha=0.3)
    ax[0, 0].set_ylabel("Radiance [W m$^{-2}$ Hz$^{-1}$ sr$^{-1}$]")
    ax[1, 0].set_ylabel("|difference| / max |I|")
    ax[0, 0].legend()
    if len(comparison.results) > 1:
        ax[1, 0].legend()
    fig.suptitle(
        f"{comparison.problem.name}\nlevel {level}, azimuth {np.degrees(ref.phi[k]):.0f} deg, "
        f"max_delta_tau = {comparison.max_delta_tau:g}"
    )
    return fig, ax


# --------------------------------------------------------------------------
# Preset problems (those run in CI by tests/core/disort/vdisort-polradtran.rt3.rt4.py)
# --------------------------------------------------------------------------
def _thermal_rayleigh(nstokes=2, nmu=8):
    """Thermal Rayleigh multilayer over a Fresnel surface (VDISORT, RT4, RT3)."""
    return Problem(
        name=f"thermal Rayleigh multilayer over Fresnel 3+0.2i, nstokes {nstokes}",
        nstokes=nstokes,
        nmu=nmu,
        nfourier=1,
        frequency=89e9,
        dz=[1.0, 2.0, 0.5, 1.0],
        gas=[0.05, 0.0, 0.2, 0.1],
        sets=[ScatteringSet(0.5, 0.45, RAYLEIGH), ScatteringSet(1.0, 1.0, RAYLEIGH)],
        set_index=[0, 1, -1, 0],
        temperature=[210.0, 230.0, 255.0, 270.0, 285.0],
        sky=2.725,
        surface_temperature=290.0,
        surface=("fresnel", 3.0 + 0.2j),
    )


def _thermal_mie(nmu=8):
    """Thermal layer of Evans' Mie particles over a Lambertian surface (VDISORT, RT4, RT3)."""
    return Problem(
        name="thermal Mie (Evans) layer over Lambertian 0.3",
        nstokes=2,
        nmu=nmu,
        nfourier=1,
        frequency=89e9,
        dz=[0.5, 1.5],
        gas=[0.1, 0.02],
        sets=[ScatteringSet(1.0, 0.9, EVANS_MIE)],
        set_index=[-1, 0],
        temperature=[220.0, 250.0, 280.0],
        sky=2.725,
        surface_temperature=285.0,
        surface=("lambertian", 0.3),
    )


def _solar_rayleigh(nmu=8):
    """Solar Rayleigh layer over a Lambertian surface, all four Stokes components (VDISORT, RT3)."""
    return Problem(
        name="solar Rayleigh layer over Lambertian 0.1, 4 Fourier modes",
        nstokes=4,
        nmu=nmu,
        nfourier=4,
        frequency=_C / 3e-6,
        dz=[1.0],
        gas=[0.0],
        sets=[ScatteringSet(0.5, 0.475, RAYLEIGH)],
        set_index=[0],
        temperature=[1.0, 1.0],
        sky=0.0,
        surface_temperature=0.0,
        surface=("lambertian", 0.1),
        thermal=False,
        beam=(1.0, 0.6),
    )


def _solar_thermal_multilayer(nmu=8):
    """Solar and thermal sources at 3 um, Rayleigh / Mie / gas layers (VDISORT, RT3)."""
    return Problem(
        name="solar + thermal at 3 um, Rayleigh / Mie / gas over Lambertian 0.25, 12 Fourier modes",
        nstokes=4,
        nmu=nmu,
        nfourier=12,
        frequency=_C / 3e-6,
        dz=[10.0, 5.0, 5.0],
        gas=[0.002, 0.01, 0.005],
        sets=[
            ScatteringSet(0.05, 0.05, RAYLEIGH),
            ScatteringSet(0.2, 0.198, EVANS_MIE),
        ],
        set_index=[0, 1, -1],
        temperature=[200.0, 270.0, 290.0, 300.0],
        sky=0.0,
        surface_temperature=300.0,
        surface=("lambertian", 0.25),
        beam=(5.0e-15, 0.5),
    )


def _all():
    """Every preset problem, as run in CI."""
    return [
        _thermal_rayleigh(2),
        _thermal_rayleigh(1),
        _thermal_mie(),
        _solar_rayleigh(),
        _solar_thermal_multilayer(),
    ]


#: Preset problems: ``thermal_rayleigh(nstokes)``, ``thermal_mie()``,
#: ``solar_rayleigh()``, ``solar_thermal_multilayer()`` and ``all()``.
cases = SimpleNamespace(
    thermal_rayleigh=_thermal_rayleigh,
    thermal_mie=_thermal_mie,
    solar_rayleigh=_solar_rayleigh,
    solar_thermal_multilayer=_solar_thermal_multilayer,
    all=_all,
)
