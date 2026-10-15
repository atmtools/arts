"""VDISORT, RT4 and RT3 on ARTS atmospheres, through their ARTS-native builders.

Compares the polarized discrete-ordinate solver VDISORT with Evans'
doubling-adding solvers RT4 (thermal, azimuthally symmetric, [I, Q]) and RT3
(solar beam and all Fourier azimuth modes, [I, Q, U, V]) of polradtran, which
ARTS builds when a Fortran compiler is found.

Every problem is ARTS data: a propagation path (one level per entry, top
first), its atmospheric points, the gas propagation matrix, ARTS scattering
species (a Rayleigh GasScatterer, a HenyeyGreensteinScatterer and liquid Mie
spheres) and a surface.  vdisort.main_data_from_path, rt4.problem_from_path
and rt3.problem_from_path turn them into each solver's input, with the same
double-Gauss streams, so the solvers solve the same problem.  RT3 and RT4
then differ from VDISORT by their first-order doubling error, below
max_delta_tau / mu0 of max I (mu0 = 1 without a beam), and every pair must
agree to 10 times that.  RT4 runs the thermal problems with nstokes <= 2.

The nstokes = 1 problem is unpolarized (isotropic and Henyey-Greenstein
scattering, unpolarized sources, a Lambertian surface), so its I does not
depend on the truncation and VDISORT's four-component I is the reference
for RT3's and RT4's scalar I.

Run it without ARTS_HEADLESS to plot each comparison: the solutions'
plots (pyarts3.plots.cppvdisort, RT4Result and RT3Result) on shared axes.
"""

import os

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts

FREQUENCY = 89e9
NMU = 8
MAX_DELTA_TAU = 1e-7
NATIVE_ANGLES = 4000  # the Mie habit's scattering-angle grid
MIE_DEGREE = 64  # the degree of the Mie habit's Legendre series

hg_ext = A.ScatteringSpeciesProperty("hg", A.ParticulateProperty.Extinction)
hg_ssa = A.ScatteringSpeciesProperty("hg", A.ParticulateProperty.SingleScatteringAlbedo)
drops = A.ScatteringSpeciesProperty("drops", A.ParticulateProperty.NumberDensity)


def rayleigh():
    return A.GasScatterer(A.ConstantGasScattering(2e-31), A.RayleighGasScattering(0.0))


def isotropic():
    return A.GasScatterer(A.ConstantGasScattering(2e-31), A.IsotropicGasScattering())


def henyey_greenstein():
    # RT3 keeps Legendre degrees up to 2 nmu - 3 = 13 (rt3.max_legendre_degree); with g = 0.2 the
    # terms it drops of p and p cos(Theta) are below 1e-8 of p, so all solvers see the same matrix
    return A.HenyeyGreensteinScatterer(hg_ext, hg_ssa, 0.2)


def mie_drops():
    """Liquid drops as a Legendre series, converted from gridded Mie data with the conversion's report.

    All solvers take the same series, so they solve the same problem whatever its accuracy; the report
    says how well it represents the gridded data, and its normalisation must suit the 1e-6 tolerance of
    the solvers' phase-function check.
    """
    x, _ = A.math.leggauss(NATIVE_ANGLES)
    angles = np.degrees(np.arccos(np.asarray(x)[::-1]))
    t_grid, f_grid = [220.0, 250.0, 280.0, 310.0], [FREQUENCY]
    habit = A.ParticleHabit.liquid_sphere(t_grid, f_grid, [1.5e-3], A.IrregularZenithAngleGrid(angles))
    habit, reports = habit.to_tro_spectral_with_report(t_grid, f_grid, MIE_DEGREE)
    (report,) = reports
    assert np.abs(np.asarray(report.normalisation_error)).max() < 1e-7, np.asarray(report.normalisation_error)
    assert np.asarray(report.tail).max() < 1e-12, "the series must have converged at its degree"
    return A.ScatteringHabit(habit, A.MonodispersePSD(drops, 0.0, 400.0), 1.0, 3.0)


def path(species):
    """Six levels from 6 km down, an exponential atmosphere and gas absorption, and the species."""
    ray_path = A.ArrayOfPropagationPathPoint()
    atm_path = A.ArrayOfAtmPoint()
    propmat = A.ArrayOfPropmatVector()
    for z in [6000.0, 4500.0, 3200.0, 2000.0, 1000.0, 0.0]:
        pp = A.PropagationPathPoint()
        pp.pos = [z, 0.0, 0.0]
        pp.los = [180.0, 0.0]
        ray_path.append(pp)
        atm = A.AtmPoint()
        atm.pressure = 1e5 * np.exp(-z / 8000.0)
        atm.temperature = 290.0 - 0.01 * z
        atm[hg_ext] = 2e-4 * np.exp(-(((z - 3000.0) / 1000.0) ** 2))
        atm[hg_ssa] = 0.9
        atm[drops] = 30.0 * np.exp(-(((z - 2500.0) / 1500.0) ** 2))
        atm_path.append(atm)
        propmat.append(A.PropmatVector(np.array([[5e-5 * np.exp(-z / 2000.0), 0, 0, 0, 0, 0, 0]])))
    s = A.ArrayOfScatteringSpecies()
    for x in species:
        s.add(x)
    return ray_path, atm_path, propmat, s


def solve(species, surface, nstokes, beam=None, phi=(0.0,)):
    """Run every solver that can represent the problem.

    Returns the solutions by solver, their [level, azimuth, stream, stokes]
    radiances by solver, and mu0.
    """
    ray_path, atm_path, propmat, species = path(species)
    freq_grid = A.AscendingGrid([FREQUENCY])
    kind, value = surface
    surf = {
        s: getattr(getattr(A, m), "LambertianSurface" if kind == "lambertian" else "FresnelSurface")(value)
        for s, m in (("vdisort", "vdisort"), ("rt3", "polradtran"), ("rt4", "polradtran"))
    }
    flux, mu0 = beam if beam else (0.0, 1.0)
    aziorder = 7 if beam else 0
    phi_rad = np.radians(phi)
    args = (ray_path, atm_path, propmat, freq_grid, 0, species)

    v = A.vdisort.main_data_from_path(
        *args,
        A.vdisort.PathSettings(
            nquad=2 * NMU,
            nfourier=aziorder + 1,
            normalisation_tolerance=1e-6,
            beam_flux=flux,
            beam_mu=mu0,
        ),
        surf["vdisort"],
        288.0,
        2.7,
    )
    u = np.asarray(v.u(tau=np.concatenate(([0.0], np.asarray(v.tau))), phi=phi_rad))
    solutions = {"VDISORT": v}
    out = {"VDISORT": (u[:, :, :NMU, :nstokes], u[:, :, NMU:, :nstokes])}

    p3 = A.rt3.problem_from_path(
        *args,
        A.rt3.PathSettings(
            nstokes=nstokes,
            nmu=NMU,
            quad=A.polradtran.QuadratureType.double_gauss,
            aziorder=aziorder,
            max_delta_tau=MAX_DELTA_TAU,
            normalisation_tolerance=1e-6,
        ),
        surf["rt3"],
        288.0,
        2.7,
    )
    if beam:
        p3.direct_flux, p3.direct_mu = flux, mu0
    r3 = A.rt3.solve(p3)
    solutions["RT3"] = r3
    out["RT3"] = tuple(np.asarray(A.rt3.azimuth_radiance(c, phi_rad)) for c in (r3.up, r3.down))

    if beam is None and nstokes <= 2:
        p4 = A.rt4.problem_from_path(
            *args,
            A.rt4.PathSettings(nstokes=nstokes, nmu=NMU, max_delta_tau=MAX_DELTA_TAU),
            surf["rt4"],
            288.0,
            2.7,
        )
        r4 = A.rt4.solve(p4)
        solutions["RT4"] = r4
        out["RT4"] = tuple(np.repeat(np.asarray(c)[:, None], len(phi), axis=1) for c in (r4.up, r4.down))

    assert np.allclose(np.asarray(r3.mu), np.asarray(v.mu)[:NMU], rtol=0, atol=1e-15), "the solvers' streams differ"
    return solutions, out, mu0


def deviation(radiances, a, b):
    """max |a - b| per Stokes component over levels, azimuths, streams and directions, / max |I_b|."""
    (ua, da), (ub, db) = radiances[a], radiances[b]
    scale = max(np.abs(ub[..., 0]).max(), np.abs(db[..., 0]).max())
    return np.maximum(np.abs(ua - ub).max(axis=(0, 1, 2)), np.abs(da - db).max(axis=(0, 1, 2))) / scale


def polarization(radiances):
    """max |Q|, |U|, |V| / max |I| of VDISORT's solution."""
    up, down = radiances["VDISORT"]
    top = np.maximum(np.abs(up).max(axis=(0, 1, 2)), np.abs(down).max(axis=(0, 1, 2)))
    return top[1:] / top[0]


# name, species, surface, nstokes, beam (flux on the horizontal, mu0), azimuths [deg],
# the solvers it must reach, and the minimum max |Q|, |U|, |V| / max |I| that shows it exercises polarization
CASES = [
    (
        "thermal, Rayleigh + Mie drops, Fresnel 3+0.2i, nstokes 2",
        [rayleigh, mie_drops],
        ("fresnel", 3.0 + 0.2j),
        2,
        None,
        (0.0,),
        ["VDISORT", "RT3", "RT4"],
        (1e-3,),
    ),
    (
        "thermal, Henyey-Greenstein, Fresnel 3+0.2i, nstokes 2",
        [henyey_greenstein],
        ("fresnel", 3.0 + 0.2j),
        2,
        None,
        (0.0,),
        ["VDISORT", "RT3", "RT4"],
        (1e-3,),
    ),
    (
        "thermal, Rayleigh + Mie drops, Lambertian 0.3, nstokes 2",
        [rayleigh, mie_drops],
        ("lambertian", 0.3),
        2,
        None,
        (0.0,),
        ["VDISORT", "RT3", "RT4"],
        (1e-3,),
    ),
    (
        "thermal, unpolarized: isotropic + Henyey-Greenstein, Lambertian 0.3, nstokes 1",
        [isotropic, henyey_greenstein],
        ("lambertian", 0.3),
        1,
        None,
        (0.0,),
        ["VDISORT", "RT3", "RT4"],
        (),
    ),
    (
        "solar (mu0 0.6) + thermal, Rayleigh + Henyey-Greenstein + Mie drops, Lambertian 0.3, nstokes 4, "
        "8 Fourier modes",
        [rayleigh, henyey_greenstein, mie_drops],
        ("lambertian", 0.3),
        4,
        (1e-15, 0.6),
        (0.0, 30.0, 75.0, 135.0, 180.0, 250.0),
        ["VDISORT", "RT3"],
        (1e-3, 1e-3, 1e-6),
    ),
]

for name, species, surface, nstokes, beam, phi, solvers, minimum in CASES:
    solutions, radiances, mu0 = solve([s() for s in species], surface, nstokes, beam, phi)
    tolerance = 10 * MAX_DELTA_TAU / mu0
    pol = polarization(radiances)
    print(f"{name}: max |Q|, |U|, |V| / max |I| = {np.round(pol[: nstokes - 1], 6)}")
    assert list(radiances) == solvers, f"{name}: ran {list(radiances)}, expected {solvers}"
    for s, least in enumerate(minimum):
        assert pol[s] > least, f"{name}: Stokes {s + 1} must be exercised (> {least} of I), got {pol[s]:.2e}"
    names = list(radiances)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            d = deviation(radiances, a, b)
            print(f"    {a:7s} vs {b:7s}  max |diff| / max |I| = {d}  tolerance {tolerance:.1e}")
            assert d.max() <= tolerance, f"{name}: {a} and {b} must agree to {tolerance:.1e} of max I, got {d.max():.3e}"

    if "ARTS_HEADLESS" not in os.environ:
        # The middle level, at the last azimuth, with each solution's own plot on shared axes
        level, azimuth, stokes = 3, phi[-1], list(range(nstokes))
        fig, ax = pyarts.plots.cppvdisort.plot(
            solutions["VDISORT"], level=level, azimuth=azimuth, stokes=stokes, label="VDISORT", color="C0", lw=2
        )
        pyarts.plots.RT3Result.plot(
            solutions["RT3"], fig=fig, ax=ax, level=level, azimuth=azimuth, stokes=stokes,
            label="RT3", color="C2", ls="", marker="x",
        )
        if "RT4" in solutions:
            pyarts.plots.RT4Result.plot(
                solutions["RT4"], fig=fig, ax=ax, level=level, stokes=stokes,
                label="RT4", color="C1", ls="", marker="o", mfc="none",
            )
        fig.suptitle(f"{name}\nlevel {level}, azimuth {azimuth:g} deg")

if "ARTS_HEADLESS" not in os.environ:
    import matplotlib.pyplot as plt

    plt.show()
