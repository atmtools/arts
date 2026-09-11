"""Compress a wide water-line spectrum to one wind-sensitive measurement.

A clear-sky LTE calculation over a 1 GHz band also sees the broad line at a
30 km tangent altitude. At 60 km,
most of that band contributes almost no wind information. ReducedOEM retains
a signed combination of channels, not the radiance of one central channel.
Prints repeated OEM-call timings and separate basis-preparation timings.
Run normally for plots; CTest sets ARTS_HEADLESS to suppress them.
"""

import os
from time import perf_counter

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
PRIOR_WIND = 10.0  # m/s, uniform eastward wind; the only retrieved parameter
TRUE_WIND = 30.0
PRIOR_STD = 80.0
NOISE_STD = 0.1  # K per independent channel; observations below are noise-free
OFFSET = np.linspace(-500e6, 500e6, 2001)
TIMING_REPEATS = 6


def run():
    ws = pyarts.Workspace()
    ws.abs_speciesSet(species=["H2O-161"])
    ws.ReadCatalogData()
    # Isolate the strong 557 GHz transition, including its pressure-broadened
    # wings throughout the much wider measurement band.
    ws.abs_bandsSelectFrequencyByLine(fmin=556.9e9, fmax=557e9)
    lines = [line for band in ws.abs_bands.values() for line in band.lines]
    assert len(lines) == 1
    line_center = float(lines[0].f0)
    ws.freq_grid = line_center + OFFSET
    ws.spectral_propmat_agendaAuto()
    ws.surf_fieldEarth(model="Sphere")
    ws.surf_field["t"] = 295.0
    ws.atm_fieldRead(
        toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
    )
    # The existing wind derivative has a known all-zero-wind limitation (see
    # core/wind/freq_grid_derivatives.py). Use a nonzero reference and verify J.
    ws.atm_field["wind_u"] = PRIOR_WIND
    ws.atm_field["wind_v"] = 0.0
    ws.atm_field["wind_w"] = 0.0
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.max_stepsize = 10e3
    ws.ray_path_observer_agendaSetGeometric(add_limb=1)
    radius = float(ws.surf_field.ellipsoid[0])
    pos = [600e3, 0.0, 0.0]

    def los(tangent_altitude):
        zenith = 180 - np.rad2deg(np.arcsin(
            (radius + tangent_altitude) / (radius + pos[0])))
        return [zenith, 90.0]  # East-looking limb: sensitive to eastward wind.

    ws.ray_pathGeometric(pos=pos, los=los(60e3), add_limb=1)
    np.testing.assert_allclose(min(p.pos[0] for p in ws.ray_path), 60e3, atol=1e-3, rtol=0)
    ws.measurement_sensorSimple(pos=pos, los=los(60e3))
    ws.RetrievalInit()
    ws.RetrievalAddWindField(component="u", matrix=[[PRIOR_STD**2]])
    ws.RetrievalFinalizeDiagonal()
    ws.model_state_vec_aprioriFromData()
    ws.measurement_vecFromSensor()
    prior_spectrum = np.array(ws.measurement_vec)
    jacobian = np.array(ws.measurement_jac)
    assert jacobian.shape == (len(OFFSET), 1)
    j = jacobian[:, 0]

    # Independent central difference through the real radiative-transfer path.
    perturbed = []
    for delta in (-0.1, 0.1):
        ws.atm_field["wind_u"] = PRIOR_WIND + delta
        ws.measurement_vecFromSensor()
        perturbed.append(np.array(ws.measurement_vec))
    numerical_jacobian = (perturbed[1] - perturbed[0]) / 0.2
    derivative_error = np.linalg.norm(j - numerical_jacobian) / np.linalg.norm(j)
    assert derivative_error < 1e-3, derivative_error

    ws.atm_field["wind_u"] = TRUE_WIND
    ws.measurement_vecFromSensor()
    observations = np.array(ws.measurement_vec)

    # Demonstrate why an instrument also viewing deeper layers has a wide band.
    ws.atm_field["wind_u"] = PRIOR_WIND
    ws.measurement_sensorSimple(pos=pos, los=los(30e3))
    ws.measurement_vecFromSensor()
    deeper_spectrum = np.array(ws.measurement_vec)
    ws.measurement_sensorSimple(pos=pos, los=los(60e3))
    ws.measurement_vec = observations
    ws.measurement_jac = jacobian
    ws.measurement_vec_error_covmatConstant(value=NOISE_STD**2)

    # Only the measurement dimension shrinks: the full state already has size 1.
    # B rescales that single wind parameter; C compresses 2001 channels to one.
    started = perf_counter()
    ws.ReducedOEMBasisCalc()
    basis_calc_seconds = perf_counter() - started
    started = perf_counter()
    ws.ReducedOEMBasisReduce()
    basis_reduce_seconds = perf_counter() - started
    B = np.array(ws.model_state_basis_mat)
    C = np.array(ws.measurement_basis_mat)
    assert B.shape == (1, 1) and C.shape == (1, len(OFFSET))
    assert float(ws.oem_basis_lost_dofs) == 0
    assert float(ws.oem_basis_lost_information_bits) == 0
    np.testing.assert_allclose(B @ B.T, [[PRIOR_STD**2]], rtol=1e-12)
    np.testing.assert_allclose(NOISE_STD**2 * (C @ C.T), [[1.]], rtol=1e-12)

    # For one state and equal independent noise, C is proportional to J^T.
    # Thus (C J)^2 = J^T Se^-1 J: compression preserves the local Fisher
    # information, despite discarding 2000 orthogonal measurement directions.
    weights = C[0] * np.sign((C @ jacobian)[0, 0])  # Orient only the plotted copy.
    np.testing.assert_allclose(weights * NOISE_STD, j / np.linalg.norm(j), atol=1e-12)
    full_information = float(j @ j / NOISE_STD**2)
    compressed_information = float((C @ jacobian)[0, 0]**2)
    np.testing.assert_allclose(compressed_information, full_information, rtol=1e-12)
    information_fraction = j**2 / (j @ j)
    core = abs(OFFSET) <= 100e6
    far_wing = abs(OFFSET) >= 250e6
    assert np.mean(~core) > 0.79  # Nearly 80% of the channels are outside the core.
    assert information_fraction[core].sum() > 0.99
    assert information_fraction[far_wing].sum() < 2e-4
    assert np.max(abs(j[far_wing])) < 0.003 * np.max(abs(j))
    assert np.mean(deeper_spectrum[far_wing]) > 3 * np.mean(prior_spectrum[far_wing])
    # Exactly at line center, a shift has little first-order effect. Both
    # sides of the feature are needed, with opposite signs in the combination.
    assert abs(j[len(j)//2]) < 0.02 * np.max(abs(j))
    assert weights.min() < 0 < weights.max()

    # This small Doppler shift needs little damping. Starting at 1 also
    # permits the convergence check once the small-step criterion is met.
    lm_settings = arts.LevenbergMarquardtSettings(initial_damping=1.0)

    def retrieve(name):
        # Every call starts from the same prior, with no cached fit/Jacobian.
        # Resetting inputs and checking/copying outputs are outside the timer.
        ws.atm_field["wind_u"] = PRIOR_WIND
        ws.model_state_vec = []
        ws.measurement_vec_fit = []
        ws.measurement_jac = arts.Matrix()
        method = getattr(ws, name)
        started = perf_counter()
        method(
            method="lm", max_iter=30, stop_dx=0.01,
            lm_ga_settings=lm_settings,
        )
        elapsed = perf_counter() - started
        assert ws.oem_diagnostics.status == arts.OptimalEstimationStatus.Converged, ws.oem_diagnostics
        assert not len(ws.oem_diagnostics.errors), ws.oem_diagnostics
        wind = float(ws.model_state_vec[0])
        fitted = np.array(ws.measurement_vec_fit)
        assert fitted.shape == observations.shape  # Output remains the full spectrum.
        assert np.asarray(ws.measurement_jac).shape == jacobian.shape
        np.testing.assert_allclose(wind, TRUE_WIND, atol=0.01, rtol=0)
        return (wind, fitted), elapsed, int(ws.oem_diagnostics.iterations)

    names = ("OEM", "ReducedOEM")
    fits = {}
    # Exclude one warm-up per method, then balance which method runs first.
    for name in names:
        fits[name], _, _ = retrieve(name)
    timings = {name: [] for name in names}
    iterations = {name: [] for name in names}
    for repeat in range(TIMING_REPEATS):
        for name in names if repeat % 2 == 0 else names[::-1]:
            (wind, fitted), elapsed, count = retrieve(name)
            np.testing.assert_allclose(wind, fits[name][0], atol=1e-8, rtol=0)
            np.testing.assert_allclose(fitted, fits[name][1], atol=1e-9, rtol=0)
            timings[name].append(elapsed)
            iterations[name].append(count)
    # A fixed basis is locally lossless, not a global nonlinear guarantee.
    # Here the small physical shift gives nearly identical nonlinear LM fits.
    np.testing.assert_allclose(fits["OEM"][0], fits["ReducedOEM"][0], atol=1e-4, rtol=0)
    np.testing.assert_allclose(fits["OEM"][1], fits["ReducedOEM"][1], atol=1e-5, rtol=0)
    print(f"H2O {line_center/1e9:.6f} GHz, tangent 60 km, band +/-500 MHz")
    print(f"Measurement dimension: {len(OFFSET)} -> {C.shape[0]}; state dimension: 1 -> 1")
    print(f"Wind Fisher information within +/-100 MHz: {100*information_fraction[core].sum():.4f}%")
    print(f"Wind Fisher information beyond +/-250 MHz: {100*information_fraction[far_wing].sum():.6f}%")
    print(f"Jacobian central-difference relative error: {derivative_error:.3g}")
    print(f"Wind: truth={TRUE_WIND:g}, prior={PRIOR_WIND:g}, "
          f"OEM={fits['OEM'][0]:.6f}, ReducedOEM={fits['ReducedOEM'][0]:.6f} m/s")
    print(f"\nTiming: {TIMING_REPEATS} runs per method after one warm-up each; alternating order.")
    print("Same prior and empty fit/Jacobian each run; full forward calculations and gain included.")
    print("Basis preparation, input resets, validation of results and plots are outside the OEM-call timer.")
    medians = {name: float(np.median(timings[name])) for name in names}
    for name in names:
        print(f"  {name}: median {1e3*medians[name]:.3f} ms, "
              f"range {1e3*min(timings[name]):.3f}-{1e3*max(timings[name]):.3f} ms, "
              f"iterations {min(iterations[name])}-{max(iterations[name])}")
    delta = medians["ReducedOEM"] - medians["OEM"]
    print(f"  ReducedOEM - OEM: {1e3*delta:+.3f} ms ({100*delta/medians['OEM']:+.2f}%; positive is slower)")
    print(f"  Basis preparation from existing Jacobian (one call each): "
          f"Calc {1e3*basis_calc_seconds:.3f} ms, Reduce {1e3*basis_reduce_seconds:.3f} ms")
    print(f"  Basis preparation + median ReducedOEM: "
          f"{1e3*(basis_calc_seconds + basis_reduce_seconds + medians['ReducedOEM']):.3f} ms")
    # Timings are observations, not pass/fail limits: scheduling and hardware vary.
    return dict(line_center=line_center, prior=prior_spectrum, deeper=deeper_spectrum,
                observations=observations, jacobian=j, numerical_jacobian=numerical_jacobian,
                weights=weights, information_fraction=information_fraction, fits=fits,
                timings=timings, basis_calc_seconds=basis_calc_seconds,
                basis_reduce_seconds=basis_reduce_seconds)


def plot(result):
    import matplotlib.pyplot as plt

    f = OFFSET / 1e6
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax = axes[0, 0]
    ax.plot(f, result["deeper"], color="0.6", label="30 km, prior wind")
    ax.plot(f, result["prior"], label="60 km, prior wind")
    ax.set_ylabel("Brightness temperature [K]")
    ax.set_title("Wide band for observing deeper layers")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(f, result["jacobian"], label="ARTS wind Jacobian")
    ax.plot(f, result["numerical_jacobian"], "--", label="Central difference")
    ax.set_xscale("symlog", linthresh=5)
    ticks = [-500, -100, -10, 0, 10, 100, 500]
    ax.set_xticks(ticks, labels=[str(tick) for tick in ticks])
    ax.set_xlim(-500, 500)
    ax.set_ylabel("Wind derivative [K / (m/s)]")
    ax.set_title("Wind sensitivity at 60 km; expanded center")
    ax.legend()

    ax = axes[1, 0]
    weights = result["weights"] / np.max(abs(result["weights"]))
    ax.plot(f, weights, color="tab:green", label="Signed weights (normalized)")
    ax.axvspan(-100, 100, color="tab:green", alpha=.1, label=">99% of wind information")
    ax.axhline(0, color="0.5", linewidth=.5)
    ax.set_ylabel("Signed relative weight")
    ax.set_title(f"{len(f)} channels -> one weighted measurement")
    ax.legend(fontsize="small")

    ax = axes[1, 1]
    ax.plot(f, result["observations"] - result["prior"], color="k", label="Synthetic truth - prior")
    for name, (wind, fitted) in result["fits"].items():
        ax.plot(f, fitted - result["prior"], "--" if name == "OEM" else ":",
                label=f"{name}: {wind:.4f} m/s")
    ax.set_xlim(-60, 60)
    ax.set_ylabel("Change from prior spectrum [K]")
    ax.set_title(f"Wind retrieval: prior {PRIOR_WIND:g}, truth {TRUE_WIND:g} m/s")
    ax.legend(fontsize="small")
    for ax in axes.flat:
        ax.set_xlabel("Frequency offset from line center [MHz]")
        ax.grid(True, alpha=.3)
    fig.suptitle(f"H2O {result['line_center']/1e9:.6f} GHz: measurement reduction at a 60 km limb")
    fig.tight_layout()
    return fig


def plot_differences(result):
    """Expose residuals and numerical differences hidden by overlapping curves."""
    import matplotlib.pyplot as plt

    f = OFFSET / 1e6
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for name, (_, fitted) in result["fits"].items():
        axes[0, 0].plot(f, fitted - result["observations"],
                        linestyle="-" if name == "OEM" else "--", label=name)
    axes[0, 0].set_title("Retrieval residuals at 60 km")
    axes[0, 0].set_ylabel("Fit − measurement [K]")
    axes[0, 0].legend()

    full_wind, full_fit = result["fits"]["OEM"]
    reduced_wind, reduced_fit = result["fits"]["ReducedOEM"]
    axes[0, 1].plot(f, (reduced_fit - full_fit) * 1e6)
    axes[0, 1].set_title("Effect of measurement compression")
    axes[0, 1].set_ylabel("ReducedOEM − OEM [µK]")

    axes[1, 0].plot(f, result["jacobian"] - result["numerical_jacobian"])
    axes[1, 0].set_title("Wind Jacobian difference")
    axes[1, 0].set_ylabel("ARTS − central difference [K / (m/s)]")

    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.set_xscale("symlog", linthresh=5)
        ticks = [-500, -100, -10, 0, 10, 100, 500]
        ax.set_xticks(ticks, labels=[str(tick) for tick in ticks])
        ax.set_xlim(-500, 500)
        ax.set_xlabel("Frequency offset from line center [MHz]; expanded center")

    ax = axes[1, 1]
    ax.plot([0, 1], np.array([full_wind, reduced_wind]) - TRUE_WIND, "o")
    ax.set_xticks([0, 1], labels=["OEM", "ReducedOEM"])
    ax.set_xlim(-.5, 1.5)
    ax.set_title("Retrieved wind error")
    ax.set_ylabel("Retrieved wind − truth [m/s]")
    ax.text(.5, .5, f"ReducedOEM − OEM = {reduced_wind - full_wind:.3g} m/s",
            transform=ax.transAxes, ha="center")
    for ax in axes.flat:
        ax.axhline(0, color="k", linewidth=.5)
        ax.grid(True, alpha=.3)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False)
    fig.suptitle("ReducedOEM water-line regression: residuals and differences")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    result = run()
    if "ARTS_HEADLESS" not in os.environ:
        plot(result)
        plot_differences(result)
        import matplotlib.pyplot as plt
        plt.show()
