"""Finite-difference tests for LTE Voigt atmospheric and line derivatives."""

import os

import numpy as np
import pyarts3 as pyarts


LINE_F0 = 118750348044.712
PLOT = "ARTS_HEADLESS" not in os.environ

PLOT_CASES = {
    ("temperature", "spectral_propmatAddLines", "VP_LTE"),
    ("pressure with line cutoff", "spectral_propmatMemoryIntensiveAddVoigtLTE", "VP_LTE"),
    ("line f0", "spectral_propmatAddLines", "VP_LTE_MIRROR"),
    ("AIR line shape G0 X0", "spectral_propmatMemoryIntensiveAddVoigtLTE", "VP_LTE"),
    ("absent normalized broadener VMR", "spectral_propmatMemoryIntensiveAddVoigtLTE", "VP_LTE"),
}
PLOT_FIGURES = []


def plot_derivative(name, implementation, lineshape, frequency, analytic, expected):
    if not PLOT or (name, implementation, lineshape) not in PLOT_CASES:
        return

    import matplotlib.pyplot as plt

    components = ("A", "B", "C", "D", "U", "V", "W")
    analytic = np.asarray(analytic)
    expected = np.asarray(expected)
    active = [
        i
        for i in range(analytic.shape[1])
        if np.any(analytic[:, i] != 0.0) or np.any(expected[:, i] != 0.0)
    ]

    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i in active:
        (line,) = axes[0].plot(
            frequency / 1e9,
            analytic[:, i],
            "-",
            label=f"{components[i]} analytic",
            zorder=1,
        )
        axes[0].plot(
            frequency / 1e9,
            expected[:, i],
            "--",
            label=f"{components[i]} finite difference",
            zorder=2,
        )
        axes[1].plot(
            frequency / 1e9,
            analytic[:, i] - expected[:, i],
            color=line.get_color(),
            label=components[i],
        )
    axes[0].set_ylabel("Derivative")
    axes[0].legend(ncols=2)
    axes[1].set_xlabel("Frequency [GHz]")
    axes[1].set_ylabel("Analytic - finite difference")
    axes[1].legend()
    figure.suptitle(f"{name}: {implementation}, {lineshape}")
    figure.tight_layout()
    PLOT_FIGURES.append(figure)


def workspace(lineshape="VP_LTE", zeeman=False, cutoff=False):
    ws = pyarts.Workspace()
    ws.abs_speciesSet(species=["O2-66"])
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=118.74e9, fmax=118.76e9)
    # Cover both the Doppler core and pressure-broadened wings.  In a core-only
    # grid, number-density and Lorentz-width pressure effects nearly cancel.
    ws.freq_grid = np.linspace(-3e9, 3e9, 81) + LINE_F0
    ws.atm_point.temperature = 280.0
    ws.atm_point.pressure = 1e5
    ws.atm_point["O2"] = 0.21
    ws.atm_point["O2-66"] = 0.995
    ws.select_species = "O2"
    ws.ray_point = pyarts.arts.PropagationPathPoint()
    ws.ray_point.los = [180.0, 0.0]
    band = next(iter(ws.abs_bands))
    ws.abs_bands[band].lineshape = lineshape
    if cutoff:
        ws.abs_bandsSetCutoff(option="ByLine", value=1e9)
    species_model = ws.abs_bands[band].lines[0].ls.single_models[
        pyarts.arts.SpeciesEnum("O2")
    ]
    for parameter in ("D0", "G", "DV"):
        species_model[parameter] = pyarts.arts.TemperatureModel("T0", [0.0])
    if zeeman:
        ws.abs_bandsSetZeeman(species="O2-66", fmin=118.74e9, fmax=118.76e9)
        ws.WignerInit()
        ws.atm_point.mag = [30e-6, 20e-6, 10e-6]
    return ws


def calculate(ws, implementation):
    ws.spectral_propmatInit()
    getattr(ws, implementation)(no_negative_absorption=False)
    return np.array(ws.spectral_propmat.K)


def target_value(ws):
    """Getter and setter of the quantity of the only Jacobian target"""
    if len(ws.jac_targets.line):
        key = ws.jac_targets.line[0].type
        return (
            lambda ws: key.get_value(ws.abs_bands),
            lambda ws, value: key.set_value(ws.abs_bands, value),
        )
    key = ws.jac_targets.atm[0].type
    return (
        lambda ws: ws.atm_point[key],
        lambda ws, value: ws.atm_point.__setitem__(key, value),
    )


def frequency_value(ws):
    """The wind derivative is the frequency derivative"""
    return (
        lambda ws: np.array(ws.freq_grid),
        lambda ws, value: setattr(ws, "freq_grid", value),
    )


def check(
    name,
    implementation,
    add_target,
    step,
    lineshape="VP_LTE",
    tolerance=5e-3,
    zeeman=False,
    cutoff=False,
    value=target_value,
):
    ws = workspace(lineshape, zeeman, cutoff)
    ws.jac_targetsInit()
    add_target(ws)
    get_value, set_value = value(ws)
    calculate(ws, implementation)
    analytic = np.array(ws.spectral_propmat.dK)[0]

    original = get_value(ws)
    set_value(ws, original + step)
    plus = calculate(ws, implementation)
    set_value(ws, original - step)
    minus = calculate(ws, implementation)
    set_value(ws, original)
    expected = (plus - minus) / (2.0 * step)

    plot_derivative(
        name,
        implementation,
        lineshape,
        np.asarray(ws.freq_grid),
        analytic,
        expected,
    )

    error = np.linalg.norm(analytic - expected)
    relative_error = error / max(np.linalg.norm(expected), 1e-300)
    assert relative_error < tolerance, (
        f"{name} using {implementation}/{lineshape}: relative finite-difference "
        f"error {relative_error:.3e} exceeds {tolerance:.3e}; analytic norm "
        f"{np.linalg.norm(analytic):.3e}, expected norm {np.linalg.norm(expected):.3e}"
    )


IMPLEMENTATIONS = (
    "spectral_propmatAddLines",
    "spectral_propmatAddVoigtLTE",
    "spectral_propmatMemoryIntensiveAddVoigtLTE",
)

ATMOSPHERIC_TARGETS = (
    ("temperature", lambda ws: ws.jac_targetsAddTemperature(), 1e-2, {}),
    ("pressure", lambda ws: ws.jac_targetsAddPressure(), 1.0, {}),
    (
        "absorber VMR",
        lambda ws: ws.jac_targetsAddSpeciesVMR(species="O2"),
        1e-5,
        {},
    ),
    (
        "isotopologue ratio",
        lambda ws: ws.jac_targetsAddSpeciesIsotopologueRatio(species="O2-66"),
        1e-5,
        {},
    ),
    (
        "frequency/wind",
        lambda ws: ws.jac_targetsAddWindField(component="u"),
        10.0,
        {"value": frequency_value},
    ),
)

for implementation in IMPLEMENTATIONS:
    for name, add_target, step, options in ATMOSPHERIC_TARGETS:
        check(name, implementation, add_target, step, **options)

# The mirrored profile has its own implementation and must not regress separately.
for name, add_target, step, options in ATMOSPHERIC_TARGETS:
    check(
        name,
        "spectral_propmatAddLines",
        add_target,
        step,
        lineshape="VP_LTE_MIRROR",
        **options,
    )

# Pressure derivatives have separate cutoff paths in each implementation.
for implementation in IMPLEMENTATIONS:
    check(
        "pressure with line cutoff",
        implementation,
        *ATMOSPHERIC_TARGETS[1][1:3],
        cutoff=True,
    )
check(
    "pressure with line cutoff",
    "spectral_propmatAddLines",
    *ATMOSPHERIC_TARGETS[1][1:3],
    lineshape="VP_LTE_MIRROR",
    cutoff=True,
)


for implementation in IMPLEMENTATIONS:
    check(
        "magnetic field u",
        implementation,
        lambda ws: ws.jac_targetsAddMagneticField(component="u"),
        1e-8,
        tolerance=1e-1,
        zeeman=True,
    )

check(
    "magnetic field u",
    "spectral_propmatAddLines",
    lambda ws: ws.jac_targetsAddMagneticField(component="u"),
    1e-8,
    lineshape="VP_LTE_MIRROR",
    tolerance=1e-1,
    zeeman=True,
)


def band(ws):
    return next(iter(ws.abs_bands))


def expect_runtime_error(name, callback):
    try:
        callback()
    except RuntimeError:
        return
    raise AssertionError(f"{name}: expected RuntimeError")


invalid = workspace()
invalid.jac_targetsInit()
expect_runtime_error(
    "negative line index",
    lambda: invalid.jac_targetsAddLineParameter(
        band=band(invalid), line=-1, parameter="f0"
    ),
)
expect_runtime_error(
    "unused line parameter",
    lambda: invalid.jac_targetsAddLineParameter(
        band=band(invalid), line=0, parameter="unused"
    ),
)
expect_runtime_error(
    "negative line-shape line index",
    lambda: invalid.jac_targetsAddLineShapeParameter(
        band=band(invalid),
        line=-1,
        species="O2",
        parameter="G0",
        coefficient="X0",
    ),
)
expect_runtime_error(
    "unused line-shape parameter",
    lambda: invalid.jac_targetsAddLineShapeParameter(
        band=band(invalid),
        line=0,
        species="O2",
        parameter="unused",
        coefficient="X0",
    ),
)
expect_runtime_error(
    "unused line-shape coefficient",
    lambda: invalid.jac_targetsAddLineShapeParameter(
        band=band(invalid),
        line=0,
        species="O2",
        parameter="G0",
        coefficient="unused",
    ),
)


# Matrix scaling must use global target positions, even when line and
# atmospheric targets are interleaved.
mixed = workspace()
mixed.jac_targetsInit()
mixed.jac_targetsAddLineParameter(band=band(mixed), line=0, parameter="f0")
mixed.jac_targetsAddPressure()
mixed.jac_targetsAddWindField(component="u")
calculate(mixed, "spectral_propmatMemoryIntensiveAddVoigtLTE")
mixed_jac = np.array(mixed.spectral_propmat.dK)
for position, add_target in enumerate(
    (
        lambda ws: ws.jac_targetsAddLineParameter(
            band=band(ws), line=0, parameter="f0"
        ),
        lambda ws: ws.jac_targetsAddPressure(),
        lambda ws: ws.jac_targetsAddWindField(component="u"),
    )
):
    single = workspace()
    single.jac_targetsInit()
    add_target(single)
    calculate(single, "spectral_propmatMemoryIntensiveAddVoigtLTE")
    np.testing.assert_allclose(
        mixed_jac[position],
        single.spectral_propmat.dK[0],
        err_msg=f"mixed target position {position}",
    )

for parameter, step, tolerance in (
    ("f0", 10.0, 5e-3),
    ("e0", 1e-24, 1e-6),
    ("a", 1e-5, 1e-6),
):
    for implementation in IMPLEMENTATIONS:
        check(
            f"line {parameter}",
            implementation,
            lambda ws, parameter=parameter: ws.jac_targetsAddLineParameter(
                band=band(ws), line=0, parameter=parameter
            ),
            step,
            tolerance=tolerance,
        )
    check(
        f"line {parameter}",
        "spectral_propmatAddLines",
        lambda ws, parameter=parameter: ws.jac_targetsAddLineParameter(
            band=band(ws), line=0, parameter=parameter
        ),
        step,
        lineshape="VP_LTE_MIRROR",
        tolerance=tolerance,
    )


for parameter, coefficient, step, tolerance in (
    ("G0", 0, 1e-2, 5e-3),
    ("D0", 0, 1e-3, 3e-2),
    ("Y", 0, 1e-10, 5e-3),
    ("G", 0, 1e-12, 5e-3),
    ("DV", 0, 1e-6, 3e-2),
):
    for implementation in IMPLEMENTATIONS:
        check(
            f"line shape {parameter} X{coefficient}",
            implementation,
            lambda ws, parameter=parameter, coefficient=coefficient: ws.jac_targetsAddLineShapeParameter(
                band=band(ws),
                line=0,
                species="O2",
                parameter=parameter,
                coefficient=f"X{coefficient}",
            ),
            step,
            tolerance=tolerance,
        )
    check(
        f"line shape {parameter} X{coefficient}",
        "spectral_propmatAddLines",
        lambda ws, parameter=parameter, coefficient=coefficient: ws.jac_targetsAddLineShapeParameter(
            band=band(ws),
            line=0,
            species="O2",
            parameter=parameter,
            coefficient=f"X{coefficient}",
        ),
        step,
        lineshape="VP_LTE_MIRROR",
        tolerance=tolerance,
    )

# The broadener selector must not silently fall back to the absorbing species.
for implementation in IMPLEMENTATIONS:
    check(
        "AIR line shape G0 X0",
        implementation,
        lambda ws: ws.jac_targetsAddLineShapeParameter(
            band=band(ws),
            line=0,
            species="AIR",
            parameter="G0",
            coefficient="X0",
        ),
        1e-2,
    )

# Exercise every non-X0 coefficient through the new typed target interface.
for coefficient in (1, 2, 3):
    for implementation in IMPLEMENTATIONS:
        check(
            f"O2 line shape Y X{coefficient}",
            implementation,
            lambda ws, coefficient=coefficient: ws.jac_targetsAddLineShapeParameter(
                band=band(ws),
                line=0,
                species="O2",
                parameter="Y",
                coefficient=f"X{coefficient}",
            ),
            1e-12,
        )


def add_normalized_optional_broadener_target(ws):
    models = ws.abs_bands[band(ws)].lines[0].ls.single_models
    del models[pyarts.arts.SpeciesEnum("AIR")]
    co2 = pyarts.arts.LineShapeSpeciesModel(
        models[pyarts.arts.SpeciesEnum("O2")]
    )
    co2["G0"]["X0"] *= 1.1
    models[pyarts.arts.SpeciesEnum("CO2")] = co2
    ws.jac_targetsAddSpeciesVMR(species="CO2")


def absent_co2_value(ws):
    """CO2 is absent from AtmPoint, which means a VMR of zero"""
    return (
        lambda ws: 0.0,
        lambda ws, value: ws.atm_point.__setitem__("CO2", value),
    )


# Without AIR/Bath, broadener mixing is normalized by the sum of the provided
# VMRs.  CO2 is deliberately absent from AtmPoint at the linearization point.
for implementation in IMPLEMENTATIONS:
    check(
        "absent normalized broadener VMR",
        implementation,
        add_normalized_optional_broadener_target,
        1e-6,
        tolerance=1e-2,
        value=absent_co2_value,
    )

if PLOT_FIGURES:
    import matplotlib.pyplot as plt

    plt.show()
