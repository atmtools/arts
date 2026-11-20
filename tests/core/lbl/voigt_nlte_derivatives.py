"""Finite-difference tests for every supported non-LTE Voigt derivative."""

import os

import numpy as np
import pyarts3 as pyarts


LEVELS = pyarts.arts.ArrayOfQuantumLevelIdentifier(
    [
        "H2O-161 J 1 Ka 0 Kc 1",
        "H2O-161 J 1 Ka 1 Kc 0",
        "H2O-161 J 2 Ka 1 Kc 2",
        "H2O-161 J 2 Ka 2 Kc 1",
        "H2O-161 J 3 Ka 0 Kc 3",
        "H2O-161 J 3 Ka 1 Kc 2",
        "H2O-161 J 3 Ka 2 Kc 1",
    ]
)

O2_UPPER = pyarts.arts.QuantumLevelIdentifier(
    "O2-66 ElecStateLabel X Lambda 0 S 1 v 0 N 1 J 1"
)
O2_LOWER = pyarts.arts.QuantumLevelIdentifier(
    "O2-66 ElecStateLabel X Lambda 0 S 1 v 0 N 1 J 0"
)
O2_BAND_LEVEL = pyarts.arts.QuantumLevelIdentifier(
    "O2-66 ElecStateLabel X Lambda 0 S 1 v 0"
)
PLOT = "ARTS_HEADLESS" not in os.environ

PLOT_CASES = {
    "temperature",
    "pressure with line cutoff",
    "magnetic field u",
    f"population {LEVELS[0]}",
    "line f0",
    "line shape G0 X0",
}
PLOT_FIGURES = []


def plot_quantity(axes, frequency, analytic, expected, components, title):
    active = [
        i
        for i in range(analytic.shape[1])
        if np.any(analytic[:, i] != 0.0) or np.any(expected[:, i] != 0.0)
    ]
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
    axes[0].set_title(title)
    axes[0].set_ylabel("Derivative")
    axes[0].legend(ncols=2)
    axes[1].set_xlabel("Frequency [GHz]")
    axes[1].set_ylabel("Analytic - finite difference")
    axes[1].legend()


def plot_derivative(name, frequency, analytic, expected):
    if not PLOT or name not in PLOT_CASES:
        return

    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(14, 8), sharex="col")
    plot_quantity(
        axes[:, 0],
        frequency,
        analytic[0],
        expected[0],
        ("A", "B", "C", "D", "U", "V", "W"),
        "Propagation matrix",
    )
    plot_quantity(
        axes[:, 1],
        frequency,
        analytic[1],
        expected[1],
        ("I", "Q", "U", "V"),
        "NLTE source vector",
    )
    figure.suptitle(name)
    figure.tight_layout()
    PLOT_FIGURES.append(figure)


def workspace(zeeman=False, cutoff=False):
    ws = pyarts.Workspace()
    ws.abs_speciesSet(species=["H2O"])
    ws.abs_bands.readxml("../nlte/nlte_lines.xml")
    ws.abs_bandsSetNonLTE()

    ws.atm_point.temperature = 150.0
    ws.atm_point.pressure = 1000.0
    # The other broadening species in the catalogue are deliberately absent:
    # atmospheric species are optional and must behave as zero VMR.
    ws.atm_point["H2O"] = 1.0
    ws.atm_point["H2O-161"] = 0.997
    for i, level in enumerate(LEVELS):
        ws.atm_point[level] = 0.20 + 0.01 * i

    ws.freq_grid = np.linspace(556.92e9, 556.952e9, 41)
    ws.select_species = "H2O"
    ws.ray_point = pyarts.arts.PropagationPathPoint()
    ws.ray_point.los = [180.0, 0.0]
    if cutoff:
        ws.abs_bandsSetCutoff(option="ByLine", value=10e6)
    if zeeman:
        ws.abs_bandsSetZeeman(
            species="H2O-161", fmin=556.92e9, fmax=556.952e9
        )
        ws.WignerInit()
        ws.atm_point.mag = [30e-6, 20e-6, 10e-6]
    return ws


def magnetic_workspace(zeeman=False, cutoff=False):
    ws = pyarts.Workspace()
    ws.abs_speciesSet(species=["O2-66"])
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=118.74e9, fmax=118.76e9)
    ws.abs_bandsSetNonLTE()
    ws.abs_bandsSetZeeman(species="O2-66", fmin=118.74e9, fmax=118.76e9)
    ws.WignerInit()
    ws.freq_grid = np.linspace(118.74e9, 118.76e9, 41)
    ws.atm_point.temperature = 280.0
    ws.atm_point.pressure = 1e5
    ws.atm_point["O2"] = 0.21
    ws.atm_point["O2-66"] = 0.995
    ws.atm_point[O2_UPPER] = 0.2
    ws.atm_point[O2_LOWER] = 0.3
    ws.atm_point[O2_BAND_LEVEL] = 0.25
    ws.atm_point.mag = [30e-6, 20e-6, 10e-6]
    ws.select_species = "O2"
    ws.ray_point = pyarts.arts.PropagationPathPoint()
    ws.ray_point.los = [160.0, 30.0]
    return ws


def calculate(ws):
    ws.spectral_propmatInit()
    ws.spectral_propmatAddLines(no_negative_absorption=False)
    return (
        np.array(ws.spectral_propmat.K),
        np.array(ws.spectral_propmat.nlte),
    )


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


def absent_species_value(species):
    """An absent species in AtmPoint has zero VMR"""
    return lambda ws: (
        lambda ws: 0.0,
        lambda ws, value: ws.atm_point.__setitem__(species, value),
    )


def check(
    name,
    add_target,
    step,
    rtol=1e-3,
    atol=1e-14,
    zeeman=False,
    cutoff=False,
    factory=workspace,
    expect_nonzero=True,
    value=target_value,
):
    ws = factory(zeeman=zeeman, cutoff=cutoff)
    ws.jac_targetsInit()
    add_target(ws)
    get_value, set_value = value(ws)
    calculate(ws)
    analytic = (
        np.array(ws.spectral_propmat.dK)[0],
        np.array(ws.spectral_propmat.dnlte)[0],
    )

    original = get_value(ws)
    set_value(ws, original + step)
    plus = calculate(ws)
    set_value(ws, original - step)
    minus = calculate(ws)
    set_value(ws, original)

    finite_difference_is_nonzero = False
    expected_values = []
    for quantity, actual, upper, lower in zip(
        ("propagation matrix", "source vector"), analytic, plus, minus
    ):
        expected = (upper - lower) / (2.0 * step)
        expected_values.append(expected)
        finite_difference_is_nonzero |= np.any(expected != 0.0)
        error = np.linalg.norm(actual - expected)
        limit = atol + rtol * np.linalg.norm(expected)
        assert error <= limit, (
            f"{name}: {quantity}: relative finite-difference error "
            f"{error / max(np.linalg.norm(expected), atol):.3e} exceeds {rtol:.3e}; "
            f"analytic norm {np.linalg.norm(actual):.3e}, expected norm "
            f"{np.linalg.norm(expected):.3e}; center "
            f"{actual[len(actual) // 2, 0]:.3e} vs "
            f"{expected[len(expected) // 2, 0]:.3e}"
        )
    if expect_nonzero:
        assert finite_difference_is_nonzero, f"{name}: finite difference is zero"
    plot_derivative(
        name,
        np.asarray(ws.freq_grid),
        analytic,
        tuple(expected_values),
    )


check("temperature", lambda ws: ws.jac_targetsAddTemperature(), 1e-2)
check("pressure", lambda ws: ws.jac_targetsAddPressure(), 1e-1, 3e-2)
check(
    "pressure with line cutoff",
    lambda ws: ws.jac_targetsAddPressure(),
    1e-1,
    3e-2,
    cutoff=True,
)
check(
    "absorber VMR",
    lambda ws: ws.jac_targetsAddSpeciesVMR(species="H2O"),
    1e-5,
)
check(
    "unrelated optional isotopologue ratio",
    lambda ws: ws.jac_targetsAddSpeciesIsotopologueRatio(species="H2O-181"),
    1e-5,
    expect_nonzero=False,
)

for component in ("u", "v", "w"):
    check(
        f"magnetic field {component}",
        lambda ws, component=component: ws.jac_targetsAddMagneticField(
            component=component
        ),
        1e-8,
        rtol=1e-1,
        zeeman=True,
        factory=magnetic_workspace,
    )
check(
    "optional broadener VMR",
    lambda ws: ws.jac_targetsAddSpeciesVMR(species="CO2"),
    1e-6,
    value=absent_species_value("CO2"),
)
check(
    "frequency/wind",
    lambda ws: ws.jac_targetsAddWindField(component="u"),
    10.0,
    5e-3,
    value=frequency_value,
)
check(
    "isotopologue ratio",
    lambda ws: ws.jac_targetsAddSpeciesIsotopologueRatio(species="H2O-161"),
    1e-5,
)

for level in LEVELS[:2]:
    check(
        f"population {level}",
        lambda ws, level=level: ws.jac_targetsAddAtmosphere(target=level),
        1e-6,
    )


def selected_band(ws):
    return min(
        ws.abs_bands,
        key=lambda key: abs(ws.abs_bands[key].lines[0].f0 - 556.936e9),
    )


for parameter, step, tolerance in (
    ("f0", 10.0, 2e-4),
    ("e0", 1e-24, 5e-5),
    ("a", 1e-7, 5e-5),
):
    check(
        f"line {parameter}",
        lambda ws, parameter=parameter: ws.jac_targetsAddLineParameter(
            band=selected_band(ws), line=0, parameter=parameter
        ),
        step,
        tolerance,
        expect_nonzero=parameter != "e0",
    )


for parameter, step in (("G0", 1.0), ("D0", 1.0)):
    check(
        f"line shape {parameter} X0",
        lambda ws, parameter=parameter: ws.jac_targetsAddLineShapeParameter(
            band=selected_band(ws),
            line=0,
            species="H2O",
            parameter=parameter,
            coefficient="X0",
        ),
        step,
    )

if PLOT_FIGURES:
    import matplotlib.pyplot as plt

    plt.show()
