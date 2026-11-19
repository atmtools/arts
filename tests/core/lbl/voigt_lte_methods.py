"""The Voigt LTE methods must agree with line cutoffs, line targets and Zeeman lines.

The references are independent of the method under test:

1. A single line with a cutoff is zero outside its cutoff range, and inside it
   the profile without cutoff less its value at the upper cutoff frequency, as
   tabulated from independent code.
2. The three methods implement the same model, so they give the same numbers
   for a band of many lines, whatever the order of the lines in the band.
3. A derivative is the limit of a perturbation of the input.
"""

import numpy as np
import pyarts3 as pyarts

METHODS = (
    "spectral_propmatAddLines",
    "spectral_propmatAddVoigtLTE",
    "spectral_propmatMemoryIntensiveAddVoigtLTE",
)

TEMPERATURE = 280.0


def require(ok, message):
    assert ok, message


def workspace(freq_grid, fmin, fmax, pressure, zeeman=None):
    ws = pyarts.Workspace()
    ws.abs_speciesSet(species=["O2-66"])
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=fmin, fmax=fmax)
    ws.freq_grid = freq_grid
    ws.atm_point.temperature = TEMPERATURE
    ws.atm_point.pressure = pressure
    ws.atm_point["O2"] = 0.21
    ws.atm_point["O2-66"] = 0.995
    ws.select_species = "O2"
    ws.ray_point = pyarts.arts.PropagationPathPoint()
    ws.ray_point.los = [140.0, 30.0]
    if zeeman is not None:
        ws.abs_bandsSetZeeman(species="O2-66", fmin=zeeman[0], fmax=zeeman[1])
        ws.WignerInit()
        ws.atm_point.mag = [30e-6, 20e-6, 10e-6]
    ws.jac_targetsInit()
    return ws


def main_band(ws):
    return max(ws.abs_bands, key=lambda band: len(ws.abs_bands[band].lines))


def propmat(ws, method):
    ws.spectral_propmatInit()
    getattr(ws, method)(no_negative_absorption=False)
    return np.array(ws.spectral_propmat.K)


def relative_error(value, expected):
    return np.linalg.norm(value - expected) / max(np.linalg.norm(expected), 1e-300)


def set_shape_parameter(line, parameter, value):
    for species in ("O2", "AIR"):
        line.ls.single_models[pyarts.arts.SpeciesEnum(species)][parameter] = (
            pyarts.arts.TemperatureModel("T0", [value])
        )


# %% 1. One symmetric, pressure-shifted line with a cutoff

CUTOFF = 1e9
PRESSURE = 1e4
SHIFT = 2e3  # Hz/Pa, the same for all broadening species
OFFSETS = np.array([-3, -1.5, -0.9, -0.5, -0.02, 0, 0.02, 0.5, 0.9, 1.5, 3]) * 1e9

# The absorption [1/m] of the line with cutoff at the offsets inside the cutoff range
CUT = np.array([
    2.2212605357095746e-06,
    2.637707049357196e-05,
    0.0003111860824983051,
    0.0003153470174290892,
    0.0003113946944913489,
    2.682271223719598e-05,
    2.2892683792862827e-06,
])


def single_line_workspace(freq_grid, cutoff):
    ws = workspace(freq_grid, 118.74e9, 118.76e9, PRESSURE)
    line = ws.abs_bands[main_band(ws)].lines[0]
    set_shape_parameter(line, "D0", SHIFT)
    for parameter in ("Y", "G", "DV"):
        set_shape_parameter(line, parameter, 0.0)
    if cutoff:
        ws.abs_bandsSetCutoff(option="ByLine", value=CUTOFF)
    return ws


ws = single_line_workspace([1.0], False)
center = ws.abs_bands[main_band(ws)].lines[0].f0 + SHIFT * PRESSURE
grid = center + OFFSETS
inside = np.abs(OFFSETS) <= CUTOFF

for method in METHODS:
    cut = propmat(single_line_workspace(grid, True), method)[:, 0]

    require(
        np.all(cut[~inside] == 0.0),
        f"{method}: a line must not contribute outside its cutoff range, got {cut[~inside]}",
    )
    require(
        np.allclose(cut[inside], CUT, rtol=1e-12, atol=0),
        f"{method}: inside the cutoff range got {cut[inside]}, expected {CUT}",
    )

# %% 2. Many lines with a cutoff, also when the band is not sorted by frequency

BAND_GRID = np.linspace(50e9, 70e9, 41)
BAND_CUTOFF = 5e9


def band_workspace(freq_grid=BAND_GRID, rotate=False, pressure=1e4, zeeman=None, cutoff=True):
    ws = workspace(freq_grid, 40e9, 130e9, pressure, zeeman)
    if cutoff:
        ws.abs_bandsSetCutoff(option="ByLine", value=BAND_CUTOFF)
    if rotate:
        for band in ws.abs_bands:
            lines = list(ws.abs_bands[band].lines)
            half = len(lines) // 2
            ws.abs_bands[band].lines = lines[half:] + lines[:half]
    return ws


reference = propmat(band_workspace(), METHODS[-1])
for rotate in (False, True):
    for method in METHODS:
        value = propmat(band_workspace(rotate=rotate), method)
        require(
            np.allclose(value, reference, rtol=1e-9, atol=1e-12 * np.abs(reference).max()),
            f"{method} with a line cutoff and {'unsorted' if rotate else 'sorted'} lines differs from "
            f"{METHODS[-1]} by up to {np.abs(value - reference).max():.3e} (largest value "
            f"{np.abs(reference).max():.3e})",
        )

# %% 3. Line targets of a line in the middle of a band with a cutoff

LINE = 20


def einstein_target(line):
    return lambda ws: ws.jac_targetsAddLineParameter(band=main_band(ws), line=line, parameter="a")


def shape_target(parameter):
    return lambda ws: ws.jac_targetsAddLineShapeParameter(
        band=main_band(ws), line=LINE, species="O2", parameter=parameter, coefficient="X0"
    )


def check_derivative(name, method, make_workspace, add_target, tolerance):
    ws = make_workspace()
    add_target(ws)
    propmat(ws, method)
    analytic = np.array(ws.spectral_propmat.dK)[0]

    key = ws.jac_targets.line[0].type
    original = key.get_value(ws.abs_bands)
    step = 1e-4 * abs(original)
    key.set_value(ws.abs_bands, original + step)
    plus = propmat(ws, method)
    key.set_value(ws.abs_bands, original - step)
    minus = propmat(ws, method)
    expected = (plus - minus) / (2 * step)

    error = relative_error(analytic, expected)
    require(
        error < tolerance,
        f"{name} using {method}: relative perturbation error {error:.3e} exceeds {tolerance:.3e}; "
        f"analytic norm {np.linalg.norm(analytic):.3e}, expected norm {np.linalg.norm(expected):.3e}",
    )


# The absorption is linear in the Einstein coefficient and in the line mixing,
# so only the line width is limited by the accuracy of the Faddeeva derivative.
for method in METHODS:
    check_derivative("Einstein A with cutoff", method, band_workspace, einstein_target(LINE), 1e-7)
    check_derivative("line shape Y with cutoff", method, band_workspace, shape_target("Y"), 1e-7)
    check_derivative("line shape G0 with cutoff", method, band_workspace, shape_target("G0"), 5e-3)

# %% 4. Line targets in a band where only some lines have Zeeman effect

ZEEMAN = (118e9, 119e9)
ws = band_workspace(zeeman=ZEEMAN, cutoff=False)
zeeman_line = [i for i, line in enumerate(ws.abs_bands[main_band(ws)].lines) if line.z.on]
require(len(zeeman_line) == 1 and zeeman_line[0] != LINE, f"unexpected Zeeman lines: {zeeman_line}")
zeeman_line = zeeman_line[0]

for name, line, width in (("plain", LINE, 20e6), ("Zeeman", zeeman_line, 3e6)):
    f0 = ws.abs_bands[main_band(ws)].lines[line].f0
    freq_grid = f0 + np.linspace(-width, width, 21)
    for method in METHODS:
        check_derivative(
            f"Einstein A of a {name} line in a partly Zeeman band",
            method,
            lambda: band_workspace(freq_grid, pressure=1e2, zeeman=ZEEMAN, cutoff=False),
            einstein_target(line),
            1e-7,
        )
