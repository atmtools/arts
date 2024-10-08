import pyarts
import numpy as np


def polyfit(x, y, plot=True):
    pfit = np.polyfit(x, y, deg=3)
    return pfit


def fieldfit(x, f, plot=False):
    fpfit = []
    for i in range(f.shape[1]):
        fpfit.append([])
        for j in range(f.shape[2]):
            fpfit[-1].append(polyfit(x, f[:, i, j], plot))
    return np.array(fpfit)


def lblvals(band, T, T0=296):
    Q = pyarts.arts.lbl.Q(T, band.key.isotopologue)
    lbl_str = []
    lbl_val = []

    f0 = 0
    for line in band.data.lines:
        assert line.f0 >= f0, f"Previous f0 {f0} larger than current f0 {line.f0}"
        f0 = line.f0

        lbl_str.append(
            line.s(T, Q) * pyarts.arts.constants.c**2 / (8 * pyarts.arts.constants.pi)
        )
        lbl_val.append([])
        for i in range(len(line.ls.single_models)):
            ls = line.ls.single_models[i]
            lbl_val[-1].append(
                line.f0
                + ls.D0(T0, T, ws.atmospheric_point.pressure)
                + 1j * ls.G0(T0, T, ws.atmospheric_point.pressure)
            )

    return np.array(lbl_str), np.array(lbl_val)


def eqvvals(band, T):
    band.data.lineshape = "VP_ECS_HARTMANN"
    for line in band.data.lines:
        line.ls.one_by_one = True

    eqv_str, eqv_val = pyarts.arts.lbl.equivalent_lines(
        band, ws.ecs_data, ws.atmospheric_point, T
    )
    eqv_val = np.array(eqv_val)
    eqv_str = np.array(eqv_str)
    a = 2
    idx = np.argsort(eqv_val.real, axis=a)
    eqv_val_sorted = np.take_along_axis(eqv_val, idx, axis=a)
    eqv_str_sorted = np.take_along_axis(eqv_str, idx, axis=a)
    return eqv_val_sorted, eqv_str_sorted


def lmvalues(lbl_str, lbl_val, eqv_val_sorted, eqv_str_sorted, p):
    lm_val = eqv_val_sorted - lbl_val.T
    lm_str = 1.0 * eqv_str_sorted
    lm_str[:, 0] /= lbl_str.T
    lm_str[:, 0] -= 1.0
    lm_str[:, 1] /= lbl_str.T

    foffset = lm_val.real / p**2
    goffset = lm_val.imag / p**3
    soffset = lm_str.real / p**2
    yoffset = lm_str.imag / p
    return foffset, goffset, soffset, yoffset


def fullfit(band, T, p):
    fo, go, so, yo = lmvalues(*lblvals(band, T), *eqvvals(band, T), p)
    return fieldfit(T, fo), fieldfit(T, go), fieldfit(T, so), fieldfit(T, yo)


def adaptband(band, T, p, second_order=False):
    ff, gf, sf, yf = fullfit(band, T, p)

    for i in range(ff.shape[0]):
        for j in range(ff.shape[1]):
            band.data.lines[j].ls.single_models[i]["Y"] = pyarts.arts.TemperatureModel(
                "POLY", yf[i, j][::-1]
            )
            if second_order:
                band.data.lines[j].ls.single_models[i]["G"] = (
                    pyarts.arts.TemperatureModel("POLY", sf[i, j][::-1])
                )
                band.data.lines[j].ls.single_models[i]["DV"] = (
                    pyarts.arts.TemperatureModel("POLY", ff[i, j][::-1])
                )

    band.data.lineshape = "VP_LTE"
    return band


# WARNING: If arts-cat-data has been updated, this maybe should change
# it offsets the band so that we hit an "adaptable" band with less than N lines
di = 1
N = 60

ws = pyarts.Workspace()
ws.absorption_speciesSet(species=["CO2-626"])
ws.ReadCatalogData()

p = 1e5
ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.atmospheric_pointInit()
ws.atmospheric_point.temperature = 295  # At room temperature
ws.atmospheric_point.pressure = p
ws.atmospheric_point[pyarts.arts.SpeciesEnum("CO2")] = 400e-6
ws.atmospheric_point[pyarts.arts.SpeciesEnum("O2")] = 0.21  # At 21% Oxygen
ws.atmospheric_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 79% Nitrogen
ws.atmospheric_point.mag = [40e-6, 20e-6, 10e-6]
ws.ray_path_point

ws.jacobian_targetsInit()
ws.WignerInit()

ws.ecs_dataInit()
ws.ecs_dataAddTran2011()
ws.ecs_dataAddRodrigues1997()
ws.ecs_dataAddMeanAir(vmrs=[0.21, 0.79], species=["O2", "N2"])

f2c = pyarts.arts.convert.freq2kaycm

# This is a speed-up thing
for i in range(len(ws.absorption_bands)):
    if len(ws.absorption_bands[i].data.lines) < N:
        band = pyarts.arts.AbsorptionBand(ws.absorption_bands[i + di])
        break

ws.absorption_bands = [band]
ws.frequency_grid = np.linspace(
    ws.absorption_bands[0].data.lines[0].f0 * 0.9,
    ws.absorption_bands[0].data.lines[-1].f0 * 1.1,
    10001,
)  # around the band

# VP LTE NO LINE MIXING
ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
pm_lte = 1.0 * ws.propagation_matrix[:].T[0]

T = np.linspace(200, 320, 8)
band = adaptband(band, T, p)
ws.absorption_bands = [band]

ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
pm_adapted_lte = 1.0 * ws.propagation_matrix[:].T[0]

band.data.lineshape = "VP_ECS_HARTMANN"
ws.absorption_bands = [band]
ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
pm_full = 1.0 * ws.propagation_matrix[:].T[0]

assert np.all(pm_adapted_lte > 0), "Adaptation failed"
