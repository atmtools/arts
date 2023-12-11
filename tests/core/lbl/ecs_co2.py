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
            line.s(T, Q)
            * pyarts.arts.constants.c**2
            / (8 * pyarts.arts.constants.pi)
        )
        lbl_val.append([])
        for i in range(len(line.ls.single_models)):
            ls = line.ls.single_models[i]
            lbl_val[-1].append(
                line.f0
                + ls.D0(T0, T, ws.atm_point.pressure)
                + 1j * ls.G0(T0, T, ws.atm_point.pressure)
            )

    return np.array(lbl_str), np.array(lbl_val)


def eqvvals(band, T):
    band.data.lineshape = "VP_ECS_HARTMANN"
    for line in band.data.lines:
        line.ls.one_by_one = True

    eqv_str, eqv_val = pyarts.arts.lbl.equivalent_lines(
        band, ws.ecs_data2, ws.atm_point, T
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
    lm_str[:, 0] -= lbl_str.T
    lm_str[:, 1] -= lbl_str.T

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
            band.data.lines[j].ls.single_models[i].data.set(
                "Y", pyarts.arts.TemperatureModel("POLY", yf[i, j][::-1])
            )
            if second_order:
                band.data.lines[j].ls.single_models[i].data.set(
                    "G", pyarts.arts.TemperatureModel("POLY", sf[i, j][::-1])
                )
                band.data.lines[j].ls.single_models[i].data.set(
                    "DV", pyarts.arts.TemperatureModel("POLY", ff[i, j][::-1])
                )
    return band


i = 4

ws = pyarts.Workspace()
ws.abs_speciesSet(species=["CO2-626"])
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")
ws.absorption_bandsFromAbsorbtionLines()

p = 1e5
ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.select_abs_species = []  # All species
ws.atm_pointInit()
ws.atm_point.temperature = 295  # At room temperature
ws.atm_point.pressure = p
ws.atm_point[ws.abs_species[0]] = 400e-6
ws.atm_point[pyarts.arts.SpeciesEnum("O2")] = 0.21  # At 21% Oxygen
ws.atm_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 79% Nitrogen
ws.atm_point.mag = [40e-6, 20e-6, 10e-6]

ws.jacobian_targetsInit()
ws.Wigner6Init()

ws.ecs_dataInitNEWNEW()
ws.ecs_dataAddTran2011NEWNEW()
ws.ecs_dataAddRodrigues1997NEWNEW()
ws.ecs_dataAddMeanAirNEWNEW(vmrs=[0.21, 0.79], species=["O2", "N2"])

f2c = pyarts.arts.convert.freq2kaycm

y = pyarts.arts.AbsorptionBands(ws.absorption_bands)

band = y[i]
ws.absorption_bands = [band]
ws.f_grid = np.linspace(
    ws.absorption_bands[0].data.lines[0].f0 * 0.9,
    ws.absorption_bands[0].data.lines[-1].f0 * 1.1,
    10001,
)  # around the band

# VP LTE NO LINE MIXING
ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()
pm_lte = 1.0 * ws.propmat_clearsky[:].T[0]

T = np.linspace(200, 320, 8)
band = adaptband(band, T, p)
ws.absorption_bands = [band]

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()
pm_full = 1.0 * ws.propmat_clearsky[:].T[0]

band.data.lineshape = "VP_LTE"
ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()
pm_adapted_lte = 1.0 * ws.propmat_clearsky[:].T[0]

assert np.all(pm_adapted_lte > 0), "Adaptation failed"
