import pyarts
import numpy as np


def plot_data(f, dx, dx_perturbed, title):
    import matplotlib.pyplot as plt

    plt.plot(f, dx)
    plt.plot(f, dx_perturbed, ":")
    plt.title(title)
    plt.show()


def assert_similarity(f, dx, dx_perturbed, title):
    assert np.allclose(dx, dx_perturbed, rtol=1e-5), title


compare_fn = assert_similarity

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["O2-66"])

ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
il = 95

ws.absorption_bandsFromAbsorbtionLines()
ws.absorption_bandsSelectFrequency(fmax=120e9)
ws.absorption_bandsKeepID(id=bandkey)

# Write some data to fields that does not exist
x = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("D0")
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("D0", x)
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("DV", x)
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("G", x)
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("D0", x)
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("DV", x)
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("G", x)

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.select_abs_species = []  # All species
ws.f_grid = np.linspace(-2e9, 2e9, 5) + 50474199538.7676
ws.atm_pointInit()
ws.atm_point.temperature = 295  # At room temperature
ws.atm_point.pressure = 1e5  # At 1 bar
ws.atm_point[ws.abs_species[0]] = 0.21  # At 21% atmospheric Oxygen


ws.jacobian_targetsInit()
ws.jacobian_targetsAddSpeciesIsotopologueRatio(species="O2-66")
ws.jacobian_targetsAddSpeciesVMR(species="O2")
ws.jacobian_targetsAddTemperature()
ws.jacobian_targetsAddWindField(component="u")
ws.jacobian_targetsAddLineParameter(id=bandkey, line_index=il, parameter="f0")
ws.jacobian_targetsAddLineParameter(id=bandkey, line_index=il, parameter="e0")
ws.jacobian_targetsAddLineParameter(id=bandkey, line_index=il, parameter="a")
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="G0", coefficient="X0", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="G0", coefficient="X1", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="G0", coefficient="X0", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="G0", coefficient="X1", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X0", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X1", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X2", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X3", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X0", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X1", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X2", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="Y", coefficient="X3", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="D0", coefficient="X0", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="D0", coefficient="X0", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="DV", coefficient="X0", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="DV", coefficient="X0", species="Bath"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="G", coefficient="X0", species="O2"
)
ws.jacobian_targetsAddLineParameter(
    id=bandkey, line_index=il, parameter="G", coefficient="X0", species="Bath"
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

dpm = ws.dpropmat_clearsky_dx * 1.0
pm = ws.propmat_clearsky * 1.0
ws.jacobian_targetsInit()

# ISOTOPOLOGUE RATIO
d = 0.0001
key = pyarts.arts.SpeciesIsotopeRecord("O2-66")
ws.atm_point[key] += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.atm_point[key] -= d

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[0][:, 0], dpm_dX[:, 0], "Isotopologue ratio")

# VMR
d = 0.0001
key = pyarts.arts.Species("O2")
ws.atm_point[key] += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.atm_point[key] -= d

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[1][:, 0], dpm_dX[:, 0], "VMR")

# Temperature
d = 1e-6
key = pyarts.arts.options.AtmKey.t
ws.atm_point[key] += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.atm_point[key] -= d

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[2][:, 0], dpm_dX[:, 0], "Temperature")

# Frequency
d = 1e5
orig = ws.f_grid * 1.0
ws.f_grid += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.f_grid = orig

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[3][:, 0], dpm_dX[:, 0], "Frequency")

# line center
d = 1e3
orig = ws.absorption_bands[0].data.lines[il].f0 * 1.0
ws.absorption_bands[0].data.lines[il].f0 += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].f0 = orig

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[4][:, 0], dpm_dX[:, 0], "Line center")

# lower state energy
d = 1e-26
orig = ws.absorption_bands[0].data.lines[il].e0 * 1.0
ws.absorption_bands[0].data.lines[il].e0 += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].e0 = orig

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[5][:, 0], dpm_dX[:, 0], "Lower state energy")

# einstein coefficient
d = 1e-14
orig = ws.absorption_bands[0].data.lines[il].a * 1.0
ws.absorption_bands[0].data.lines[il].a += d

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].a = orig

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[6][:, 0], dpm_dX[:, 0], "Einstein coefficient")

# O2 G0 X0
d = 1
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("G0")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "G0", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("G0", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[7][:, 0], dpm_dX[:, 0], "O2 G0 X0")

# O2 G0 X1
d = 1e-4
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("G0")
data = orig.data
data[1] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "G0", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("G0", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[8][:, 0], dpm_dX[:, 0], "O2 G0 X1")

# Bath G0 X0
d = 1
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("G0")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "G0", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("G0", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[9][:, 0], dpm_dX[:, 0], "Bath G0 X0")

# Bath G0 X1
d = 1e-4
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("G0")
data = orig.data
data[1] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "G0", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("G0", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[10][:, 0], dpm_dX[:, 0], "Bath G0 X1")

# O2 Y X0
d = 1e-10
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[11][:, 0], dpm_dX[:, 0], "O2 Y X0")

# O2 Y X1
d = 1e-12
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
data = orig.data
data[1] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[12][:, 0], dpm_dX[:, 0], "O2 Y X1")


# O2 Y X2
d = 1e-14
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
data = orig.data
data[2] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[13][:, 0], dpm_dX[:, 0], "O2 Y X2")

# O2 Y X3
d = 1e-16
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
data = orig.data
data[3] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[14][:, 0], dpm_dX[:, 0], "O2 Y X3")

# Bath Y X0
d = 1e-10
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[15][:, 0], dpm_dX[:, 0], "Bath Y X0")

# Bath Y X1
d = 1e-12
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
data = orig.data
data[1] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[16][:, 0], dpm_dX[:, 0], "Bath Y X1")


# O2 Y X2
d = 1e-14
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
data = orig.data
data[2] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[17][:, 0], dpm_dX[:, 0], "Bath Y X2")

# O2 Y X3
d = 1e-16
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
data = orig.data
data[3] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "Y", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("Y", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[18][:, 0], dpm_dX[:, 0], "Bath Y X3")

# O2 D0 X0
d = 1
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("D0")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "D0", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("D0", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[19][:, 0], dpm_dX[:, 0], "O2 D0 X0")

# Bath D0 X0
d = 1
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("D0")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "D0", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("D0", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[20][:, 0], dpm_dX[:, 0], "Bath D0 X0")

# O2 DV X0
d = 1e-4
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("DV")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "DV", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("DV", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[21][:, 0], dpm_dX[:, 0], "O2 DV X0")

# Bath DV X0
d = 1e-4
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("DV")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "DV", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("DV", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[22][:, 0], dpm_dX[:, 0], "Bath DV X0")

# O2 G X0
d = 1
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("G")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
    "G", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set("G", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[23][:, 0], dpm_dX[:, 0], "O2 G X0")

# Bath G X0
d = 1
orig = ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("G")
data = orig.data
data[0] += d
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
    "G", pyarts.arts.TemperatureModel(orig.type, data)
)

ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
ws.propmat_clearskyAddLines2()

pm_d = ws.propmat_clearsky * 1.0
ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set("G", orig)

dpm_dX = (pm_d - pm) / d
compare_fn(ws.f_grid, dpm[24][:, 0], dpm_dX[:, 0], "Bath G X0")
