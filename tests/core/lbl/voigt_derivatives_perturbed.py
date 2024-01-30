import pyarts
import numpy as np


class Setting:
    def __init__(self, pressure, zeeman, one_by_one, frange):
        self.pressure = pressure
        self.zeeman = zeeman
        self.one_by_one = one_by_one
        self.frequency_grid = frange

    def apply(self, ws, il):
        ws.atmospheric_point.pressure = self.pressure
        ws.absorption_bands[0].data.lines[il].z.on = self.zeeman
        for x in ws.absorption_bands[0].data.lines:
            x.ls.one_by_one = self.one_by_one
        ws.frequency_grid = self.frequency_grid

    def title(self, title):
        return (
            f"{title}; "
            f"pressure: {self.pressure}; "
            f"zeeman: {self.zeeman}; "
            f"one-by-one: {self.one_by_one}; "
            f"frequency_grid {self.frequency_grid[0]}-{self.frequency_grid[-1]}; "
        )


def plot_data(f, dx, dx_perturbed, title, setting):
    import matplotlib.pyplot as plt

    plt.figure(1, figsize=(8, 8))
    plt.plot(f / 1e9, dx)
    plt.plot(f / 1e9, dx_perturbed, ":")
    plt.title(setting.title(title))
    plt.show()


def assert_similarity(f, dx, dx_perturbed, title, setting):
    assert np.allclose(dx, dx_perturbed, rtol=1e-4, atol=1e-7), setting.title(
        title
    )


def plot_then_assert(f, dx, dx_perturbed, title, setting):
    plot_data(f, dx, dx_perturbed, title, setting)
    assert_similarity(f, dx, dx_perturbed, title, setting)


def pass_fn(*args):
    pass


compare_fn = assert_similarity

lc = 50474199538.7676
f1 = np.linspace(-2e6, 2e6, 10) + lc  # around the line
f2 = np.linspace(40e9, 70e9, 10)  # around the band

settings = [
    Setting(1e5, False, False, f1),
    Setting(1e0, False, False, f1),
    Setting(1e5, True, False, f1),
    Setting(1e0, True, False, f1),
    Setting(1e5, False, True, f1),
    Setting(1e0, False, True, f1),
    Setting(1e5, True, True, f1),
    Setting(1e0, True, True, f1),
    Setting(1e5, False, False, f2),
    Setting(1e0, False, False, f2),
    Setting(1e5, True, False, f2),
    Setting(1e0, True, False, f2),
    Setting(1e5, False, True, f2),
    Setting(1e0, False, True, f2),
    Setting(1e5, True, True, f2),
    Setting(1e0, True, True, f2),
]

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["O2-66"])

ws.abs_lines_per_species = pyarts.arts.ArrayOfArrayOfAbsorptionLines()
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(
    ws.abs_lines_per_species, basename="lines/"
)

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
il = 95

ws.absorption_bandsFromAbsorbtionLines(
    abs_lines_per_species=ws.abs_lines_per_species
)
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

assert np.allclose(ws.absorption_bands[0].data.lines[il].f0, lc), (
    "Line has changed, test is broken! "
    + f"lc should be set to {ws.absorption_bands[0].data.lines[il].f0}"
)

ws.WignerInit(symbol_type=3)

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.atmospheric_pointInit()
ws.atmospheric_point.temperature = 295  # At room temperature
ws.atmospheric_point[pyarts.arts.SpeciesEnum("Oxygen")] = 0.21  # At 21% atmospheric Oxygen
ws.atmospheric_point.mag = [40e-6, 20e-6, 10e-6]

for setting in settings:
    print(setting.title("Running test"))
    setting.apply(ws, il)

    ws.jacobian_targetsInit()
    ws.jacobian_targetsAddSpeciesIsotopologueRatio(species="O2-66")
    ws.jacobian_targetsAddSpeciesVMR(species="O2")
    ws.jacobian_targetsAddTemperature()
    ws.jacobian_targetsAddWindField(component="u")
    ws.jacobian_targetsAddLineParameter(
        id=bandkey, line_index=il, parameter="f0"
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey, line_index=il, parameter="e0"
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey, line_index=il, parameter="a"
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="G0",
        coefficient="X0",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="G0",
        coefficient="X1",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="G0",
        coefficient="X0",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="G0",
        coefficient="X1",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X0",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X1",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X2",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X3",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X0",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X1",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X2",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="Y",
        coefficient="X3",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="D0",
        coefficient="X0",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="D0",
        coefficient="X0",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="DV",
        coefficient="X0",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="DV",
        coefficient="X0",
        species="Bath",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="G",
        coefficient="X0",
        species="O2",
    )
    ws.jacobian_targetsAddLineParameter(
        id=bandkey,
        line_index=il,
        parameter="G",
        coefficient="X0",
        species="Bath",
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    dpm = ws.propagation_matrix_jacobian * 1.0
    pm = ws.propagation_matrix * 1.0
    ws.jacobian_targetsInit()

    # ISOTOPOLOGUE RATIO
    d = 0.0001
    key = pyarts.arts.SpeciesIsotopeRecord("O2-66")
    ws.atmospheric_point[key] += d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.atmospheric_point[key] -= d

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid,
        dpm[0][:, 0],
        dpm_dX[:, 0],
        "Isotopologue ratio",
        setting,
    )

    # VMR
    d = 0.00001
    key = pyarts.arts.SpeciesEnum("O2")
    ws.atmospheric_point[key] += d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.atmospheric_point[key] -= d

    dpm_dX = (pm_d - pm) / d
    compare_fn(ws.frequency_grid, dpm[1][:, 0], dpm_dX[:, 0], "VMR", setting)

    # Temperature
    d = 1e-6
    key = pyarts.arts.options.AtmKey.t
    ws.atmospheric_point[key] += d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.atmospheric_point[key] -= d

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[2][:, 0], dpm_dX[:, 0], "Temperature", setting
    )

    # Frequency
    d = 1e3
    orig = ws.frequency_grid * 1.0
    ws.frequency_grid = orig + d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.frequency_grid = orig

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[3][:, 0], dpm_dX[:, 0], "Frequency", setting
    )

    # line center
    d = 1e3
    orig = ws.absorption_bands[0].data.lines[il].f0 * 1.0
    ws.absorption_bands[0].data.lines[il].f0 += d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].f0 = orig

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[4][:, 0], dpm_dX[:, 0], "Line center", setting
    )

    # lower state energy
    d = 1e-26
    orig = ws.absorption_bands[0].data.lines[il].e0 * 1.0
    ws.absorption_bands[0].data.lines[il].e0 += d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].e0 = orig

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid,
        dpm[5][:, 0],
        dpm_dX[:, 0],
        "Lower state energy",
        setting,
    )

    # einstein coefficient
    d = 1e-14
    orig = ws.absorption_bands[0].data.lines[il].a * 1.0
    ws.absorption_bands[0].data.lines[il].a += d

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].a = orig

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid,
        dpm[6][:, 0],
        dpm_dX[:, 0],
        "Einstein coefficient",
        setting,
    )

    # O2 G0 X0
    d = 1e-1
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[0]
        .data.get("G0")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "G0", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "G0", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[7][:, 0], dpm_dX[:, 0], "O2 G0 X0", setting
    )

    # O2 G0 X1
    d = 1e-4
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[0]
        .data.get("G0")
    )
    data = orig.data
    data[1] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "G0", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "G0", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[8][:, 0], dpm_dX[:, 0], "O2 G0 X1", setting
    )

    # Bath G0 X0
    d = 1e-3
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[1]
        .data.get("G0")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "G0", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "G0", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[9][:, 0], dpm_dX[:, 0], "Bath G0 X0", setting
    )

    # Bath G0 X1
    d = 1e-3
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[1]
        .data.get("G0")
    )
    data = orig.data
    data[1] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "G0", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "G0", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[10][:, 0], dpm_dX[:, 0], "Bath G0 X1", setting
    )

    # O2 Y X0
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[11][:, 0], dpm_dX[:, 0], "O2 Y X0", setting
    )

    # O2 Y X1
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
    )
    data = orig.data
    data[1] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[12][:, 0], dpm_dX[:, 0], "O2 Y X1", setting
    )

    # O2 Y X2
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
    )
    data = orig.data
    data[2] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[13][:, 0], dpm_dX[:, 0], "O2 Y X2", setting
    )

    # O2 Y X3
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("Y")
    )
    data = orig.data
    data[3] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[14][:, 0], dpm_dX[:, 0], "O2 Y X3", setting
    )

    # Bath Y X0
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[15][:, 0], dpm_dX[:, 0], "Bath Y X0", setting
    )

    # Bath Y X1
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
    )
    data = orig.data
    data[1] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[16][:, 0], dpm_dX[:, 0], "Bath Y X1", setting
    )

    # O2 Y X2
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
    )
    data = orig.data
    data[2] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[17][:, 0], dpm_dX[:, 0], "Bath Y X2", setting
    )

    # O2 Y X3
    d = 1e-10
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("Y")
    )
    data = orig.data
    data[3] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "Y", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[18][:, 0], dpm_dX[:, 0], "Bath Y X3", setting
    )

    # O2 D0 X0
    d = 1e-3
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[0]
        .data.get("D0")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "D0", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "D0", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[19][:, 0], dpm_dX[:, 0], "O2 D0 X0", setting
    )

    # Bath D0 X0
    d = 1e-3
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[1]
        .data.get("D0")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "D0", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "D0", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[20][:, 0], dpm_dX[:, 0], "Bath D0 X0", setting
    )

    # O2 DV X0
    d = 1e-1 / ws.atmospheric_point.pressure
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[0]
        .data.get("DV")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "DV", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "DV", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[21][:, 0], dpm_dX[:, 0], "O2 DV X0", setting
    )

    # Bath DV X0
    d = 1e4 / ws.atmospheric_point.pressure**2
    orig = (
        ws.absorption_bands[0]
        .data.lines[il]
        .ls.single_models[1]
        .data.get("DV")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "DV", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "DV", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[22][:, 0], dpm_dX[:, 0], "Bath DV X0", setting
    )

    # O2 G X0
    d = 1e-3
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.get("G")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "G", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[0].data.set(
        "G", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[23][:, 0], dpm_dX[:, 0], "O2 G X0", setting
    )

    # Bath G X0
    d = 1
    orig = (
        ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.get("G")
    )
    data = orig.data
    data[0] += d
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "G", pyarts.arts.TemperatureModel(orig.type, data)
    )

    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines(no_negative_absorption=False)

    pm_d = ws.propagation_matrix * 1.0
    ws.absorption_bands[0].data.lines[il].ls.single_models[1].data.set(
        "G", orig
    )

    dpm_dX = (pm_d - pm) / d
    compare_fn(
        ws.frequency_grid, dpm[24][:, 0], dpm_dX[:, 0], "Bath G X0", setting
    )
