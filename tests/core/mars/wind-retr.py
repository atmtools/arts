import pyarts

import time
import copy
import tqdm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def covmat1d_markov(y, sigma, scale=1.0):
    y = np.atleast_1d(y)
    sigma = np.atleast_1d(sigma)
    scale = np.atleast_1d(scale)
    n = len(y)

    if len(sigma) == 1:
        sigma = np.ones(n) * sigma[0]
    elif len(sigma) != n:
        raise ValueError("sigma must be a scalar or an array of the same length as y")
    
    if len(scale) == 1:
        scale = np.ones(n) * scale[0]
    elif len(scale) != n:
        raise ValueError("scale must be a scalar or an array of the same length as y")

    out = np.stack([y] * n)
    out = np.abs(np.subtract(out.T, y)).T
    out = np.einsum('ij,i->ij', out, 1.0 / sigma)
    return np.einsum('ij,i->ij', np.exp(-out), scale)

def gaussgrid(sigma, N=11, lim=0.01):
    assert lim <= 0.5 and lim > 0.0, "lim must be (0, 0.5]"
    assert N > 1, "N must be greater than 1 to have a central point"
    assert N % 2, "N must be odd to keep central point"

    X = np.sqrt(-2 * np.log(lim))

    return np.linspace(-X, X, N) * sigma


ws = pyarts.Workspace()

ws.absorption_speciesSet(
    species=[
        "H2O-161",
        "CO-26",
        "O2-66",
        "PH3-1111",
        "H2O2-1661",
        "HCN-124",
        "NH3-4111",
        "CO2-626",
        "CO2-CIA-CO2",
    ]
)

ws.surface_fieldMars()
ws.surface_field["t"] = pyarts.arts.GriddedField2.fromxml(
    "planets/Mars/Ls0.day.dust-medium/surface_temperature.xml"
)
ws.surface_field["h"] = pyarts.arts.GriddedField2.fromxml(
    "planets/Mars//surface_elevation.xml"
)
ws.surface_field["t"].set_extrapolation("Nearest")
ws.surface_field["h"].set_extrapolation("Nearest")

species = [
    "H2O-161",
    "CO-26",
    "O2-66",
    "PH3-1111",
    "H2O2-1661",
    "HCN-124",
    "NH3-4111",
    "CO2-626",
]

# x = {}
ws.absorption_bandsReadSpeciesSplitARTSCAT(basename="spectroscopy/Perrin/")
for band in ws.absorption_bands:
    ws.absorption_bands[band].cutoff = "ByLine"
    ws.absorption_bands[band].cutoff_value = 50.0e9
#     if str(f"{band.isotopologue}") in species:
#         x[band] = ws.absorption_bands[band]
# ws.absorption_bands = x

ws.absorption_cia_dataReadSpeciesSplitCatalog(
    basename="cia/",
)

# Remove D0 lines from all absorption because it looks weird.  Remember question about it!
for x in ws.absorption_bands:
    for line in ws.absorption_bands[x].lines:
        line.ls.remove("D0")

print("LBL data:", [x.isotopologue for x in ws.absorption_bands], sep="\n")

print("CIA data:", [x.specs for x in ws.absorption_cia_data], sep="\n")

print()

ws.propagation_matrix_agendaAuto(T_extrapolfac=10)
ws.atmospheric_fieldRead(
    toa=150e3,
    basename="planets/Mars/Ls0.day.dust-medium/Ls0.day.dust-medium.sol-avg/",
    missing_is_zero=True,
    extrapolation="Nearest",
)
ws.atmospheric_fieldAppendLineIsotopologueData(
    basename="planets/Mars/isotopologue_ratios/", replace_existing=True
)

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.spectral_radiance_observer_agendaSet(option="Emission")
ws.ray_path_observer_agendaSetGeometric(max_step=10e3)


@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def inversion_iterate_agenda(ws):
    ws.UpdateModelStates()
    ws.measurement_vectorFromSensor()
    ws.measurement_vector_fittedFromMeasurement()


NP = 31
NF = 1001
NL = 4
pos = [1000e3, 0, 0]
limb_alt = np.linspace(10, 40, NL) * 1e3
antenna_size = 30.0e-2  # 0.30 m
delta_t = 120.0
TB_600 = 1500.0  # K
TB_1200 = 1.5 * TB_600  # K  ???
limb_zas = [
    pyarts.arts.path.geometric_tangent_zenith(pos, ws.surface_field.ellipsoid, alt, 180)
    for alt in limb_alt
]
limb_los = [pyarts.arts.path.mirror([za, 180.0]) for za in limb_zas]

df = 0.025
wdf = 2 * df
fp = (
    [
        547,
        547.66 - 2 * wdf,
        547.693 + 2 * wdf,
    ],
    [
        548,
        548.555 - df,
        548.555 + df,
        548.831 - df,
        548.831 + df,
    ],
    [
        550,
        550.152 - df,
        550.152 + df,
        550.923 - wdf,
        550.929 + wdf,
    ],
    [
        551.5,
        551.672 - df,
        551.672 + df,
        552.021 - wdf,
        552.021 + wdf,
    ],
    [
        561,
        561.115 - wdf,
        561.115 + wdf,
        561.713 - df,
        561.713 + df,
    ],
    [
        599.5,
        599.727 - df,
        599.727 + df,
        599.922 - 2 * wdf,
        599.932 + 2 * wdf,
    ],
    [
        620,
        620.682 - 2 * wdf,
        620.72 + 2 * wdf,
    ],
    [
        1101,
        1101.342 - wdf,
        1101.358 + wdf,
        1101.68 - 2 * wdf,
        1101.714 + 2 * wdf,
        1101.899 - 2 * df,
        1101.899 + 2 * df,
    ],
    [
        1107,
        1107.167 - wdf,
        1107.167 + wdf,
        1107.43 - df,
        1107.43 + df,
    ],
    [
        1120,
        1120.715 - df,
        1120.715 + df,
    ],
    [
        1122.5,
        1122.902 - df,
        1122.902 + df,
    ],
    [
        1146,
        1146.621 - wdf,
        1146.621 + wdf,
    ],
    [
        1148.5,
        1148.969 - wdf,
        1148.981 + wdf,
    ],
    [
        1167.7,
        1167.897 - df,
        1167.897 + df,
        1168.136 - df,
        1168.136 + df,
        1168.359 - df,
        1168.359 + df,
        1168.607 - df,
        1168.607 + df,
    ],
    [
        1172,
        1172.524 - wdf,
        1172.528 + wdf,
    ],
    [
        1179.5,
        1179.524 - df,
        1179.524 + df,
        1179.651 - df,
        1179.651 + df,
        1179.879 - df,
        1179.879 + df,
        1180.324 - df,
        1180.324 + df,
    ],
    [
        1181,
        1181.39 - 2 * wdf,
        1181.4 + 2 * wdf,
    ],
    [
        1188.5,
        1188.856 - 2 * wdf,
        1188.873 + 2 * wdf,
        1189.421 - df,
        1189.421 + df,
    ],
    [
        1198.5,
        1198.991 - 2 * wdf,
        1199.017 + 2 * wdf,
    ],
    [
        1211,
        1211.322 - wdf,
        1211.337 + wdf,
        1211.799 - df,
        1211.799 + df,
    ],
    [
        1212.5,
        1212.98 - 2 * wdf,
        1212.98 + 2 * wdf,
    ],
    [
        1217,
        1217.242 - 2 * wdf,
        1217.274 + 2 * wdf,
    ],
)

wind_v = copy.copy(ws.atmospheric_field["wind_v"])
wind_v.data.data = np.zeros((NP, 1, 1)) + 10
wind_v.data.grids = (np.linspace(0, 60e3, NP), [0], [0])  # Altitude from 0 to 50 km
wind_v.set_extrapolation("Nearest")

z = wind_v.data.grids[0] / 1e3  # Convert to km
scale = 100.0
makarov = covmat1d_markov(z, sigma=5.0, scale=scale)
makarov[makarov < (scale * 0.01)] = 0.0
makarov = sp.sparse.csc_matrix(makarov)
print("Elements:", makarov.nnz)

for f0r in fp:  # tqdm.tqdm(fp):
    f0 = f0r[0]
    t = time.time()
    print(f"Begin with {f0}")

    ws.frequency_grid = (f0 + np.linspace(0.0, 1.0, NF)) * 1e9
    if len(f0r) > 1:
        franges = []
        for i in range(1, len(f0r), 2):
            franges.append(
                1.0
                * ws.frequency_grid[
                    np.logical_and(
                        ws.frequency_grid >= f0r[i] * 1e9,
                        ws.frequency_grid <= f0r[i + 1] * 1e9,
                    )
                ]
            )
        ws.frequency_grid = np.concatenate(franges)
    else:
        franges = [1.0 * ws.frequency_grid]

    # ws.measurement_sensorSimple(pos=pos, los=limb_los[0])

    hphw = np.rad2deg(
        pyarts.arts.convert.freq2wavelen(np.mean(ws.frequency_grid))
        * sp.special.jn_zeros(1, 1)
        / antenna_size
        / np.pi
    )
    hphw = hphw[0]
    sigma = hphw / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    zas = gaussgrid(sigma=sigma, lim=0.01, N=5)

    res = 1 / (sigma * np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * zas**2 / sigma**2)

    sensor_raw = pyarts.arts.SortedGriddedField1(
        name="Sensor Resolution", data=res, grids=[zas], grid_names=["dza"]
    )

    ws.measurement_sensorInit()
    for i in range(NL):
        ws.measurement_sensorAddRawSensor(
            pos=pos,
            los=limb_los[i],
            raw_sensor_perturbation=sensor_raw,
            normalize=True,
        )

    ws.measurement_sensorMakeExhaustive()
    print(f"measurement_sensorMakeExhaustive done after {time.time() - t:.2f} seconds")

    ws.jacobian_targetsOff()
    ws.measurement_vectorFromSensor()
    print(f"measurement_vectorFromSensor done after {time.time() - t:.2f} seconds")

    ws.atmospheric_field["wind_v"] = wind_v

    ws.RetrievalInit()
    ws.RetrievalAddWindField(component="v", matrix=makarov)
    ws.RetrievalFinalizeDiagonal()

    TB = TB_600 if f0 < 1000 else TB_1200
    noise = TB / np.sqrt(1e9* delta_t / NF)

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.atmospheric_field["wind_v"].data += 10.0
    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
    NM = len(ws.measurement_vector)
    ws.measurement_vector += np.random.normal(0, noise, NM)

    ws.OEM(method="gn")
    print(f"OEM done after {time.time() - t:.2f} seconds")
    ws.measurement_averaging_kernelCalc()
    print(f"measurement_averaging_kernelCalc done after {time.time() - t:.2f} seconds")
    ws.measurement_vector_error_covariance_matrix_observation_systemCalc()
    print(
        f"measurement_vector_error_covariance_matrix_observation_systemCalc done after {time.time() - t:.2f} seconds"
    )

    plt.figure(1, figsize=(10, 15))
    plt.clf()

    plt.subplot(7, 1, 1)
    for i in range(len(franges)):
        fs = franges[i] / 1e9
        i0 = sum([len(x) for x in franges[:i]])
        for j in range(NL):
            plt.plot(
                fs,
                ws.measurement_vector.reshape(NL, NM // NL)[j, i0 : i0 + len(fs)],
                "b",
                label="Measurement Vector",
            )
            plt.plot(
                fs,
                ws.measurement_vector_fitted.reshape(NL, NM // NL)[
                    j, i0 : i0 + len(fs)
                ],
                "r",
                label="Fitted Measurement Vector",
            )
            if i == 0 and j == 0:
                plt.legend()

    plt.subplot(7, 1, 2)
    for i in range(len(franges)):
        fs = franges[i] / 1e9
        i0 = sum([len(x) for x in franges[:i]])
        for j in range(NL):
            plt.plot(
                fs,
                ws.measurement_vector.reshape(NL, NM // NL)[j, i0 : i0 + len(fs)]
                - ws.measurement_vector_fitted.reshape(NL, NM // NL)[
                    j, i0 : i0 + len(fs)
                ],
            )
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Residual [K]")

    J = 1.0 * ws.measurement_jacobian
    J[np.abs(J) < np.max(np.abs(J)) * 0.001] = 0.0  # Set small values to zero
    plt.subplot(7, 1, 3)
    plt.contourf(J.T, 50)
    plt.xlabel("Measurement [#]")
    plt.ylabel("Altitude Index [#]")
    cbar = plt.colorbar()
    cbar.set_label("Jacobian [K / m/s]")

    plt.subplot(7, 1, 4)
    plt.plot(wind_v.data.flatten(), z, label="True Wind V")
    plt.plot(ws.model_state_vector, z, label="Retrieved Wind V")
    plt.plot(ws.model_state_vector_apriori, z, label="Apriori Wind V")
    plt.xlabel("Wind V [m/s]")
    plt.ylabel("Altitude [km]")
    plt.legend()

    plt.subplot(7, 2, 9)
    plt.plot(ws.measurement_averaging_kernel, z)
    plt.plot(ws.measurement_averaging_kernel @ np.ones(len(z)), z, "k")
    plt.ylabel("Altitude [km]")
    plt.xlabel("Averaging Kernel [-]")

    plt.subplot(7, 2, 10)
    plt.plot(ws.measurement_averaging_kernel.T, z)
    plt.plot(ws.measurement_averaging_kernel.T @ np.ones(len(z)), z, "k")
    plt.ylabel("Altitude [km]")
    plt.xlabel("Averaging Kernel Transpose [-]")

    cv = 1.0 * ws.measurement_vector_error_covariance_matrix_observation_system
    cv[np.abs(cv) < np.max(np.abs(cv)) * 0.001] = 0.0  # Set small values to zero
    plt.subplot(7, 1, 6)
    plt.contourf(z, z, cv, 50)
    plt.xlabel("Altitude [km]")
    plt.ylabel("Altitude [km]")
    cbar = plt.colorbar()
    cbar.set_label("Covariance [K^2]")

    plt.subplot(7, 2, 13)
    plt.plot(
        ws.measurement_averaging_kernel
        @ (ws.model_state_vector - ws.model_state_vector_apriori),
        z,
        "k",
    )
    plt.ylabel("Altitude [km]")
    plt.xlabel("Smoothing Error [m/s]")

    plt.subplot(7, 2, 14)
    plt.plot(
        ws.measurement_gain_matrix
        @ (ws.measurement_vector - ws.measurement_vector_fitted),
        z,
        "k",
    )
    plt.ylabel("Altitude [km]")
    plt.xlabel("Measurement Error [m/s]")

    plt.tight_layout()
    plt.savefig(f"/home/richard/Work/test_figs/{f0}.png")

    print(f"Done with {f0} - Took {time.time() - t:.2f} seconds")
