"""LM through workspace covariance methods: diagonal, mixed and correlated.

The forward model is deliberately affine, so the converged LM state is the
exact linear OEM solution, tabulated below for each pattern. Dense/Sparse
storage and supplied inverses must not change that solution. Changing
correlations can, and the plots show that difference.
Run normally for plots; ARTS_HEADLESS suppresses them in the test suite.
"""

import os

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
K = np.array([[1, .3, .2, 0], [0, 1, .4, .2], [.2, 0, 1, .3], [.4, .2, 0, 1]])
SCALE = np.array([.002, .002, .03, .03])
OBSERVED = np.array([1., -.5, .3, .8])
START = np.array([1., -1., 1., -1.])
PATTERNS = {
    "Diagonal": (0., 0.),
    "Mixed": (0., .6),
    "Correlated": (.4, .6),
}
# The exact linear MAP and gain in standardized coordinates (Sa = I), to 12
# decimals, of each pattern: x = (I + K^T Se^-1 K)^-1 K^T Se^-1 y.
EXPECTED = {
    "Diagonal": (
        [0.941171708459, -0.388454324853, -0.0289045677, 0.379294550562],
        [
            [0.759143428093, -0.217575791669, -0.03099934523, 0.103175235126],
            [0.097772334874, 0.759143428093, -0.272851564989, -0.03099934523],
            [-0.03099934523, 0.103175235126, 0.759143428093, -0.217575791669],
            [-0.272851564989, -0.03099934523, 0.097772334874, 0.759143428093],
        ],
    ),
    "Mixed": (
        [0.942772991361, -0.368139607258, -0.083412933314, 0.367060360229],
        [
            [0.760337375657, -0.215163295463, -0.046349704238, 0.110948599055],
            [0.089439267864, 0.753600759851, -0.278447664431, 0.003444755166],
            [-0.013910289711, 0.109933558152, 0.816613493659, -0.324399890781],
            [-0.254857539289, -0.007684785304, -0.003383709074, 0.773863274485],
        ],
    ),
    "Correlated": (
        [1.018384900177, -0.422650131251, -0.09342461928, 0.348212311658],
        [
            [0.795279164573, -0.283781207597, -0.035028431878, 0.114654576711],
            [0.024066364712, 0.761661883471, -0.267802609613, 0.01806903582],
            [-0.024005080308, 0.11338431816, 0.817673765035, -0.322536886752],
            [-0.259450597742, 0.013662336406, -0.008132708624, 0.771167362738],
        ],
    ),
}


def workspace():
    ws = pyarts.Workspace()
    for name, group in (
        ("atm_field", arts.AtmField),
        ("surf_field", arts.SurfaceField),
        ("subsurf_field", arts.SubsurfaceField),
        ("abs_bands", arts.AbsorptionBands),
        ("measurement_sensors", arts.ArrayOfSensor),
    ):
        setattr(ws, name, group())
    ws.surf_field.ellipsoid = [1, 1]
    for species, values in (("H2O", [.01, .005]), ("O2", [.21, .20])):
        ws.atm_field[species] = arts.GriddedField3(
            data=np.array(values).reshape(2, 1, 1),
            grid_names=["Altitude", "Latitude", "Longitude"],
            grids=[[0, 1000], [0], [0]],
        )
    ws.oemInit()
    ws.oemAddSpeciesVMR(species="H2O", matrix=np.eye(2) * SCALE[0] ** 2)
    ws.oemAddSpeciesVMR(species="O2", matrix=np.eye(2) * SCALE[2] ** 2)
    ws.jac_targetsFinalize()
    ws.measurement_vec = OBSERVED
    ws.oem.uncheck()
    ws.oem.model_state_vec = []

    ws.measurement_jac = arts.Matrix()
    ws.model_state_vecFromData()
    ws.measurement_sensorsSimple(freq_grid=np.arange(1.0, len(OBSERVED) + 1), pos=[0, 0, 0], los=[0, 0])
    ws.oemSetApriori()
    ws.oemSetMeasurement()
    ws.oemMeasurementCovmatConstant(value=1.)
    ws.oemFinalizeDiagonal()
    prior = np.array(ws.oem.model_state_vec_apriori, copy=True)

    def forward(atm, bands, sensors, surf, subsurf, state_targets, targets, state):
        fit = K @ ((np.asarray(state) - prior) / SCALE)
        jac = K / SCALE if targets.x_size() else np.empty((0, 0))
        return atm, bands, sensors, surf, subsurf, arts.Vector(fit), arts.Matrix(jac)

    ws.inversion_iterate_agendaSetOperator(f=forward)
    return ws, prior


def measurement_covariance(ws, correlations, sparse, supplied_inverse):
    # Build two measurement covariance blocks through the owning workspace
    # helpers, varying storage and supplied inverses without editing caches.
    ws.oemMeasurementCovmatInit()
    for rho in correlations:
        block = .25 * np.array([[1., rho], [rho, 1.]])
        storage = arts.Sparse if sparse else arts.Matrix
        options = {"inverse": storage(np.linalg.solve(
            block, np.eye(2)))} if supplied_inverse else {}
        ws.oemMeasurementCovmatAdd(
            matrix=storage(block), **options,
        )


def run():
    ws, prior = workspace()
    results = []
    for name, correlations in PATTERNS.items():
        reference = None
        for sparse in (False, True):
            for supplied_inverse in (False, True):
                measurement_covariance(ws, correlations, sparse, supplied_inverse)
                expected, expected_gain = map(np.array, EXPECTED[name])
                # First call prepares Se; the second reuses it. Both start from
                # the same deliberately perturbed state, not the previous fit.
                for repeat in range(2):
                    ws.oem.uncheck()
                    ws.oem.model_state_vec = prior + SCALE * START
                    ws.oem.measurement_vec_fit = []
                    ws.oem.measurement_jac = arts.Matrix()
                    ws.oemCheck()
                    ws.oemCalc(settings=arts.OptimalEstimationSettings(method="lm", lm=arts.LevenbergMarquardtSettings(),
                           max_iter=50, stop_dx=1e-12))
                    assert not len(ws.oem.diagnostics.errors), str(
                        ws.oem.diagnostics.errors)
                    assert ws.oem.diagnostics.status == pyarts.arts.OptimalEstimationStatus.Converged, ws.oem.diagnostics
                    fitted = (np.array(ws.oem.model_state_vec) - prior) / SCALE
                    fit_y = np.array(ws.oem.measurement_vec_fit, copy=True)
                    np.testing.assert_allclose(fitted, expected, rtol=0, atol=2e-7)
                    np.testing.assert_allclose(fit_y, K @ expected, rtol=0, atol=2e-7)
                    np.testing.assert_allclose(
                        np.array(ws.oem.measurement_gain_mat) / SCALE[:, None],
                        expected_gain, rtol=0, atol=1e-10,
                    )
                    if reference is None:
                        reference = fitted.copy()
                    np.testing.assert_allclose(fitted, reference, rtol=0, atol=2e-7)
                label = f"{name}: {'Sparse' if sparse else 'Matrix'}, {'inverse' if supplied_inverse else 'automatic'}"
                results.append((name, label, fitted.copy(), fit_y))
                print(label, "->", fitted)
    print(pyarts.retrieval.information_from_workspace(ws))

    return results


def plot(results):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    coordinate = np.arange(4)
    axes[0].plot(coordinate, np.zeros(4), "k--", label="A priori")
    axes[0].plot(coordinate, START, "k:", marker="x", label="Manipulated start")
    axes[1].plot(coordinate, OBSERVED, "kx", ms=9, label="Measurements")
    axes[1].plot(coordinate, K @ START, "k:", label="Manipulated start")
    colors = dict(zip(PATTERNS, ("tab:blue", "tab:orange", "tab:green")))
    for i, (name, label, state, fit_y) in enumerate(results):
        # Four equivalent representations overlap for each covariance pattern.
        for ax, values in zip(axes, (state, fit_y)):
            ax.plot(coordinate, values, color=colors[name], marker=("o", "s", "^", "D")[i % 4],
                    fillstyle="none", linestyle="-", label=name if i % 4 == 0 else None)
    axes[0].set_xticks(coordinate, ["H2O[0]", "H2O[1]", "O2[0]", "O2[1]"])
    axes[0].set_ylabel("State departure from prior / prior standard deviation")
    axes[1].set_xticks(coordinate, ["Channel 0", "Channel 1", "Channel 2", "Channel 3"])
    axes[1].set_ylabel("Synthetic measurement")
    for ax in axes:
        ax.grid(True)
        ax.legend(fontsize="small")
    fig.suptitle(
        "Workspace LM: covariance representations agree (four overlapping markers per pattern)")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    results = run()
    if "ARTS_HEADLESS" not in os.environ:
        plot(results)
        import matplotlib.pyplot as plt
        plt.show()
