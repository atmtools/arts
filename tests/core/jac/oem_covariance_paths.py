"""LM through workspace covariance methods: diagonal, mixed and correlated.

The forward model is deliberately affine, so NumPy supplies an independent
exact OEM solution. Dense/Sparse storage and supplied inverses must not change
that solution. Changing correlations can, and the plots show that difference.
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


def workspace():
    ws = pyarts.Workspace()
    for name, group in (
        ("atm_field", arts.AtmField),
        ("surf_field", arts.SurfaceField),
        ("subsurf_field", arts.SubsurfaceField),
        ("abs_bands", arts.AbsorptionBands),
        ("measurement_sensor", arts.ArrayOfSensorObsel),
    ):
        setattr(ws, name, group())
    ws.surf_field.ellipsoid = [1, 1]
    for species, values in (("H2O", [.01, .005]), ("O2", [.21, .20])):
        ws.atm_field[species] = arts.GriddedField3(
            data=np.array(values).reshape(2, 1, 1),
            grid_names=["Altitude", "Latitude", "Longitude"],
            grids=[[0, 1000], [0], [0]],
        )
    ws.RetrievalInit()
    ws.RetrievalAddSpeciesVMR(species="H2O", matrix=np.eye(2) * SCALE[0] ** 2)
    ws.RetrievalAddSpeciesVMR(species="O2", matrix=np.eye(2) * SCALE[2] ** 2)
    ws.RetrievalFinalizeDiagonal()
    ws.model_state_vec_aprioriFromData()
    prior = np.array(ws.model_state_vec_apriori, copy=True)

    def forward(local):
        x = np.asarray(local.get("model_state_vec"))
        fit = K @ ((x - prior) / SCALE)
        jac = K / SCALE if local.get("jac_targets").x_size() else np.empty((0, 0))
        local.get("measurement_vec_fit").value = arts.Vector(fit)
        local.get("measurement_jac").value = arts.Matrix(jac)

    agenda = arts.Agenda("inversion_iterate_agenda")
    agenda.add(arts.Method("affine_measurements", arts.CallbackOperator(
        forward, ["model_state_vec", "jac_targets"],
        ["measurement_vec_fit", "measurement_jac"],
    )))
    agenda.finalize(True)
    ws.inversion_iterate_agenda = agenda
    ws.measurement_vec = OBSERVED
    return ws, prior


def measurement_covariance(ws, correlations, sparse, supplied_inverse):
    # This toy problem has two measurement groups matching the two retrieval
    # blocks in size and offset. Redirect the existing workspace block builder
    # to Se. No private cache manipulation or direct .blocks edits are needed.
    ws.model_state_covmatInit(model_state_covmat=ws.measurement_vec_error_covmat)
    se = np.zeros((4, 4))
    for i, (species, rho) in enumerate(zip(("H2O", "O2"), correlations)):
        block = .25 * np.array([[1., rho], [rho, 1.]])
        se[2*i:2*i+2, 2*i:2*i+2] = block
        storage = arts.Sparse if sparse else arts.Matrix
        options = {"inverse": storage(np.linalg.solve(block, np.eye(2)))} if supplied_inverse else {}
        ws.model_state_covmatAddSpeciesVMR(
            model_state_covmat=ws.measurement_vec_error_covmat,
            species=species, matrix=storage(block), **options,
        )
    return se


def run():
    ws, prior = workspace()
    results = []
    for name, correlations in PATTERNS.items():
        reference = None
        for sparse in (False, True):
            for supplied_inverse in (False, True):
                se = measurement_covariance(ws, correlations, sparse, supplied_inverse)
                # In standardized coordinates Sa=I. Solve the independent
                # analytic normal equations; do not use an ARTS inverse/cache.
                weighted_k = np.linalg.solve(se, K)
                expected = np.linalg.solve(np.eye(4) + K.T @ weighted_k,
                                           K.T @ np.linalg.solve(se, OBSERVED))
                expected_gain = np.linalg.solve(np.eye(4) + K.T @ weighted_k,
                                                weighted_k.T)
                # First call prepares Se; the second reuses it. Both start from
                # the same deliberately perturbed state, not the previous fit.
                for repeat in range(2):
                    ws.model_state_vec = prior + SCALE * START
                    ws.measurement_vec_fit = []
                    ws.measurement_jac = arts.Matrix()
                    ws.OEM(method="lm", lm_ga_settings=arts.LevenbergMarquardtSettings(),
                           max_iter=50, stop_dx=1e-12)
                    assert not len(ws.errors), str(ws.errors)
                    assert ws.oem_diagnostics[0] == 0, ws.oem_diagnostics
                    fitted = (np.array(ws.model_state_vec) - prior) / SCALE
                    fit_y = np.array(ws.measurement_vec_fit, copy=True)
                    np.testing.assert_allclose(fitted, expected, rtol=0, atol=2e-7)
                    np.testing.assert_allclose(fit_y, K @ expected, rtol=0, atol=2e-7)
                    np.testing.assert_allclose(
                        np.array(ws.measurement_gain_mat) / SCALE[:, None],
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
    fig.suptitle("Workspace LM: covariance representations agree (four overlapping markers per pattern)")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    results = run()
    if "ARTS_HEADLESS" not in os.environ:
        plot(results)
        import matplotlib.pyplot as plt
        plt.show()
