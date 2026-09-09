"""Matching-grid correlations and recovery of a prior from a perturbed start."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def field(values, altitude=(0, 4000, 12000)):
    return arts.GriddedField3(
        data=np.asarray(values).reshape(3, 1, 1),
        grid_names=["Altitude", "Latitude", "Longitude"],
        grids=[list(altitude), [0], [0]],
    )


def setup(sparse=False, third=False):
    ws = pyarts.workspace.Workspace()
    ws.atm_fieldInit(toa=12000.0)
    ws.surf_fieldPlanet(option="Earth")
    ws.surf_field["t"] = 295.0
    ws.atm_field["t"] = field([290.0, 265.0, 225.0])
    ws.atm_field["p"] = field([100000.0, 60000.0, 20000.0])
    ws.atm_field["H2O"] = field([0.012, 0.003, 0.0002])
    ws.atm_field["O2"] = field([0.2095, 0.2095, 0.2095])
    ws.abs_bands = arts.AbsorptionBands()
    ws.measurement_sensor = arts.ArrayOfSensorObsel()

    ws.RetrievalInit()

    def matrix(variance):
        value = np.diag(np.full(3, variance))
        return arts.Sparse(value) if sparse else value

    # PWR98 uses numerical derivatives. Perturbations are in physical
    # kelvin/VMR units even when the retrieval coordinate is logarithmic.
    ws.RetrievalAddTemperature(matrix=matrix(9.0), d=1e-3)
    ws.RetrievalAddSpeciesVMR(species="H2O", matrix=matrix(0.04), d=1e-7)
    if third:
        ws.RetrievalAddSpeciesVMR(species="O2", matrix=matrix(0.01))
        ws.RetrievalAddPressure(matrix=matrix(10000.0))
    ws.RetrievalFinalizeDiagonal()
    ws.jac_targetsToggleLogarithmicAtmTarget(key="H2O")
    return ws


def blocks(ws):
    return {
        tuple(b.indices): np.asarray(
            b.matrix.tocsr().toarray()
            if isinstance(b.matrix, arts.Sparse)
            else b.matrix
        ).copy()
        for b in ws.model_state_covmat.blocks
    }


def correlate(ws, rho, first="temperature", second="H2O"):
    ws.model_state_covmatCorrelate(target1=first, target2=second, correlation=rho)


def fails_unchanged(ws, rho, first="temperature", second="H2O"):
    before = blocks(ws)
    try:
        correlate(ws, rho, first, second)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Invalid correlation accepted")
    after = blocks(ws)
    assert before.keys() == after.keys()
    for key in before:
        np.testing.assert_array_equal(before[key], after[key])


for sparse in (False, True):
    ws = setup(sparse)
    marginal = blocks(ws)
    for rho in (0.6, -0.4):
        correlate(ws, rho)
        result = blocks(ws)
        for key in marginal:
            np.testing.assert_array_equal(result[key], marginal[key])
        np.testing.assert_allclose(result[(0, 1)], np.eye(3) * rho * 0.6)
        # Reversed arguments must replace the same upper-triangular block.
        correlate(ws, rho, "H2O", "temperature")
        assert len(blocks(ws)) == 3
    for rho in (-1.0, 1.0, np.nan, np.inf):
        fails_unchanged(ws, rho)
    fails_unchanged(ws, 0.2, "temperature", "temperature")
    fails_unchanged(ws, 0.2, "temperature", "O2")
    correlate(ws, 0.0)
    assert blocks(ws).keys() == marginal.keys()
    ws.atm_field["H2O"] = field([0.012, 0.003, 0.0002], (0, 5000, 12000))
    fails_unchanged(ws, 0.2)

# Individually valid pair coefficients can form an invalid joint covariance.
ws = setup(third=True)
correlate(ws, 0.9)
fails_unchanged(ws, 0.9, "temperature", "O2")
# Also exercise the two same-enum overloads with distinct targets.
correlate(ws, 0.1, "temperature", "pressure")
correlate(ws, 0.1, "H2O", "O2")
np.testing.assert_allclose(blocks(ws)[(0, 3)], np.eye(3) * 30.0)
np.testing.assert_allclose(blocks(ws)[(1, 2)], np.eye(3) * 0.002)

# Exercise the generated workspace dispatch (Python direct calls bypass it).
agenda = arts.Agenda("correlation_setup")
for name, value in (
    ("cross_first", arts.AtmKey.temperature),
    ("cross_second", arts.SpeciesEnum.H2O),
    ("cross_rho", arts.Numeric(0.5)),
):
    agenda.add(arts.Method(name, value))
agenda.add(
    arts.Method(
        "model_state_covmatCorrelate",
        [],
        {
            "target1": "cross_first",
            "target2": "cross_second",
            "correlation": "cross_rho",
        },
    )
)
ws = setup()
agenda.execute(ws)
np.testing.assert_allclose(blocks(ws)[(0, 1)], np.eye(3) * 0.3)

# Off-diagonal marginal entries are outside the helper's contract.
for sparse in (False, True):
    ws = setup(sparse)
    entries = ws.model_state_covmat.blocks
    dense = np.eye(3) * 9.0
    dense[0, 1] = dense[1, 0] = 0.1
    entries[0].matrix = arts.Sparse(dense) if sparse else arts.Matrix(dense)
    ws.model_state_covmat.blocks = entries
    fails_unchanged(ws, 0.2)

# A full radiative-transfer retrieval in T and ln(H2O) coordinates.
ws = setup()
correlate(ws, 0.6)
ws.freq_grid = np.linspace(20e9, 200e9, 61)
ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98"])
ws.ReadCatalogData()
ws.spectral_propmat_agendaAuto()
ws.spectral_rad_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()
ws.measurement_sensorSimple(pos=[0.0, 0.0, 0.0], los=[0.0, 0.0])
ws.model_state_vec_aprioriFromData()
prior = np.array(ws.model_state_vec_apriori)
ws.measurement_vecFromSensor()
observations = np.array(ws.measurement_vec)
ws.measurement_vec_error_covmatConstant(value=0.01)

# Change the starting atmosphere, keeping both prior and observations fixed.
ws.atm_field["t"] = field([294.0, 270.0, 231.0])
ws.atm_field["H2O"] = field(np.array([0.012, 0.003, 0.0002]) * 1.3)
ws.model_state_vecFromData()
manipulated = np.array(ws.model_state_vec, copy=True)
assert np.linalg.norm(manipulated - prior) > 1
ws.measurement_vecFromSensor()
assert np.linalg.norm(np.array(ws.measurement_vec) - observations) > 0.1
ws.measurement_vec = observations
ws.measurement_vec_fit = []
ws.measurement_jac = arts.Matrix()
ws.OEM(method="lm", lm_ga_settings=arts.OEMLMSettings(), max_iter=100, stop_dx=1e-12)
assert ws.oem_diagnostics[0] == 0, ws.oem_diagnostics
np.testing.assert_allclose(
    (np.array(ws.model_state_vec) - prior) / np.array([3.0] * 3 + [0.2] * 3),
    0.0,
    atol=1e-5,
)
np.testing.assert_allclose(ws.measurement_vec_fit, observations, atol=1e-5)
np.testing.assert_allclose(
    ws.atm_field["t"].data.data.value.reshape(-1), [290.0, 265.0, 225.0], atol=1e-5
)
np.testing.assert_allclose(
    ws.atm_field["H2O"].data.data.value.reshape(-1), [0.012, 0.003, 0.0002], rtol=1e-5
)


def snapshot():
    with TemporaryDirectory() as directory:
        path = Path(directory) / "covariance.xml"
        ws.model_state_covmat.savexml(str(path))
        return path.read_text()


# Updating a covariance after OEM must discard any cached inverse.
before = snapshot()
assert len(ET.fromstring(before).find("CovarianceMatrix")[1]) > 0
fails_unchanged(ws, 1.0)
assert snapshot() == before
correlate(ws, 0.3)
ws.model_state_covmat.validate()
assert len(ET.fromstring(snapshot()).find("CovarianceMatrix")[1]) == 0
print("Correlated T/log-water retrieval recovered the prior:", ws.oem_diagnostics)


if "ARTS_HEADLESS" not in os.environ:
    import matplotlib.pyplot as plt

    fitted = np.array(ws.model_state_vec)
    altitude = np.asarray(ws.atm_field["t"].data.grids[0]) / 1000
    nlevels = len(altitude)
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(9, 5))
    for ax, elements, label in zip(
        axes,
        (slice(0, nlevels), slice(nlevels, 2 * nlevels)),
        ("Temperature [K]", "ln(H2O VMR)"),
    ):
        ax.plot(manipulated[elements], altitude, "o--", label="Manipulated")
        ax.plot(prior[elements], altitude, "o-", label="A priori")
        # Open markers keep the fitted state visible where it overlaps the prior.
        ax.plot(fitted[elements], altitude, "s:", fillstyle="none", label="Fitted")
        ax.set_xlabel(label)
        ax.grid(True)
        ax.legend()
    axes[0].set_ylabel("Altitude [km]")
    fig.suptitle("Temperature / log-water retrieval (prior correlation = 0.6)")
    fig.tight_layout()
    plt.show()

pyarts.retrieval.information_from_workspace(ws)
