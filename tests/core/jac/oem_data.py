"""Ownership, auxiliary cleanup and reuse through the public OEM workspace API."""
import tempfile
from pathlib import Path

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
ws = pyarts.Workspace()
ws.abs_bands = arts.AbsorptionBands()
ws.subsurf_field = arts.SubsurfaceField()
ws.measurement_sensor = arts.ArrayOfSensorObsel()
ws.surf_field.ellipsoid = [1., 1.]
ws.atm_field[arts.AtmKey.t] = arts.GriddedField3(
    name="Temperature", data=np.array([1., 2.]).reshape(2, 1, 1),
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0., 1.], [0.], [0.]],
)
ws.oemInit()
ws.oemAddTemperature(matrix=np.diag([2., 2.]), inverse=np.diag([.5, .5]))
assert len(ws.oem.covmat_diagonal_blocks) == 1
assert not ws.has("covmat_diagonal_blocks")
with tempfile.TemporaryDirectory() as directory:
    filename = str(Path(directory) / "pending.xml")
    ws.oem.savexml(filename)
    restored = arts.OptimalEstimationData.fromxml(filename)
    assert len(restored.covmat_diagonal_blocks) == 1
# Incomplete numerical inputs must not leave an apparently checked setup.
try:
    ws.oemFinalizeDiagonal()
except RuntimeError as error:
    assert "model_state_vec_apriori" in str(error)
    assert "measurement_vec" in str(error)
    assert "measurement_vec_error_covmat" in str(error)
else:
    raise AssertionError("Incomplete setup was finalized without validation")
assert not ws.oem.checked
# Pending covariance blocks and assembled covariances are owned by oem.
ws.measurement_sensor = arts.ArrayOfSensorObsel([arts.SensorObsel()] * 3)
J = np.array([[1., 0.], [0., 1.], [1., 1.]])
prior = np.array([1., 2.])
y = J @ [3., 4.]
ws.model_state_vecFromData()
ws.measurement_vec = y
ws.measurement_vec_fit = J @ prior
ws.measurement_jac = J
x_storage = np.asarray(ws.model_state_vec).__array_interface__["data"][0]
y_storage = np.asarray(ws.measurement_vec).__array_interface__["data"][0]
j_storage = np.asarray(ws.measurement_jac).__array_interface__["data"][0]
ws.oemSetApriori()
ws.oemSetMeasurement()
ws.oemMeasurementCovmatConstant(value=1.)
ws.oemFinalizeDiagonal()
assert ws.oem.checked
assert len(ws.oem.model_state_covmat.blocks) == 1
ws.oemFinalizeDiagonal()
assert ws.oem.checked
assert len(ws.oem.model_state_covmat.blocks) == 1
# Finalization must revalidate even a previously checked object.
first_measurement = float(ws.oem.measurement_vec[0])
ws.oem.measurement_vec[0] = np.nan
try:
    ws.oemFinalizeDiagonal()
except RuntimeError as error:
    assert "measurement_vec[0]" in str(error)
else:
    raise AssertionError("Invalid observations passed finalization")
assert not ws.oem.checked
ws.oem.measurement_vec[0] = first_measurement
ws.oemFinalizeDiagonal()
assert ws.oem.checked
sa_storage = np.asarray(ws.oem.model_state_covmat.blocks[0].matrix).__array_interface__["data"][0]
sa = np.asarray(ws.oem.model_state_covmat.blocks[0].matrix).copy()
se = ws.oem.measurement_vec_error_covmat.blocks[0].matrix.tocsr().toarray()
# Same-sized setters must preserve both validation and covariance storage.
ws.model_state_vec = prior
ws.measurement_vec = y
x_storage = np.asarray(ws.model_state_vec).__array_interface__["data"][0]
y_storage = np.asarray(ws.measurement_vec).__array_interface__["data"][0]
ws.oemSetApriori()
ws.oemSetMeasurement()
assert ws.oem.checked
assert np.asarray(ws.oem.model_state_covmat.blocks[0].matrix).__array_interface__["data"][0] == sa_storage
assert len(ws.oem.covmat_diagonal_blocks) == 1
np.testing.assert_array_equal(ws.oem.model_state_covmat.blocks[0].matrix, sa)
np.testing.assert_array_equal(ws.oem.measurement_vec_error_covmat.blocks[0].matrix.tocsr().toarray(), se)
assert np.asarray(ws.model_state_vec).size == 0
assert np.asarray(ws.measurement_vec).size == 0
np.testing.assert_array_equal(ws.measurement_vec_fit, J @ prior)
np.testing.assert_array_equal(ws.measurement_jac, J)
assert np.asarray(ws.oem.model_state_vec_apriori).__array_interface__["data"][0] == x_storage
assert np.asarray(ws.oem.measurement_vec).__array_interface__["data"][0] == y_storage
assert np.asarray(ws.measurement_jac).__array_interface__["data"][0] == j_storage
assert np.asarray(ws.oem.measurement_jac).size == 0
assert np.asarray(ws.oem.measurement_vec_fit).size == 0
np.testing.assert_array_equal(ws.oem.model_state_vec_apriori, prior)


# Initialized measurement dimensions constrain subsequent covariance edits.
try:
    ws.oemMeasurementCovmatAdd(matrix=[[1.]])
except RuntimeError:
    pass
else:
    raise AssertionError("Covariance block exceeds initialized measurements")
assert len(ws.oem.measurement_vec_error_covmat.blocks) == 1


evaluations = []

def forward(local):
    evaluations.append(True)
    x = np.asarray(local.get("model_state_vec"))
    local.get("measurement_vec_fit").value = arts.Vector(J @ x)
    local.get("measurement_jac").value = arts.Matrix(
        J if local.get("jac_targets").x_size() else np.empty((0, 0)))


callback = arts.CallbackOperator(forward, ["model_state_vec", "jac_targets"],
                                 ["measurement_vec_fit", "measurement_jac"])
agenda = arts.Agenda("inversion_iterate_agenda")
agenda.add(arts.Method("linear_observation", callback))
agenda.finalize(True)
ws.inversion_iterate_agenda = agenda
expected = prior + np.linalg.solve(.5 * np.eye(2) + J.T @ J, J.T @ (y - J @ prior))
# Data validation must run before covariance preparation,
# forward evaluations or clearing diagnostics/results, for both entry points.
ws.oem.uncheck()
assert not ws.oem.checked
ws.oem.diagnostics.iterations = 123
ws.oem.measurement_vec[0] = np.nan
for calculate in (ws.oemCalc, ws.oemCalcReduced):
    try:
        calculate(settings="gn")
    except RuntimeError as error:
        assert "measurement_vec[0]" in str(error)
        assert not ws.oem.checked
    else:
        raise AssertionError("Invalid OEM data was accepted")
    assert ws.oem.diagnostics.iterations == 123
    assert not evaluations
    assert np.asarray(ws.oem.model_state_vec).size == 0
    assert np.asarray(ws.oem.measurement_vec_fit).size == 0
    assert np.asarray(ws.oem.measurement_jac).size == 0
    np.testing.assert_array_equal(ws.oem.measurement_vec, np.r_[np.nan, y[1:]])
    np.testing.assert_array_equal(ws.oem.model_state_vec_apriori, prior)
ws.oem.measurement_vec[0] = y[0]

# The workspace check is also recordable/executable in an agenda. It changes
# validation status without consuming numerical inputs or running the model.
@pyarts.arts_agenda
def check_oem_setup(ws):
    ws.oemCheck()

check_oem_setup.execute(ws)
assert ws.oem.checked
assert not evaluations
assert ws.oem.diagnostics.iterations == 123
try:
    ws.oemCheck(jac_targets=arts.JacobianTargets())
except RuntimeError as error:
    assert "jac_targets must be finalized" in str(error)
else:
    raise AssertionError("Workspace check accepted unfinalized targets")
assert not ws.oem.checked
np.testing.assert_array_equal(ws.oem.measurement_vec, y)
np.testing.assert_array_equal(ws.oem.model_state_vec_apriori, prior)

ws.oemFinalizeDiagonal()
assert ws.oem.checked
ws.oem.uncheck()
ws.oemCalc(settings=arts.OptimalEstimationSettings(method="gn", stop_dx=1e-10))
assert ws.oem.checked
np.testing.assert_allclose(ws.oem.model_state_vec, expected, atol=1e-10)
ws.oemAveragingKernelCalc()
ws.oemObservationErrorCalc()
ws.oemSmoothingErrorCalc()
ws.oemBasisCalc()
ws.oemBasisReduce(rank=2)
b = np.array(ws.oem.model_state_basis_mat)
c = np.array(ws.oem.measurement_basis_mat)
status = ws.oem.diagnostics.status
assert len(ws.oem.covmat_diagonal_blocks) == 1
ws.oemCheck()
ws.oem.clear_auxiliary()
assert ws.oem.checked
assert len(ws.oem.covmat_diagonal_blocks) == 0
for field in ("measurement_jac", "measurement_vec_fit", "measurement_gain_mat",
              "measurement_averaging_kernel", "observation_error_covmat",
              "smoothing_error_covmat", "basis_singular_values"):
    assert np.asarray(getattr(ws.oem, field)).size == 0, field
assert ws.oem.diagnostics.status == status
np.testing.assert_array_equal(ws.oem.model_state_basis_mat, b)
np.testing.assert_array_equal(ws.oem.measurement_basis_mat, c)
np.testing.assert_allclose(ws.oem.model_state_vec, expected, atol=1e-10)
ws.oem.uncheck()
ws.oemCalcReduced(settings=arts.OptimalEstimationSettings(method="lm", stop_dx=1e-10))
assert ws.oem.checked
np.testing.assert_allclose(ws.oem.model_state_vec, expected, atol=1e-9)

# Both default data and a completed retrieval are printable and XML storable.
with tempfile.TemporaryDirectory() as directory:
    for index, value in enumerate((arts.OptimalEstimationData(), ws.oem)):
        filename = str(Path(directory) / f"oem-{index}.xml")
        value.savexml(filename)
        restored = arts.OptimalEstimationData.fromxml(filename)
        np.testing.assert_array_equal(restored.measurement_vec, value.measurement_vec)
        np.testing.assert_array_equal(restored.measurement_jac, value.measurement_jac)
        assert "OptimalEstimationData" in str(restored)

# Restoring the prior updates only the physical model; prior results are retained.
ws.atm_field[arts.AtmKey.t] = arts.GriddedField3(
    name="Temperature", data=np.array([9., 10.]).reshape(2, 1, 1),
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0., 1.], [0.], [0.]],
)
previous_state = np.array(ws.oem.model_state_vec, copy=True)
previous_jac = np.array(ws.oem.measurement_jac, copy=True)
previous_status = ws.oem.diagnostics.status
ws.oemRestoreApriori()
np.testing.assert_array_equal(ws.oem.model_state_vec, previous_state)
np.testing.assert_array_equal(ws.oem.measurement_jac, previous_jac)
assert ws.oem.diagnostics.status == previous_status
ws.model_state_vecFromData()
np.testing.assert_array_equal(ws.model_state_vec, prior)
np.testing.assert_array_equal(ws.oem.measurement_vec, y)
assert ws.oem.diagnostics.status == previous_status
ws.oemCalc(settings=arts.OptimalEstimationSettings(method="gn", stop_dx=1e-10))
np.testing.assert_allclose(ws.oem.model_state_vec, expected, atol=1e-10)

# New workspaces must never share a mutable OEM default object.
assert not pyarts.Workspace().has("oem")

# Initialization also works without physical fields or a target mapping.
standalone = pyarts.Workspace()
standalone.model_state_vec = prior
standalone.measurement_vec = y
standalone.measurement_vec_fit = J @ prior
standalone.measurement_jac = J
standalone.model_state_covmat = arts.CovarianceMatrix(ws.oem.model_state_covmat)
standalone.measurement_vec_error_covmat = arts.CovarianceMatrix(ws.oem.measurement_vec_error_covmat)
standalone.oemInitFromData()
np.testing.assert_array_equal(standalone.oem.model_state_vec_apriori, prior)
assert np.asarray(standalone.oem.model_state_vec).size == 0
# Initialization leaves modeling results outside OEM. Information analysis
# explicitly supplies its linearization, independently of input initialization.
np.testing.assert_array_equal(standalone.measurement_jac, J)
np.testing.assert_array_equal(standalone.measurement_vec_fit, J @ prior)
assert np.asarray(standalone.oem.measurement_jac).size == 0
assert np.asarray(standalone.oem.measurement_vec_fit).size == 0
standalone.oem.measurement_jac = standalone.measurement_jac
standalone.oemBasisCalc()
standalone.oemBasisReduce(rank=2)
assert standalone.oem.model_state_basis_mat.shape == (2, 2)
report = pyarts.retrieval.information_from_workspace(standalone)
assert report.degrees_of_freedom > 0
assert len(standalone.model_state_covmat.blocks) == 0
assert len(standalone.measurement_vec_error_covmat.blocks) == 0
# Consuming standalone inputs must not consume another OEM object's members.
assert len(ws.oem.model_state_covmat.blocks) == 1
assert len(ws.oem.measurement_vec_error_covmat.blocks) == 1

# Invalid initialization does not take ownership or destroy a previous result.
ws.model_state_covmat = arts.CovarianceMatrix(ws.oem.model_state_covmat)
ws.measurement_vec_error_covmat = arts.CovarianceMatrix(ws.oem.measurement_vec_error_covmat)
ws.model_state_vec = [9.]
ws.measurement_vec = y
ws.measurement_vec_fit = []
ws.measurement_jac = np.empty((0, 0))
try:
    ws.oemInitFromData()
except RuntimeError:
    pass
else:
    raise AssertionError("Incompatible starting state accepted")
assert len(ws.model_state_covmat.blocks) == 1
assert len(ws.measurement_vec_error_covmat.blocks) == 1
np.testing.assert_array_equal(ws.model_state_vec, [9.])
np.testing.assert_array_equal(ws.measurement_vec, y)
np.testing.assert_allclose(ws.oem.model_state_vec, expected, atol=1e-9)
# A wrong measurement covariance must also fail before any input is moved.
ws.model_state_vec = prior
ws.measurement_vec = [0., 0.]
ws.measurement_vec_error_covmatConstant(value=1.)
ws.measurement_vec = y
try:
    ws.oemInitFromData()
except RuntimeError:
    pass
else:
    raise AssertionError("Incompatible measurement covariance accepted")
np.testing.assert_array_equal(ws.model_state_vec, prior)
np.testing.assert_array_equal(ws.measurement_vec, y)
assert len(ws.model_state_covmat.blocks) == 1
assert len(ws.measurement_vec_error_covmat.blocks) == 1
np.testing.assert_allclose(ws.oem.model_state_vec, expected, atol=1e-9)

ws.oem.uncheck()
ws.oem.clear()
assert np.asarray(ws.oem.measurement_vec).size == 0
assert np.asarray(ws.oem.model_state_vec_apriori).size == 0
assert ws.oem.diagnostics.status == arts.OptimalEstimationStatus.NotRun

# Attribute replacement is guarded after checking; reading and in-place value
# edits intentionally retain the existing reference semantics.
checked_data = standalone.oem
checked_data.check(ws.jac_targets)
assert checked_data.checked
try:
    checked_data.check(arts.JacobianTargets())
except RuntimeError as error:
    assert "jac_targets" in str(error)
else:
    raise AssertionError("Unfinalized targets passed check()")
assert not checked_data.checked
checked_data.check()
assert checked_data.checked
for name in (
    "measurement_vec", "model_state_vec_apriori", "model_state_covmat",
    "measurement_vec_error_covmat", "model_state_vec", "measurement_vec_fit",
    "measurement_jac", "model_state_basis_mat", "measurement_basis_mat",
    "model_state_covmat_normalization", "measurement_vec_normalization",
    "measurement_gain_mat", "measurement_averaging_kernel", "observation_error_covmat",
    "smoothing_error_covmat", "diagnostics", "basis_singular_values",
    "basis_lost_dofs", "basis_lost_information_bits", "covmat_diagonal_blocks",
):
    value = getattr(checked_data, name)
    try:
        setattr(checked_data, name, value)
    except RuntimeError as error:
        assert name in str(error) and "uncheck()" in str(error), str(error)
    else:
        raise AssertionError(f"Checked member {name} was replaceable")
try:
    checked_data.checked = False
except AttributeError:
    pass
else:
    raise AssertionError("checked must be read-only")
for operation in (checked_data.clear,):
    try:
        operation()
    except RuntimeError as error:
        assert "uncheck()" in str(error)
    else:
        raise AssertionError("Checked data could be cleared manually")
checked_data.clear_auxiliary()
assert checked_data.checked
assert np.asarray(checked_data.measurement_jac).size == 0
checked_data.check()
checked_data.measurement_vec[0] += 1.
assert checked_data.checked
with tempfile.TemporaryDirectory() as directory:
    filename = str(Path(directory) / "checked.xml")
    checked_data.savexml(filename)
    restored = arts.OptimalEstimationData.fromxml(filename)
    assert not restored.checked
    restored.check()
    assert restored.checked
    try:
        checked_data.readxml(filename)
    except RuntimeError as error:
        assert "uncheck()" in str(error)
    else:
        raise AssertionError("Checked data could be replaced from XML")
checked_data.uncheck()
assert not checked_data.checked
checked_data.measurement_vec = [np.nan]
checked_data.model_state_vec_apriori = []
try:
    ws.oemCheck(oem=checked_data)
except RuntimeError as error:
    report = str(error)
    assert "model_state_vec_apriori" in report and "measurement_vec" in report
    assert "measurement_vec_error_covmat" in report
else:
    raise AssertionError("Invalid data passed check()")
assert not checked_data.checked

# Trusted input setters preserve checked status at the same size and do not
# import unrelated workspace modeling results or reset existing OEM results.
updates = pyarts.Workspace()
updates.oem = restored
updates.oem.check()
old_prior = np.array(updates.oem.model_state_vec_apriori)
old_y = np.array(updates.oem.measurement_vec)
updates.measurement_jac = [[123.]]
updates.measurement_vec_fit = [456.]
old_basis = np.array(updates.oem.model_state_basis_mat)
updates.measurement_vec = old_y + 2.
updates.oemSetMeasurement()
assert updates.oem.checked
assert np.asarray(updates.measurement_vec).size == 0
np.testing.assert_array_equal(updates.oem.measurement_vec, old_y + 2.)
np.testing.assert_array_equal(updates.oem.model_state_vec_apriori, old_prior)
updates.model_state_vec = old_prior + 1.
updates.oemSetApriori()
assert updates.oem.checked
assert np.asarray(updates.model_state_vec).size == 0
np.testing.assert_array_equal(updates.oem.model_state_vec_apriori, old_prior + 1.)
np.testing.assert_array_equal(updates.oem.measurement_vec, old_y + 2.)
np.testing.assert_array_equal(updates.oem.model_state_basis_mat, old_basis)
np.testing.assert_array_equal(updates.measurement_jac, [[123.]])
np.testing.assert_array_equal(updates.measurement_vec_fit, [456.])
assert np.asarray(updates.oem.measurement_jac).size == 0
for source, method, original in (
    ("measurement_vec", updates.oemSetMeasurement, old_y),
    ("model_state_vec", updates.oemSetApriori, old_prior),
):
    setattr(updates, source, [0.])
    method()
    assert not updates.oem.checked
    try:
        updates.oem.check()
    except RuntimeError:
        pass
    else:
        raise AssertionError("Changed dimension passed with old covariance")
    setattr(updates, source, original)
    method()
    assert not updates.oem.checked
    updates.oem.check()
    assert updates.oem.checked

# Full initialization can consume the previous object's own input members.
updates.oemInitFromData(
    model_state_vec=updates.oem.model_state_vec_apriori,
    measurement_vec=updates.oem.measurement_vec,
    model_state_covmat=updates.oem.model_state_covmat,
    measurement_vec_error_covmat=updates.oem.measurement_vec_error_covmat,
)
assert not updates.oem.checked
np.testing.assert_array_equal(updates.oem.model_state_vec_apriori, old_prior)
np.testing.assert_array_equal(updates.oem.measurement_vec, old_y)
updates.oem.check()

ws.oem = restored
assert ws.oem.checked
ws.oemInit()
assert not ws.oem.checked
