"""Named OEM settings and diagnostics preserve numerical results."""

import copy
import pickle

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

# Named diagnostics expose enum status, integer iterations, and owned history/messages.
assert arts.OptimalEstimationStatus(
    "Converged") == arts.OptimalEstimationStatus.Converged
assert set(arts.OptimalEstimationStatus.get_options_as_strings()) == {
    "NotRun", "Converged", "IterationLimit", "DampingLimit", "Error", "StartCostLimit"
}
diagnostics = arts.OptimalEstimationDiagnostics()
assert diagnostics.status == arts.OptimalEstimationStatus.NotRun
diagnostics.status = arts.OptimalEstimationStatus.Converged
diagnostics.initial_cost = 2
diagnostics.final_cost = 1
diagnostics.measurement_cost = 0.5
diagnostics.iterations = 3
assert diagnostics.iterations == 3
assert isinstance(diagnostics.iterations, int)
assert len(diagnostics.lm_ga_history) == 0
assert len(diagnostics.errors) == 0
assert "measurement_cost" in str(diagnostics)
assert copy.copy(diagnostics).final_cost == 1

FIELDS = (
    "initial_damping",
    "decrease_factor",
    "increase_factor",
    "maximum_damping",
    "damping_threshold",
    "convergence_damping_limit",
)
DEFAULTS = [10, 2, 2, 100, 1, 0]
DISTINCT = [12, 3, 2, 1e6, 0.01, 10]
METHODS = ("lm", "ml", "lm_cg", "ml_cg")


def named(values):
    return arts.LevenbergMarquardtSettings(**dict(zip(FIELDS, values)))


def rejects(operation, *names):
    try:
        operation()
    except (RuntimeError, ValueError, TypeError) as error:
        for name in names:
            assert name in str(error), (name, str(error))
    else:
        raise AssertionError("Invalid damping settings were accepted")


def test_value_object():
    defaults = arts.LevenbergMarquardtSettings()
    defaults.validate()
    np.testing.assert_array_equal([getattr(defaults, field)
                                  for field in FIELDS], DEFAULTS)
    for field, value in zip(FIELDS, DEFAULTS):
        assert getattr(defaults, field) == value
        assert f"{field}=" in repr(defaults)
    assert "LevenbergMarquardtSettings" in repr(defaults)
    description = defaults.describe()
    assert repr(defaults) in description
    assert "stop_dx" in description
    assert "covariance" in description.lower()
    assert len(description) > len(repr(defaults))
    rejects(lambda: arts.LevenbergMarquardtSettings(10))  # Constructor is keyword-only.

    for values in (DISTINCT, tuple(DISTINCT), arts.Vector(DISTINCT), np.array(DISTINCT)):
        imported = arts.LevenbergMarquardtSettings(values)
        for field, expected in zip(FIELDS, DISTINCT):
            assert getattr(imported, field) == expected
    for count in (0, 5, 7):
        rejects(lambda: arts.LevenbergMarquardtSettings([1] * count))
    rejects(lambda: arts.LevenbergMarquardtSettings(
        [10, 1, 2, 100, 1, 0]), "decrease_factor")

    settings = named(DISTINCT)
    duplicates = [
        copy.copy(settings),
        copy.deepcopy(settings),
        pickle.loads(pickle.dumps(settings)),
    ]
    settings.initial_damping = 0
    assert settings.initial_damping == 0
    for duplicate in duplicates:
        np.testing.assert_array_equal([getattr(duplicate, field)
                                      for field in FIELDS], DISTINCT)
        assert duplicate is not settings
        for field in FIELDS:
            assert f"{field}=" in repr(duplicate)


def test_validation():
    for field in FIELDS:
        for invalid in (-1, np.nan, np.inf, -np.inf):
            rejects(lambda: arts.LevenbergMarquardtSettings(**{field: invalid}), field)
            changed = arts.LevenbergMarquardtSettings()
            rejects(lambda: setattr(changed, field, invalid), field)
            # A failed edit leaves a usable configuration and its diagnostics
            # identify the field before implicit conversion can obscure them.
            changed.validate()
            np.testing.assert_array_equal(
                [getattr(changed, field) for field in FIELDS], DEFAULTS)

    for field, invalid in (
        ("decrease_factor", 0),
        ("decrease_factor", 1),
        ("increase_factor", 0),
        ("increase_factor", 1),
        ("maximum_damping", 0),
        ("damping_threshold", 0),
    ):
        rejects(lambda: arts.LevenbergMarquardtSettings(**{field: invalid}), field)
        changed = arts.LevenbergMarquardtSettings()
        rejects(lambda: setattr(changed, field, invalid), field)
        np.testing.assert_array_equal([getattr(changed, field)
                                      for field in FIELDS], DEFAULTS)

    for field in ("initial_damping", "damping_threshold"):
        rejects(
            lambda: arts.LevenbergMarquardtSettings(**{field: 101}),
            field,
            "maximum_damping",
        )
        changed = arts.LevenbergMarquardtSettings()
        rejects(lambda: setattr(changed, field, 101), field, "maximum_damping")
        np.testing.assert_array_equal([getattr(changed, field)
                                      for field in FIELDS], DEFAULTS)

    # These are useful supported boundaries, not additional ordering constraints.
    for overrides in (
        {"initial_damping": 0},
        {"initial_damping": 0.5},  # Initial damping can be below its threshold.
        {"initial_damping": 100, "damping_threshold": 100},
        {"initial_damping": 0, "damping_threshold": 0.1, "maximum_damping": 0.1},
        {"convergence_damping_limit": 101},  # May exceed maximum damping.
    ):
        settings = arts.LevenbergMarquardtSettings(**overrides)
        settings.validate()

    # Coupled settings can be replaced together or edited in a valid order.
    changed = arts.LevenbergMarquardtSettings()
    changed.maximum_damping = 200
    changed.initial_damping = 200
    replacement = arts.LevenbergMarquardtSettings(
        initial_damping=200, maximum_damping=200)
    np.testing.assert_array_equal([getattr(changed, field) for field in FIELDS], [
                                  getattr(replacement, field) for field in FIELDS])


def covariance(values):
    matrix = arts.Matrix(values)
    size = np.asarray(matrix).shape[0]
    result = arts.CovarianceMatrix()
    result.blocks = [
        arts.Block(arts.Range(0, size), arts.Range(0, size), (0, 0), matrix)
    ]
    return result


def workspace(nonlinear=False):
    ws = pyarts.Workspace()
    for variable, group in (
        ("model_state_vec", arts.Vector),
        ("measurement_vec_fit", arts.Vector),
        ("measurement_jac", arts.Matrix),
        ("atm_field", arts.AtmField),
        ("abs_bands", arts.AbsorptionBands),
        ("measurement_sensor", arts.ArrayOfSensorObsel),
        ("surf_field", arts.SurfaceField),
        ("subsurf_field", arts.SubsurfaceField),
        ("jac_targets", arts.JacobianTargets),
    ):
        setattr(ws, variable, group())

    # Only target dimensions are used by the synthetic forward model. Creating
    # tiny temperature grids satisfies the public agenda's target-size checks.
    size = 1 if nonlinear else 2
    ws.surf_field.ellipsoid = [1, 1]
    ws.atm_field[arts.AtmKey.t] = arts.GriddedField3(
        name="Synthetic state metadata",
        data=np.ones((size, 1, 1)),
        grid_names=["Altitude", "Latitude", "Longitude"],
        grids=[np.arange(size, dtype=float), [0], [0]],
    )
    ws.jac_targetsAddTemperature()
    ws.jac_targetsFinalize()

    if nonlinear:
        ws.model_state_vec = [0.1]
        ws.model_state_vec_apriori = [1]
        ws.measurement_vec = [4]
        ws.model_state_covmat = covariance([[4]])
        ws.measurement_vec_error_covmat = covariance([[0.25]])
    else:
        ws.model_state_vec_apriori = [0.5, -0.25]
        ws.measurement_vec = [2, -1, 1.5]
        ws.model_state_covmat = covariance([[4, 1], [1, 2]])
        ws.measurement_vec_error_covmat = covariance(
            [[1, 0.2, 0], [0.2, 2, 0.3], [0, 0.3, 0.5]]
        )
    jacobian = np.array([[1, 2], [2, -1], [1, 1]], dtype=float)

    def forward(local):
        state = np.asarray(local.get("model_state_vec"))
        if nonlinear:
            fit = np.array([state[0] ** 2])
            derivative = np.array([[2 * state[0]]])
        else:
            fit = jacobian @ state + [0.25, -0.5, 1]
            derivative = jacobian
        if not local.get("jac_targets").x_size():
            derivative = np.empty((0, 0))
        # Preserve the existing output objects shared with the outer workspace.
        local.get("measurement_vec_fit").value = arts.Vector(fit)
        local.get("measurement_jac").value = arts.Matrix(derivative)

    callback = arts.CallbackOperator(
        forward,
        ["model_state_vec", "jac_targets"],
        ["measurement_vec_fit", "measurement_jac"],
    )
    agenda = arts.Agenda("inversion_iterate_agenda")
    agenda.add(arts.Method("synthetic_forward_model", callback))
    agenda.finalize(True)
    ws.inversion_iterate_agenda = agenda
    return ws


def retrieve(method, settings, nonlinear=False, max_iter=40, stop_dx=1e-9):
    ws = workspace(nonlinear)
    ws.OEM(
        method=method,
        **({"lm_ga_settings": settings} if settings is not None else {}),
        max_iter=max_iter,
        stop_dx=stop_dx,
        display_progress=0,
    )
    assert len(ws.oem_diagnostics.errors) == 0, str(ws.oem_diagnostics.errors)
    return {
        key: copy.deepcopy(ws.oem_diagnostics) if key == "oem_diagnostics" else np.array(
            ws.get(key), copy=True)
        for key in (
            "model_state_vec",
            "measurement_vec_fit",
            "measurement_jac",
            "measurement_gain_mat",
            "oem_diagnostics",
        )
    }


def equivalent(left, right):
    for output in left:
        if output == "oem_diagnostics":
            a, b = left[output], right[output]
            assert a.status == b.status
            assert a.iterations == b.iterations
            assert list(a.errors) == list(b.errors)
            for field in ("initial_cost", "final_cost", "measurement_cost", "lm_ga_history"):
                np.testing.assert_allclose(getattr(a, field), getattr(
                    b, field), rtol=0, atol=1e-12, equal_nan=True)
            continue
        np.testing.assert_allclose(
            left[output],
            right[output],
            rtol=0,
            atol=1e-12,
            equal_nan=True,
            err_msg=output,
        )


def test_retrieval_equivalence():
    for method in METHODS:
        equivalent(retrieve(method, None), retrieve(
            method, arts.LevenbergMarquardtSettings()))
        # Correlated prior distinguishes diag(Sa^-1) from inverse(diag(Sa)).
        values = [10, 3, 2, 1e8, 0.1, 0]
        result = retrieve(method, named(values))
        equivalent(result, retrieve(method, copy.deepcopy(named(values))))
        equivalent(result, retrieve(method, values))
        equivalent(result, retrieve(method, tuple(values)))
        equivalent(result, retrieve(method, arts.Vector(values)))

        np.testing.assert_allclose(
            result["model_state_vec"],
            [3112 / 18575, 21727 / 37150],
            rtol=0,
            atol=2e-7,
        )
        np.testing.assert_allclose(
            result["measurement_gain_mat"],
            [
                [1276 / 18575, 1092 / 3715, 4586 / 18575],
                [4323 / 18575, -784 / 3715, 4328 / 18575],
            ],
            rtol=0,
            atol=1e-11,
        )
        np.testing.assert_allclose(
            result["oem_diagnostics"].final_cost, 245309 / 891600, rtol=0, atol=1e-10
        )
        assert result["oem_diagnostics"].status == arts.OptimalEstimationStatus.Converged

        # All six distinct named controls reach both the direct and CG solver.
        result = retrieve(method, named(DISTINCT), max_iter=1, stop_dx=1e3)
        equivalent(
            result, retrieve(method, copy.deepcopy(
                named(DISTINCT)), max_iter=1, stop_dx=1e3)
        )
        np.testing.assert_allclose(
            result["model_state_vec"],
            [31732 / 69551, 24127 / 139102],
            rtol=0,
            atol=1e-11,
        )
        np.testing.assert_array_equal(result["oem_diagnostics"].lm_ga_history, [12, 4])
        assert result["oem_diagnostics"].status == arts.OptimalEstimationStatus.Converged

        # Seven rejected quadratic-model trials restart/increase damping before
        # gamma=6.4 accepts x=1.8. This independently tests nonlinear behavior.
        values = [0, 3, 2, 100, 0.1, 0]
        result = retrieve(method, named(values), nonlinear=True, max_iter=1)
        equivalent(
            result, retrieve(method, copy.deepcopy(
                named(values)), nonlinear=True, max_iter=1)
        )
        np.testing.assert_allclose(result["model_state_vec"], [1.8], rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            result["measurement_vec_fit"], [3.24], rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            result["measurement_jac"], [[3.6]], rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            result["oem_diagnostics"].lm_ga_history, [0, 6.4], rtol=0, atol=1e-12
        )
        assert result["oem_diagnostics"].status == arts.OptimalEstimationStatus.IterationLimit


def test_agenda_capture():
    # arts_agenda resolves captures in the defining module's globals, then
    # constructs the named settings workspace type.
    global captured_settings
    captured_settings = named(DISTINCT)

    @pyarts.arts_agenda
    def captured_retrieval(ws):
        ws.OEM(method="lm", lm_ga_settings=captured_settings)

    captured = [
        method.val
        for method in captured_retrieval.methods
        if method.name == "@lm_ga_settings"
    ]
    assert len(captured) == 1
    assert isinstance(captured[0], arts.LevenbergMarquardtSettings)
    np.testing.assert_array_equal([getattr(captured[0], field)
                                  for field in FIELDS], DISTINCT)
    captured_settings.initial_damping = 0
    np.testing.assert_array_equal([getattr(captured[0], field)
                                  for field in FIELDS], DISTINCT)

    captured_settings = DISTINCT.copy()

    @pyarts.arts_agenda
    def shorthand_retrieval(ws):
        ws.OEM(method="lm", lm_ga_settings=captured_settings)
    shorthand = next(
        method.val for method in shorthand_retrieval.methods if method.name == "@lm_ga_settings")
    assert isinstance(shorthand, arts.LevenbergMarquardtSettings)
    for field, expected in zip(FIELDS, DISTINCT):
        assert getattr(shorthand, field) == expected


def test_damping_outcomes():
    for method in METHODS:
        # These named values used to produce false convergence because the
        # failure sentinel maximum_damping + 1 rounds back to maximum_damping.
        result = retrieve(
            method,
            arts.LevenbergMarquardtSettings(initial_damping=1e20, maximum_damping=1e20),
        )
        assert result["oem_diagnostics"].status == arts.OptimalEstimationStatus.DampingLimit, result["oem_diagnostics"]
        np.testing.assert_array_equal(result["model_state_vec"], [0.5, -0.25])
        np.testing.assert_array_equal(
            result["oem_diagnostics"].lm_ga_history[:2], [1e20, 1e20])
        assert result["oem_diagnostics"].final_cost == result["oem_diagnostics"].initial_cost

        # An accurate affine solution must not turn into damping exhaustion
        # when reduction ratios become dominated by cost-evaluation roundoff.
        for tolerance in (1e-12, 1e-16, 1e-20):
            result = retrieve(method, named([10, 3, 2, 1e8, 0.1, 0]), stop_dx=tolerance)
            assert result["oem_diagnostics"].status == arts.OptimalEstimationStatus.Converged, result["oem_diagnostics"]
            np.testing.assert_allclose(
                result["model_state_vec"],
                [3112 / 18575, 21727 / 37150],
                rtol=0,
                atol=2e-7,
            )

        # The nonlinear fixture really rejects every permitted step. The fix
        # must still report exhaustion and preserve the last accepted state.
        result = retrieve(method, named([0, 3, 2, 1, 0.1, 0]), nonlinear=True)
        assert result["oem_diagnostics"].status == arts.OptimalEstimationStatus.DampingLimit, result["oem_diagnostics"]
        np.testing.assert_array_equal(result["model_state_vec"], [0.1])
        np.testing.assert_allclose(
            result["measurement_vec_fit"], [0.01], rtol=0, atol=1e-16
        )


if __name__ == "__main__":
    test_value_object()
    test_validation()
    test_retrieval_equivalence()
    test_agenda_capture()
    test_damping_outcomes()
