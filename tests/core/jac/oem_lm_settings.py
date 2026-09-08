"""Named OEM damping controls preserve the legacy numerical interface."""

import copy
import pickle

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
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
    return arts.OEMLMSettings(**dict(zip(FIELDS, values)))


def rejects(operation, *names):
    try:
        operation()
    except (RuntimeError, ValueError, TypeError) as error:
        for name in names:
            assert name in str(error), (name, str(error))
    else:
        raise AssertionError("Invalid damping settings were accepted")


def test_value_object():
    defaults = arts.OEMLMSettings()
    defaults.validate()
    np.testing.assert_array_equal(defaults.as_vector(), DEFAULTS)
    for field, value in zip(FIELDS, DEFAULTS):
        assert getattr(defaults, field) == value
        assert f"{field}=" in repr(defaults)
    assert "OEMLMSettings" in repr(defaults)
    description = defaults.describe()
    assert repr(defaults) in description
    assert "stop_dx" in description
    assert "covariance" in description.lower()
    assert len(description) > len(repr(defaults))
    rejects(lambda: arts.OEMLMSettings(10))  # Constructor is keyword-only.

    settings = named(DISTINCT)
    vector = settings.as_vector()
    assert isinstance(vector, arts.Vector)
    np.testing.assert_array_equal(vector, DISTINCT)
    np.testing.assert_array_equal(arts.Vector(settings), DISTINCT)
    imported = arts.OEMLMSettings.from_vector(vector)
    np.testing.assert_array_equal(imported.as_vector(), DISTINCT)
    np.testing.assert_array_equal(
        arts.OEMLMSettings.from_vector(DISTINCT).as_vector(), DISTINCT
    )

    # Conversion and serialization produce independent values.
    vector[0] = 13
    assert settings.initial_damping == 12
    assert imported.initial_damping == 12
    duplicates = [
        copy.copy(settings),
        copy.deepcopy(settings),
        pickle.loads(pickle.dumps(settings)),
    ]
    settings.initial_damping = 0
    assert settings.as_vector()[0] == 0
    for duplicate in duplicates:
        np.testing.assert_array_equal(duplicate.as_vector(), DISTINCT)
        assert duplicate is not settings
        for field in FIELDS:
            assert f"{field}=" in repr(duplicate)


def test_validation():
    for field in FIELDS:
        for invalid in (-1, np.nan, np.inf, -np.inf):
            rejects(lambda: arts.OEMLMSettings(**{field: invalid}), field)
            changed = arts.OEMLMSettings()
            rejects(lambda: setattr(changed, field, invalid), field)
            # A failed edit leaves a usable configuration and its diagnostics
            # identify the field before implicit conversion can obscure them.
            changed.validate()
            np.testing.assert_array_equal(changed.as_vector(), DEFAULTS)
            np.testing.assert_array_equal(arts.Vector(changed), DEFAULTS)

    for field, invalid in (
        ("decrease_factor", 0),
        ("decrease_factor", 1),
        ("increase_factor", 0),
        ("increase_factor", 1),
        ("maximum_damping", 0),
        ("damping_threshold", 0),
    ):
        rejects(lambda: arts.OEMLMSettings(**{field: invalid}), field)
        changed = arts.OEMLMSettings()
        rejects(lambda: setattr(changed, field, invalid), field)
        np.testing.assert_array_equal(changed.as_vector(), DEFAULTS)

    for field in ("initial_damping", "damping_threshold"):
        rejects(
            lambda: arts.OEMLMSettings(**{field: 101}),
            field,
            "maximum_damping",
        )
        changed = arts.OEMLMSettings()
        rejects(lambda: setattr(changed, field, 101), field, "maximum_damping")
        np.testing.assert_array_equal(changed.as_vector(), DEFAULTS)

    for count in (0, 5, 7):
        rejects(
            lambda: arts.OEMLMSettings.from_vector(arts.Vector([1] * count)),
            "6",
        )
    for i, field in enumerate(FIELDS):
        values = DISTINCT.copy()
        values[i] = np.nan
        rejects(lambda: arts.OEMLMSettings.from_vector(arts.Vector(values)), field)

    # These are useful supported boundaries, not additional ordering constraints.
    for overrides in (
        {"initial_damping": 0},
        {"initial_damping": 0.5},  # Initial damping can be below its threshold.
        {"initial_damping": 100, "damping_threshold": 100},
        {"initial_damping": 0, "damping_threshold": 0.1, "maximum_damping": 0.1},
        {"convergence_damping_limit": 101},  # May exceed maximum damping.
    ):
        settings = arts.OEMLMSettings(**overrides)
        settings.validate()
        np.testing.assert_array_equal(
            arts.OEMLMSettings.from_vector(settings.as_vector()).as_vector(),
            settings.as_vector(),
        )

    # Coupled settings can be replaced together or edited in a valid order.
    changed = arts.OEMLMSettings()
    changed.maximum_damping = 200
    changed.initial_damping = 200
    replacement = arts.OEMLMSettings(initial_damping=200, maximum_damping=200)
    np.testing.assert_array_equal(changed.as_vector(), replacement.as_vector())


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
        if not int(local.get("do_jac")):
            derivative = np.empty((0, 0))
        # Preserve the existing output objects shared with the outer workspace.
        local.get("measurement_vec_fit").value = arts.Vector(fit)
        local.get("measurement_jac").value = arts.Matrix(derivative)

    callback = arts.CallbackOperator(
        forward,
        ["model_state_vec", "do_jac"],
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
        lm_ga_settings=settings,
        max_iter=max_iter,
        stop_dx=stop_dx,
        display_progress=0,
    )
    assert len(ws.errors) == 0, str(ws.errors)
    return {
        key: np.array(ws.get(key), copy=True)
        for key in (
            "model_state_vec",
            "measurement_vec_fit",
            "measurement_jac",
            "measurement_gain_mat",
            "oem_diagnostics",
            "lm_ga_history",
        )
    }


def equivalent(left, right):
    for output in left:
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
        # Correlated prior distinguishes diag(Sa^-1) from inverse(diag(Sa)).
        values = [10, 3, 2, 1e8, 0.1, 0]
        result = retrieve(method, named(values))
        equivalent(result, retrieve(method, arts.Vector(values)))
        equivalent(result, retrieve(method, values))  # Legacy lists remain valid.
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
            result["oem_diagnostics"][2], 245309 / 891600, rtol=0, atol=1e-10
        )
        assert result["oem_diagnostics"][0] == 0

        # All six distinct named controls reach both the direct and CG solver.
        result = retrieve(method, named(DISTINCT), max_iter=1, stop_dx=1e3)
        equivalent(
            result, retrieve(method, arts.Vector(DISTINCT), max_iter=1, stop_dx=1e3)
        )
        np.testing.assert_allclose(
            result["model_state_vec"],
            [31732 / 69551, 24127 / 139102],
            rtol=0,
            atol=1e-11,
        )
        np.testing.assert_array_equal(result["lm_ga_history"], [12, 4])
        assert result["oem_diagnostics"][0] == 0

        # Seven rejected quadratic-model trials restart/increase damping before
        # gamma=6.4 accepts x=1.8. This independently tests nonlinear behavior.
        values = [0, 3, 2, 100, 0.1, 0]
        result = retrieve(method, named(values), nonlinear=True, max_iter=1)
        equivalent(
            result, retrieve(method, arts.Vector(values), nonlinear=True, max_iter=1)
        )
        np.testing.assert_allclose(result["model_state_vec"], [1.8], rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            result["measurement_vec_fit"], [3.24], rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            result["measurement_jac"], [[3.6]], rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            result["lm_ga_history"], [0, 6.4], rtol=0, atol=1e-12
        )
        assert result["oem_diagnostics"][0] == 1


def test_agenda_capture():
    # arts_agenda resolves captures in the defining module's globals, then
    # constructs a Vector for the existing workspace argument type.
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
    assert isinstance(captured[0], arts.Vector)
    np.testing.assert_array_equal(captured[0], DISTINCT)
    captured_settings.initial_damping = 0
    np.testing.assert_array_equal(captured[0], DISTINCT)


if __name__ == "__main__":
    test_value_object()
    test_validation()
    test_retrieval_equivalence()
    test_agenda_capture()
