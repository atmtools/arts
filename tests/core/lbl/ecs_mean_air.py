"""Check the legacy mean-air coefficient average without catalogue data.

This tests its defined coefficient averaging, not equivalence to mixing full
collision matrices or spectra.
"""

import math

import pyarts3 as pyarts


arts = pyarts.arts
O2 = arts.SpeciesEnum("O2")
N2 = arts.SpeciesEnum("N2")
HE = arts.SpeciesEnum("He")
AIR = arts.SpeciesEnum("AIR")
FIELDS = ("scaling", "beta", "lambda_", "collisional_distance")


def preset_workspace():
    ws = pyarts.Workspace()
    assert callable(ws.abs_ecs_dataAddMeanAir)
    ws.abs_ecs_dataInit()
    ws.abs_ecs_dataAddMakarov2020()
    ws.abs_ecs_dataAddRodrigues1997()
    ws.abs_ecs_dataAddTran2011()
    assert len(ws.abs_ecs_data) == 4
    return ws


def snapshot(ws, include_bath=True):
    """Stable exact comparison, including NaN coefficients in invalid inputs."""
    return sorted(
        (
            str(isot),
            str(partner),
            field,
            str(getattr(params, field).type),
            tuple(float(x).hex() for x in getattr(params, field).data),
        )
        for isot, partners in ws.abs_ecs_data.items()
        for partner, params in partners.items()
        if include_bath or partner != AIR
        for field in FIELDS
    )


# The 0.21 O2 and 0.79 N2 mean of the preset coefficients, confirmed with
# independent code.  The CO2 isotopologues share their coefficients.
CO2_MEAN = {
    "scaling": ("T1", [5251.138953450777, 0.7765]),
    "beta": ("T0", [0.00779]),
    "lambda_": ("T1", [0.8121, -0.007101999999999997]),
    "collisional_distance": ("T0", [2.242e-10]),
}
O2_N2_MEAN = {
    "CO2-626": CO2_MEAN,
    "CO2-628": CO2_MEAN,
    "CO2-636": CO2_MEAN,
    "O2-66": {
        "scaling": ("T0", [1.0]),
        "beta": ("T0", [0.567]),
        "lambda_": ("T0", [0.39]),
        "collisional_distance": ("T0", [6.1e-11]),
    },
}


def coefficients_of(ws, partner):
    """The coefficients of a partner, for an average that only it contributes to"""
    return {
        str(isot): {
            field: (
                str(getattr(partners[partner], field).type),
                [float(x) for x in getattr(partners[partner], field).data],
            )
            for field in FIELDS
        }
        for isot, partners in ws.abs_ecs_data.items()
    }


def check_average(ws, weights, species, expected):
    unchanged = snapshot(ws, include_bath=False)
    ws.abs_ecs_dataAddMeanAir(vmrs=weights, species=species)
    assert snapshot(ws, include_bath=False) == unchanged
    assert sorted(str(isot) for isot in ws.abs_ecs_data) == sorted(expected)
    for isot, partners in ws.abs_ecs_data.items():
        for field in FIELDS:
            actual = getattr(partners[AIR], field)
            model_type, coefficients = expected[str(isot)][field]
            assert str(actual.type) == model_type
            assert len(actual.data) == len(coefficients)
            for actual_value, expected_value in zip(actual.data, coefficients):
                assert math.isclose(
                    float(actual_value), expected_value, rel_tol=2e-14, abs_tol=1e-30
                ), (str(isot), field, actual_value, expected_value)


def fails_without_mutation(ws, weights, species):
    before = snapshot(ws)
    try:
        ws.abs_ecs_dataAddMeanAir(vmrs=weights, species=species)
    except Exception:
        pass
    else:
        raise AssertionError(("Invalid mean-air input accepted", weights, species))
    assert snapshot(ws) == before, "A rejected mean-air update modified ECS data"


def rollback_workspace():
    ws = preset_workspace()
    check_average(ws, [0.21, 0.79], [O2, N2], O2_N2_MEAN)
    # Distinct existing bath values expose partial updates before a later error.
    for partners in ws.abs_ecs_data.values():
        partners[AIR].scaling = arts.TemperatureModel("T0", [321.0])
    return ws, list(ws.abs_ecs_data)[-1]


ws = preset_workspace()
check_average(ws, [0.21, 0.79], [O2, N2], O2_N2_MEAN)
# Accepted roundoff is normalized away.
check_average(ws, [1.00001], [O2], coefficients_of(ws, O2))
# No data are needed for an absent partner.
check_average(ws, [1.0, 0.0], [O2, HE], coefficients_of(ws, O2))

for weights, species in (
    ([], []),
    ([1.0], [O2, N2]),
    ([-0.1, 1.1], [O2, N2]),
    ([math.nan, 1.0], [O2, N2]),
    ([math.inf, 1.0], [O2, N2]),
    ([0.0, 0.0], [O2, N2]),
    ([0.2, 0.2], [O2, N2]),
    ([1.0], [AIR]),
    ([0.5, 0.5], [O2, HE]),
):
    fails_without_mutation(ws, weights, species)

for field in FIELDS:
    ws, isot = rollback_workspace()
    oxygen = getattr(ws.abs_ecs_data[isot][O2], field)
    other = (
        arts.TemperatureModel("T1", [1.0, 0.5])
        if str(oxygen.type) == "T0"
        else arts.TemperatureModel("T0", [1.0])
    )
    setattr(ws.abs_ecs_data[isot][N2], field, other)
    fails_without_mutation(ws, [0.21, 0.79], [O2, N2])

for values in ([1.0, 2.0], []):
    ws, isot = rollback_workspace()
    ws.abs_ecs_data[isot][O2].scaling = arts.TemperatureModel("POLY", [1.0])
    ws.abs_ecs_data[isot][N2].scaling = arts.TemperatureModel("POLY", values)
    fails_without_mutation(ws, [0.21, 0.79], [O2, N2])

for invalid in (math.nan, math.inf):
    ws, isot = rollback_workspace()
    model = ws.abs_ecs_data[isot][N2].beta
    coefficients = [float(x) for x in model.data]
    coefficients[0] = invalid
    ws.abs_ecs_data[isot][N2].beta = arts.TemperatureModel(model.type, coefficients)
    fails_without_mutation(ws, [0.21, 0.79], [O2, N2])

print("Mean-air coefficient averaging, input validation, and rollback passed")
