"""Regression tests for the tested ARTS2 particle-size distributions."""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts

rain_mass = A.ScatteringSpeciesProperty("rain", A.ParticulateProperty.MassDensity)
ice_mass = A.ScatteringSpeciesProperty("ice", A.ParticulateProperty.MassDensity)
ice_mean = A.ScatteringSpeciesProperty("ice", A.ParticulateProperty.DVeq)
ice_number = A.ScatteringSpeciesProperty("ice", A.ParticulateProperty.NumberDensity)
ice_intercept = A.ScatteringSpeciesProperty(
    "ice", A.ParticulateProperty.InterceptParameter
)

sizes = np.arange(5e-6, 5e-3 + 2.5e-6, 5e-6)
water_content = 1e-4
mean_size = 1e-4
number_density = 1e6


def point(prop, value, temperature=273.0):
    result = A.AtmPoint()
    result.temperature = temperature
    result[prop] = value
    return result


def mass_integral(values, a):
    return np.trapezoid(values * a * sizes**3, sizes)


def assert_mass_closure(
    model, atm_point, a, b=3.0, tolerance=1e-2, mass_property=ice_mass
):
    values = np.asarray(model.evaluate(atm_point.temperature, atm_point.ssprops, sizes, a, b))
    assert np.all(np.isfinite(values))
    assert np.isclose(
        mass_integral(values, a),
        atm_point[mass_property],
        rtol=tolerance,
    )
    return values


# The complete set exercised by ARTS2 TestPsds.arts.
rain = point(rain_mass, water_content)
assert_mass_closure(
    A.MGDSingleMoment(rain_mass, "Abel12", 0, 400, False),
    rain,
    523.0,
    mass_property=rain_mass,
)
assert_mass_closure(
    A.MGDSingleMoment(rain_mass, "Wang16", 0, 400, False),
    rain,
    523.0,
    mass_property=rain_mass,
)

ice = point(ice_mass, water_content)
delanoe = A.DelanoeEtAl14(ice_mass)
assert_mass_closure(delanoe, ice, 480.0)
ice_with_dm = point(ice_mass, water_content)
ice_with_dm[ice_mean] = mean_size
delanoe_with_dm = A.DelanoeEtAl14(ice_mass, ice_intercept, ice_mean)
assert_mass_closure(delanoe_with_dm, ice_with_dm, 480.0)

field_tropical = A.FieldEtAl07(ice_mass, "TR", 0, 400)
field_midlatitude = A.FieldEtAl07(ice_mass, "ML", 0, 400)
field_tr = assert_mass_closure(field_tropical, ice, 480.0)
field_ml = assert_mass_closure(field_midlatitude, ice, 480.0)
assert not np.allclose(field_tr, field_ml)

mh97 = A.McFarquharHeymsfield97(ice_mass, 0, 400)
assert_mass_closure(mh97, ice, 480.0)

# Mass-constrained MGD: infer N0 from a fixed lambda, then infer lambda from a
# fixed N0, matching the two ARTS2 test configurations.
mgd_n0 = A.MGDMass(ice_mass, np.nan, 1.0, 1e5, 1.0, 0, 400)
mgd_lambda = A.MGDMass(ice_mass, 1e11, 0.0, np.nan, 1.0, 0, 400)
assert_mass_closure(mgd_n0, ice, 480.0)
assert_mass_closure(mgd_lambda, ice, 480.0)

# Two-moment MGD constrained by mass-weighted mean size.
ice[ice_mean] = mean_size
mgd_mean = A.MGDTwoMoment(
    ice_mass,
    ice_mean,
    A.MGDTwoMomentType.MassMeanSize,
    2.0,
    1.0,
    0,
    400,
)
mean_values = assert_mass_closure(mgd_mean, ice, 480.0)
calculated_mean = np.trapezoid(mean_values * 480.0 * sizes**4, sizes) / mass_integral(
    mean_values, 480.0
)
assert np.isclose(calculated_mean, mean_size, rtol=1e-2)

# Two-moment MGD constrained by mass and total number.
ice[ice_number] = number_density
mgd_number = A.MGDTwoMoment(
    ice_mass,
    ice_number,
    A.MGDTwoMomentType.MassNumberDensity,
    2.0,
    1.0,
    0,
    400,
)
number_values = assert_mass_closure(mgd_number, ice, 480.0)
assert np.isclose(np.trapezoid(number_values, sizes), number_density, rtol=1e-2)


def check_derivative(model, atm_point, prop, a, b=3.0):
    result = model.evaluate_with_derivatives(atm_point.temperature, atm_point.ssprops, sizes, a, b)
    derivative = np.asarray(result.derivative(prop))
    value = atm_point[prop]
    step = 2e-6 * max(abs(value), 1e-8)
    plus = A.AtmPoint(atm_point)
    minus = A.AtmPoint(atm_point)
    plus[prop] = value + step
    minus[prop] = value - step
    finite_difference = (
        np.asarray(model.evaluate(plus.temperature, plus.ssprops, sizes, a, b))
        - np.asarray(model.evaluate(minus.temperature, minus.ssprops, sizes, a, b))
    ) / (2 * step)
    scale = max(np.linalg.norm(finite_difference), 1.0)
    assert np.linalg.norm(derivative - finite_difference) / scale < 2e-5


# Derivatives are exposed by atmospheric property for future F06 use.
for model, atm_point, a in [
    (mgd_n0, ice, 480.0),
    (delanoe, point(ice_mass, water_content), 480.0),
    (field_tropical, ice, 480.0),
    (mh97, ice, 480.0),
]:
    check_derivative(model, atm_point, ice_mass, a)

check_derivative(mgd_mean, ice, ice_mean, 480.0)
check_derivative(mgd_number, ice, ice_number, 480.0)
check_derivative(delanoe_with_dm, ice_with_dm, ice_mean, 480.0)

# Non-strict temperature limits return zero; strict limits reject the call.
cold_rain = point(rain_mass, water_content, temperature=250.0)
non_strict = A.MGDSingleMoment(rain_mass, "Abel12", 273, 373, False)
assert np.all(np.asarray(non_strict.evaluate(cold_rain.temperature, cold_rain.ssprops, sizes, 523.0, 3.0)) == 0.0)

strict = A.MGDSingleMoment(rain_mass, "Abel12", 273, 373, True)
try:
    strict.evaluate(cold_rain.temperature, cold_rain.ssprops, sizes, 523.0, 3.0)
except RuntimeError:
    pass
else:
    raise AssertionError("A strict PSD accepted a temperature outside its range")
