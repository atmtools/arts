"""Moving owned C++ operator results must preserve Python-owned values."""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
ws = pyarts.Workspace()
for name, group in (
    ("atm_field", arts.AtmField),
    ("abs_bands", arts.AbsorptionBands),
    ("measurement_sensor", arts.ArrayOfSensorObsel),
    ("surf_field", arts.SurfaceField),
    ("subsurf_field", arts.SubsurfaceField),
    ("jac_targets", arts.JacobianTargets),
):
    setattr(ws, name, group())
ws.surf_field.ellipsoid = [2, 1]
ws.model_state_vec = [0.5, -0.25]
ws.inversion_iterate_agenda_counter = 0

fit_values = np.array([1.0, 2.0, 3.0])
jac_values = np.array([[1.0, 2.0], [2.0, -1.0], [1.0, 1.0]])
retained_fit = arts.Vector(fit_values)
retained_jac = arts.Matrix(jac_values)
retained_inputs = []


def forward(atm, bands, sensor, surf, subsurf, targets, state, do_jac, counter):
    np.testing.assert_array_equal(state, [0.5, -0.25])
    np.testing.assert_array_equal(surf.ellipsoid, [2, 1])
    retained_inputs.append((state, surf))
    jac = retained_jac if int(do_jac) else arts.Matrix()
    # The same Python objects can be returned repeatedly and kept by the caller.
    return atm, bands, sensor, surf, subsurf, retained_fit, jac


operator = arts.inversion_iterate_agendaOperator(forward)
for do_jac in (1, 1, 0):
    ws.inversion_iterate_agendaExecuteOperator(
        do_jac=do_jac, inversion_iterate_agenda_operator=operator
    )
    np.testing.assert_array_equal(ws.measurement_vec_fit, fit_values)
    if do_jac:
        np.testing.assert_array_equal(ws.measurement_jac, jac_values)
        ws.measurement_jac[0, 0] = 99
    else:
        assert np.asarray(ws.measurement_jac).shape == (0, 0)

    # Output storage must be independent of retained Python return values.
    ws.measurement_vec_fit[0] = 99
    np.testing.assert_array_equal(retained_fit, fit_values)
    np.testing.assert_array_equal(retained_jac, jac_values)
    np.testing.assert_array_equal(ws.model_state_vec, [0.5, -0.25])

# Passing an inout field through the tuple must not consume Python references
# to the original argument, including when the output is subsequently changed.
ws.surf_field.ellipsoid = [4, 3]
for state, surf in retained_inputs:
    np.testing.assert_array_equal(state, [0.5, -0.25])
    np.testing.assert_array_equal(surf.ellipsoid, [2, 1])
