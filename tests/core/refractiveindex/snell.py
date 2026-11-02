import numpy as np
import pyarts3 as pyarts

# ray_pointPastRefractive should follow n sin(za) = const at each point,
# and then take the same geometric step as ray_pointPastGeometric

ws = pyarts.Workspace()
ws.atm_fieldInit(toa=100e3)
ws.surf_fieldEarth()
ws.max_stepsize = 1.0
ws.single_propmat = pyarts.arts.Propmat()


def step(option, za, n_old, n_new):
    start = pyarts.arts.PropagationPathPoint()
    start.pos_type = "atm"
    start.los_type = "atm"
    start.pos = [50e3, 0.0, 0.0]
    start.los = [za, 0.0]
    start.nreal = n_old
    ws.ray_path = [start]
    ws.single_dispersion = n_new - 1.0
    ws.ray_point_back_propagation_agendaSet(option=option)
    ws.ray_point_back_propagation_agendaExecute()
    return ws.ray_point


def refract(za, n_old, n_new):
    point = step("RefractiveStepwise", za, n_old, n_new)
    assert np.isclose(point.nreal, n_new, rtol=0, atol=1e-15)
    return point.los[0]


def expected(za_refracted):
    return step("GeometricStepwise", za_refracted, 1.0, 1.0).los[0]


def snell(za, n_old, n_new):
    s = n_old / n_new * np.sin(np.radians(za))
    return 180 - np.degrees(np.arcsin(s)) if za > 90 else np.degrees(np.arcsin(s))


# Vertical rays do not bend
for n_old, n_new in [(1.0, 1.001), (1.001, 1.0)]:
    assert np.isclose(refract(180.0, n_old, n_new), 180.0, rtol=0, atol=1e-10)
    assert np.isclose(refract(0.0, n_old, n_new), 0.0, rtol=0, atol=1e-10)

# Oblique rays bend towards the normal into higher n and away from it into lower n
for za in [20.0, 80.0, 100.0, 160.0]:
    for n_old, n_new in [(1.0, 1.001), (1.001, 1.0)]:
        za_snell = snell(za, n_old, n_new)
        za_new = refract(za, n_old, n_new)
        assert np.isclose(za_new, expected(za_snell), rtol=0, atol=1e-10), (za, n_old, n_new, za_new)
        assert (za_snell > 90) == (za > 90)
        if n_new > n_old:
            assert abs(za_snell - 90) > abs(za - 90)
        else:
            assert abs(za_snell - 90) < abs(za - 90)

# Total internal reflection at grazing incidence into lower n
assert np.isclose(refract(89.0, 1.001, 1.0), expected(91.0), rtol=0, atol=1e-10)
assert np.isclose(refract(91.0, 1.001, 1.0), expected(89.0), rtol=0, atol=1e-10)
