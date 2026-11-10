"""The kernels of measurement_vecFromSensor against a direct sum of the weights

A callback stands in for radiative transfer: the spectral radiance and its
Jacobian are analytical functions of the geometry, the frequency, and the
Stokes component.  Each measurement element is then the dense sum of its
weights times these functions, computed here with numpy, which no kernel's
sparse summation takes part in.

One sensor has elements with weights at several geometries (an antenna, and
one element given weights at a geometry of another placement as well), the
other elements with weights at one geometry only (a pencil beam), so both ways
the kernels add up the elements are covered.  The sensors have different
frequency grids, and the weights are polarized.
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
sensors = arts.sensors

# Rounding of summing the products of weights and radiance in another order
ROUNDING = 16 * np.finfo(float).eps


def radiance(poslos, freq):
    """Spectral radiance [geometry, frequency, Stokes], different in each"""
    alt, lat, lon, zen, azi = (poslos[:, i, None, None] for i in range(5))
    f = freq[None, :, None] / 1e9
    s = np.arange(4)[None, None, :]
    return (
        (s + 1) * (1 + alt / 1e6)
        + (s + 2) * 0.1 * lat
        - 0.05 * lon
        + (s - 1.5) * np.cos(np.radians(zen))
        + 0.01 * (s + 0.5) * np.sin(np.radians(azi))
        + (0.3 + 0.1 * s) * f
    )


def radiance_jac(poslos, freq, nx):
    """Its derivatives [target, geometry, frequency, Stokes], different for each target"""
    i = radiance(poslos, freq)
    f = freq[None, :, None] / 1e9
    return np.stack([(j + 1) * i - j * f for j in range(nx)]) if nx else np.zeros((0, *i.shape))


def forward(local):
    poslos = np.concatenate([np.array(local.get("obs_pos")), np.array(local.get("obs_los"))])[None]
    freq = np.array(local.get("freq_grid"))
    nx = local.get("jac_targets").x_size()
    local.get("spectral_rad").value = arts.StokvecVector(radiance(poslos, freq)[0])
    local.get("spectral_rad_jac").value = arts.StokvecMatrix(radiance_jac(poslos, freq, nx)[:, 0])
    local.get("ray_path").clear()
    local.get("ray_path").append(arts.PropagationPathPoint())


def weights(sensor):
    """The dense weights [element, geometry, frequency, Stokes]"""
    return np.array([np.array(e.dense) for e in sensor.elements])


def expected(sensor_list, nx):
    """The measurement and its Jacobian summed directly, and the rounding allowed for each"""
    y, y_tol, jac, jac_tol = [], [], [], []
    for sensor in sensor_list:
        poslos = np.array(sensor.geometry)
        freq = np.array(sensor.freq_grid)
        w = weights(sensor)
        i = radiance(poslos, freq)
        di = radiance_jac(poslos, freq, nx)
        y.append(np.einsum("egfs,gfs->e", w, i))
        y_tol.append(ROUNDING * np.einsum("egfs,gfs->e", abs(w), abs(i)))
        jac.append(np.einsum("egfs,jgfs->ej", w, di))
        jac_tol.append(ROUNDING * np.einsum("egfs,jgfs->ej", abs(w), abs(di)))
    return np.concatenate(y), np.concatenate(y_tol), np.concatenate(jac), np.concatenate(jac_tol)


def assert_within(name, got, expected, tol):
    err = abs(np.asarray(got) - expected)
    print(f"{name}: largest error {err.max(initial=0.0):.3e}, largest error over tolerance "
          f"{(err / np.where(tol > 0, tol, 1)).max(initial=0.0):.3f}")
    assert np.all(err <= tol), f"{name}: errors {err} exceed {tol}"


ell = arts.planets.Earth.ellipsoid

# Each element of an antenna has weights at the lines of sight its placement samples
vertical = sensors.BoxChannel(arts.AscendingGrid([89e9, 89.5e9, 90e9]))
vertical.polarization.state = [1.0, 1.0, 0.0, 0.0]
antenna = arts.SensorBuilder(
    [vertical, sensors.DiracChannel(183e9)],
    sensors.GaussianAntenna([0.0, 0.5, 1.0], 0.4, azi_grid_size=4),
)(
    [arts.Vector3([800e3, 10.0 + i, 20.0]) for i in range(3)],
    [arts.Vector2([150.0 + 5 * i, 30.0]) for i in range(3)],
    ell,
    name="antenna",
)

# The first element also has a weight at the last geometry, of the last placement
last = len(antenna.geometry) - 1
f183 = int(np.argmin(abs(np.array(antenna.freq_grid) - 183e9)))
without_extra_weight = weights(antenna)[0]
elements = antenna.elements
w = arts.SparseStokvecMatrix(elements[0])
w[last, f183] = arts.Stokvec([0.25, -0.125, 0.0625, 0.03125])
elements[0] = w
antenna.elements = elements
assert np.array(antenna.elements[0].dense)[last, f183, 0] == 0.25

# Each element of a pencil beam has weights at the line of sight of its placement only
pencil = arts.SensorBuilder(
    [sensors.DiracChannel(150e9), sensors.BoxChannel(arts.AscendingGrid([160e9, 161e9]))],
    sensors.PencilBeamAntenna(),
)(
    [arts.Vector3([600e3, -5.0 + 2 * i, 40.0 + i]) for i in range(4)],
    [arts.Vector2([140.0 + 3 * i, 10.0 * i]) for i in range(4)],
    ell,
    name="pencil",
)

# The test covers both: elements with weights at several geometries, and at one only
geometries_per_element = [
    [np.count_nonzero(np.any(e != 0, axis=(1, 2))) for e in weights(s)] for s in (antenna, pencil)
]
assert min(geometries_per_element[0]) > 1, geometries_per_element
assert max(geometries_per_element[1]) == 1, geometries_per_element

ws = pyarts.Workspace()
ws.surf_fieldEarth()
ws.surf_field["t"] = 280.0
ws.atm_fieldInit(toa=0.0)
ws.abs_bands = arts.AbsorptionBands()
ws.measurement_sensors = [antenna, pencil]
ws.spectral_rad_transform_operatorSet(option="1")

agenda = arts.Agenda("spectral_rad_observer_agenda")
agenda.add(
    arts.Method(
        "analytical_radiance",
        arts.CallbackOperator(
            forward,
            ["freq_grid", "jac_targets", "obs_pos", "obs_los"],
            ["spectral_rad", "spectral_rad_jac", "ray_path"],
        ),
    )
)
agenda.finalize(True)
ws.spectral_rad_observer_agenda = agenda

# Three targets, so that the Jacobian has rows to mix up
ws.jac_targetsInit()
ws.jac_targetsAddSurface(target="t")
ws.jac_targetsAddErrorPolyFit(t=np.arange(len(pencil.elements), dtype=float), sensor_elem=1, polyorder=1)
ws.jac_targetsFinalize()
targets = ws.jac_targets
assert targets.x_size() == 3

for jac_targets in (targets, arts.JacobianTargets()):
    ws.jac_targets = jac_targets
    nx = jac_targets.x_size()
    y, y_tol, jac, jac_tol = expected([antenna, pencil], nx)

    # The weight added to the first element is seen by the comparison
    y_without = np.einsum("gfs,gfs->", without_extra_weight,
                          radiance(np.array(antenna.geometry), np.array(antenna.freq_grid)))
    assert abs(y_without - y[0]) > y_tol[0]

    for kernel in ("Low Memory", "High Performance", "Indexed"):
        ws.measurement_vecFromSensor(kernel=kernel)
        assert np.shape(ws.measurement_jac) == (len(y), nx), kernel
        assert_within(f"{kernel}, {nx} targets, measurement_vec", ws.measurement_vec, y, y_tol)
        assert_within(f"{kernel}, {nx} targets, measurement_jac", ws.measurement_jac, jac, jac_tol)

# A weight outside the geometries of its sensor is refused by the indexed kernel
bad = arts.Sensor(pencil)
elements = bad.elements
w = arts.SparseStokvecMatrix(elements[0])
w[len(bad.geometry), 0] = arts.Stokvec([1.0, 0.0, 0.0, 0.0])
elements[0] = w
bad.elements = elements
ws.measurement_sensors = [bad]
ws.jac_targets = arts.JacobianTargets()
try:
    ws.measurement_vecFromSensor(kernel="Indexed")
except RuntimeError:
    pass
else:
    raise AssertionError("The indexed kernel must refuse a weight outside the geometries of its sensor")

print("done!")
