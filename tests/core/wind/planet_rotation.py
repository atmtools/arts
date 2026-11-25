"""
Tests the rotation of the planet, rot_u of atm_field, as set by atm_fieldRigidPlanetRotation.

The rotation is an AtmRigidPlanetRotation, so the field can be saved to XML and
read back, and the same rotation can be set from Python.

The gas moves with the wind and the rotation, so freq_gridWindShift shifts the
frequencies by their sum, but the frequency wind shift Jacobian is that of the
wind.

For a rigid rotation with angular speed Omega, the velocity of the gas at r
along the direction n is Omega z . (r x n), away from the observer.  So the
frequencies of a path point at r are scaled by 1 + Omega z . (r x n) / c, at
every latitude.  This is tested from 15 Jupiter radii, on rays across the disk
and past the limbs, the polar limb included.  Along a straight ray from an
observer at (d, 0, 0), z . (r x n) = d n_y is the same at every point.  The
points of a geometric path are a few metres off that line, so the test takes
r from the positions of the path, and checks the constancy to that accuracy.

Profile2Path takes the wind and rotation of each path point from atm_field, as
atm_pathFromPath does, so the two give the same path frequencies.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

C = pyarts.arts.constants.c
PERIOD = 9.9249125 * 3600  # Jupiter [s]
OMEGA = 2 * np.pi / PERIOD

# %% The rotation is added to the wind, but the Jacobian is of the wind

ws = pyarts.Workspace()
ws.ray_point.los = [63.0, 137.0]


def shift(wind, rot):
    ws.atm_point.wind = wind
    ws.atm_point.rot_u = rot
    ws.freq_grid = [1e12]
    ws.freq_gridWindShift()
    return ws.freq_grid[0], np.array(ws.freq_wind_shift_jac)


for wind, rot, total in [([0.0, 0.0, 0.0], 12.5e3, [12.5e3, 0.0, 0.0]), ([100.0, -40.0, 3.0], 12.5e3, [12.6e3, -40.0, 3.0])]:
    f_rot, jac_rot = shift(wind, rot)
    f_total, jac_total = shift(total, 0.0)
    assert f_rot == f_total, (f_rot, f_total)
    assert np.array_equal(jac_rot, jac_total), (jac_rot, jac_total)

# %% Rigid rotation along straight rays from far away

ws = pyarts.Workspace()
ws.surf_fieldJupiter(model="Ellipsoid", surf_elevation=0.0)
a, b = ws.surf_field.ellipsoid
E2 = 1 - (b / a) ** 2
TOA = 1000e3
ws.atm_field = pyarts.arts.AtmField(toa=TOA)
ws.atm_fieldRigidPlanetRotation(planet_rotation_period=PERIOD)

# %% The rotation is saved with the field, and is the same as from Python

ALT, LAT, LON = np.meshgrid([0.0, 500e3, TOA], [-89.0, -45.0, 0.0, 30.0, 90.0], [-120.0, 0.0, 77.0], indexing="ij")
rot = ws.atm_field["rot_u"].data(ALT, LAT, LON)

with tempfile.TemporaryDirectory() as tmp:
    filename = os.path.join(tmp, "atm_field.xml")
    ws.atm_field.savexml(filename)
    read = pyarts.arts.AtmField.fromxml(filename)
assert np.array_equal(read["rot_u"].data(ALT, LAT, LON), rot)

python = pyarts.arts.AtmField(toa=TOA)
python["rot_u"] = pyarts.arts.NumericTernaryOperator(pyarts.arts.AtmRigidPlanetRotation(period=PERIOD, ell=[a, b]))
assert np.array_equal(python["rot_u"].data(ALT, LAT, LON), rot)
assert np.allclose(rot[:, 2, :], OMEGA * (a + ALT[:, 2, :]), rtol=1e-15, atol=0), "The equator rotates at Omega (a + h)"
ws.max_stepsize = 100e3
ws.ray_path_observer_agendaSetGeometric(max_step_option="step")
ws.freq_grid = [1e9, 2e9]

d = 15 * a
POS = [d - a, 0.0, 0.0]

# Points of the sky plane through the centre of Jupiter [m], east (y) and north (z): the disk, and the equatorial
# and polar limbs 300 km above the surface
SKY = [(-0.5 * a, 0.0), (-0.5 * a, 0.8 * b), (0.3 * a, -0.9 * b), (-0.2 * a, 0.95 * b), (-(a + 300e3), 0.0), (0.0, b + 300e3)]


def direction(y, z):
    n = np.array([-d, y, z])
    return n / np.linalg.norm(n)


def los_of(n):
    """The los at POS, where up, east and north are x, y and z"""
    return [np.degrees(np.arccos(n[0])), np.degrees(np.arctan2(n[1], n[2])) % 360]


def ecef(alt, lat, lon):
    lat, lon = np.radians(lat), np.radians(lon)
    N = a / np.sqrt(1 - E2 * np.sin(lat) ** 2)
    return np.array(
        [(N + alt) * np.cos(lat) * np.cos(lon), (N + alt) * np.cos(lat) * np.sin(lon), (N * (1 - E2) + alt) * np.sin(lat)]
    )


paths = []
lats = []
for y, z in SKY:
    n = direction(y, z)
    ws.ray_path_observer_agendaExecute(obs_pos=POS, obs_los=los_of(n))
    ws.atm_pathFromPath()
    ws.freq_grid_pathFromPath()
    paths.append(ws.ray_path)

    away = np.array([OMEGA * np.cross(ecef(*p.pos), n)[2] for p in ws.ray_path])
    scale = np.array([f[0] / ws.freq_grid[0] for f in ws.freq_grid_path])
    v_err = (scale - (1 + away / C)) * C
    line_err = away - OMEGA * d * n[1]
    lats.extend(p.pos[1] for p in ws.ray_path)
    print(f"Ray through ({y / a:+.3f}, {z / a:+.3f}) RJ: {len(scale)} points, "
          f"rotation {-OMEGA * d * n[1]:9.2f} m/s towards the observer, largest error {np.abs(v_err).max():.1e} m/s, "
          f"of the path off the line {np.abs(line_err).max():.1e} m/s")
    assert len(scale) > 10
    assert np.abs(v_err).max() < 1e-6, v_err
    assert np.abs(line_err).max() < 0.1, line_err

assert min(lats) < -60 and max(lats) > 85, "The rays must cover high latitudes, where the rotation is slow"

# %% Profile2Path, against atm_pathFromPath and freq_grid_pathFromPath


@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def spectral_propmat_agenda(ws):
    ws.spectral_propmatInit()


ws.spectral_propmat_agenda = spectral_propmat_agenda
ws.alt_grid = np.unique([0.0, TOA] + [p.pos[0] for path in paths for p in path])
ws.lat = 0.0
ws.lon = 0.0
ws.ProfileFromAltitude()

beta = OMEGA * (a + TOA) / C
assert ws.freq_grid_profile[0] <= ws.freq_grid[0] * (1 - beta), "The profile must cover the rotation"
assert ws.freq_grid_profile[-1] >= ws.freq_grid[-1] * (1 + beta), "The profile must cover the rotation"

for path in paths:
    ws.ray_path = path
    ws.atm_pathFromPath()
    ws.freq_grid_pathFromPath()
    ref = [np.array(f) for f in ws.freq_grid_path]
    rot_ref = [p.rot_u for p in ws.atm_path]

    ws.Profile2Path()
    assert all(np.array_equal(f, r) for f, r in zip(ws.freq_grid_path, ref))
    assert [p.rot_u for p in ws.atm_path] == rot_ref
