"""The radiance leaving a subsurface by VDISORT (spectral_radSubsurfaceVDisort) against closed forms.

The subsurface absorbs by its refractive index n' and loss tangent, has no
scattering, and its interface refracts and reflects by Fresnel's equations.
With n = sqrt(n'^2 (1 + i tan(delta))), the absorption coefficient is
4 pi f Im(n) / c and the interface that of Re(n).  The reflectances rv and rh
of a ray leaving at theta from the facet normal give, for the radiance L_p and
L_s within at the refracted angle and the sky's B_c from the specular direction,

  I = (1 - (rv + rh) / 2) L + (rv + rh) / 2 B_c,
  Q = -(rv - rh) / 2 L + (rv - rh) / 2 B_c

for unpolarized L within.
"""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts
C = 299792458.0
FREQ = np.array([280e9, 490e9])
N_PRIME = np.array([1.8, 1.9])
LOSS_TANGENT = np.array([0.004, 0.006])
DEPOLARIZATION = np.array([0.1, 0.25])
T0 = 250.0
DEPTH = np.linspace(0.0, -0.03, 31)

ws = pyarts.Workspace()
ws.freq_grid = FREQ
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("h")] = 0.0
ws.surf_field[A.SurfaceKey("t")] = T0
ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)
ws.abs_speciesSet(species=[])
ws.abs_bands = {}
ws.spectral_propmat_agendaAuto()
ws.ray_path_observer_agendaSetGeometric()
ws.spectral_rad_observer_agendaSet(option="Emission")
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")

ws.subsurf_field.bottom_depth = -10.0
ws.subsurf_field[A.SubsurfaceKey("t")] = T0
ws.subsurf_field[A.SubsurfaceKey("rho")] = 1500.0


def spectral(values, grids, names):
    return getattr(A, f"SpectralGeodeticField{len(grids) + 1}")(
        data=np.reshape(values, (len(FREQ),) + (1,) * len(grids)),
        grid_names=["frequency", *names],
        grids=[FREQ, *grids],
    )


ws.subsurf_spectral_field["refractive index"] = spectral(N_PRIME, [[0.0], [0.0], [0.0]], ["alt", "lat", "lon"])
ws.subsurf_spectral_field["loss tangent"] = spectral(LOSS_TANGENT, [[0.0], [0.0], [0.0]], ["alt", "lat", "lon"])
ws.surf_spectral_field["depolarization factor"] = spectral(DEPOLARIZATION, [[0.0], [0.0]], ["lat", "lon"])

ws.quadrature_dimension = 16
ws.fourier_mode_dimension = 1
ws.legendre_polynomial_dimension = 1
ws.measurement_sensors = []
ws.jac_targetsInit()
ws.jac_targetsFinalize()

n = np.sqrt(N_PRIME**2 * (1 + 1j * LOSS_TANGENT))
alpha = 4 * np.pi * FREQ * n.imag / C
B0 = np.array([A.physics.planck(f, T0) for f in FREQ])
BC = np.array([A.physics.planck(f, A.constants.cosmic_microwave_background_temperature) for f in FREQ])


def reflectances(theta):
    """|Rv|^2 and |Rh|^2 of a ray leaving the interface at theta [deg] from its normal"""
    c1 = np.cos(np.radians(theta))
    c2 = np.sqrt(1 - (np.sin(np.radians(theta)) / n.real) ** 2)
    rv = (n.real * c1 - c2) / (n.real * c1 + c2)
    rh = (c1 - n.real * c2) / (c1 + n.real * c2)
    return rv**2, rh**2


def radiance(theta, azimuth=0.0, lat=0.0):
    point = A.PropagationPathPoint()
    point.pos = [0.0, lat, 0.0]
    point.los = [theta, azimuth]
    point.pos_type = point.los_type = A.PathPositionType("surface")
    ws.ray_point = point
    ws.spectral_radSubsurfaceVDisort(depth_profile=DEPTH)
    return np.asarray(ws.spectral_rad).copy()


def expected(theta, inside):
    rv, rh = reflectances(theta)
    mean, diff = 0.5 * (rv + rh), 0.5 * (rv - rh)
    return np.column_stack([(1 - mean) * inside + mean * BC, -diff * inside + diff * BC, 0 * BC, 0 * BC])


# %% An isothermal subsurface is Kirchhoff's, for any optical thickness
ws.subsurf_vdisort_settings_agendaSetup()
for theta in [0.0, 25.0, 50.0, 70.0]:
    got = radiance(theta)
    assert np.allclose(got, expected(theta, B0), rtol=1e-10, atol=1e-30), (theta, got, expected(theta, B0))

# %% Without emission within, the radiation from the deepest level is absorbed along the refracted ray
ws.subsurf_vdisort_settings_agendaSetup(layer_emission_setting="None", sky_setting="None")
for theta in [0.0, 40.0, 65.0]:
    mu_inside = np.sqrt(1 - (np.sin(np.radians(theta)) / n.real) ** 2)
    inside = B0 * np.exp(-alpha * abs(DEPTH[-1]) / mu_inside)
    got = radiance(theta)
    assert np.allclose(got, expected(theta, inside), rtol=1e-9, atol=1e-30), (theta, got, expected(theta, inside))

# %% Every step can be left out: nothing within leaves only the reflected sky
ws.subsurf_vdisort_settings_agendaSetup(layer_emission_setting="None", deep_setting="None", interface_setting="None")
got = radiance(40.0)
assert np.allclose(got, expected(40.0, 0 * B0), rtol=1e-10, atol=1e-30)

# %% Depolarization scales the linear polarization, and its derivative is that of the scaling
ws.subsurf_vdisort_settings_agendaSetup()
ws.jac_targetsInit()
ws.jac_targetsAddSurfaceSpectral(target="depolarization factor")
ws.jac_targetsFinalize()
polarized = radiance(55.0)
ws.spectral_radSurfaceDepolarization()
depolarized = np.asarray(ws.spectral_rad).copy()
jac = np.asarray(ws.spectral_rad_jac).copy()
assert np.allclose(depolarized[:, 0], polarized[:, 0], rtol=1e-14)
assert np.allclose(depolarized[:, 1], (1 - DEPOLARIZATION) * polarized[:, 1], rtol=1e-14)
assert np.allclose(jac[[0, 1], [0, 1], 1], -polarized[:, 1], rtol=1e-14)
assert np.all(jac[[0, 1], [1, 0]] == 0.0)
ws.jac_targetsInit()
ws.jac_targetsFinalize()

# %% A facet tilted along the ray's azimuth only changes the angle from its normal
SLOPE = np.tan(np.radians(10.0))
lat = np.array([-0.01, 0.01])
ws.surf_field[A.SurfaceKey("h")] = A.GeodeticField2(
    data=np.repeat((SLOPE * np.radians(lat) * ws.surf_field.ellipsoid[0])[:, None], 2, axis=1),
    grid_names=["lat", "lon"],
    grids=[lat, [-0.01, 0.01]],
)


def enu(los):
    za, aa = np.radians(los)
    return np.array([np.sin(za) * np.sin(aa), np.sin(za) * np.cos(aa), np.cos(za)])


# The normal of a surface point is the line of sight into the surface along it
normal = -enu(ws.surf_field(0.0, 0.0).normal)
# About 10 deg: the slope is along the meridian, whose radius of curvature is not the equatorial radius
assert abs(np.degrees(np.arccos(normal[2])) - 10.0) < 0.1, ws.surf_field(0.0, 0.0).normal



def check_tilted(los):
    """The radiance along los from the tilted facet is the flat one at the angle from its normal, turned by psi

    psi turns the facet's basis (p, s), s across the plane of incidence, into the ray's meridional basis (v, h).
    Returns sin(2 psi).
    """
    k = enu(los)
    flat = expected(np.degrees(np.arccos(k @ normal)), B0)
    got = radiance(*los)
    h = np.cross(k, [0, 0, 1])
    h /= np.linalg.norm(h)
    v = np.cross(h, k)
    s = np.cross(-(2 * (k @ normal) * normal - k), normal)
    s /= np.linalg.norm(s)
    p = np.cross(s, k)
    cos_psi, sin_psi = v @ p, v @ s
    cos2, sin2 = cos_psi**2 - sin_psi**2, 2 * sin_psi * cos_psi
    assert np.allclose(got[:, 0], flat[:, 0], rtol=1e-9), (los, got, flat)
    assert np.allclose(got[:, 1], cos2 * flat[:, 1], rtol=1e-8), (los, got, flat, cos2)
    assert np.allclose(np.abs(got[:, 2]), np.abs(sin2 * flat[:, 1]), rtol=1e-6, atol=1e-30), (los, got, flat, sin2)
    assert np.all(got[:, 3] == 0.0)
    return sin2


# Along the ray's azimuth, the polarization is not turned (the normal is that of the elevation, to rounding)
for azimuth in [0.0, 180.0]:
    assert abs(check_tilted([30.0, azimuth])) < 1e-3

# %% A facet tilted across the ray's azimuth also turns the polarization by psi into the ray's basis
assert abs(check_tilted([30.0, 90.0])) > 0.05, "the tilt must turn the polarization for the test to mean anything"
ws.surf_field[A.SurfaceKey("h")] = 0.0

# %% The sun in the specular direction is reflected as a glint
theta = 30.0
without = radiance(theta)
ws.lat = -theta
ws.lon = 0.0
ws.sunBlackbody()
ws.suns = [ws.sun]
ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")
with_sun = radiance(theta)
sun = np.asarray(ws.sun.spectrum)[:, 0] / np.pi
rv, rh = reflectances(theta)
assert np.allclose(with_sun[:, 0] - without[:, 0], 0.5 * (rv + rh) * (sun - BC), rtol=1e-9)
assert np.allclose(with_sun[:, 1] - without[:, 1], 0.5 * (rv - rh) * (sun - BC), rtol=1e-9)
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")

# %% Scattering by the species of the subsurface, a Henyey-Greenstein scatterer of its own properties
ws.subsurf_spectral_field["scattering extinction"] = spectral(np.array([60.0, 90.0]), [[0.0], [0.0], [0.0]], ["alt", "lat", "lon"])
ws.subsurf_spectral_field["single scattering albedo"] = spectral(np.array([0.9, 0.8]), [[0.0], [0.0], [0.0]], ["alt", "lat", "lon"])
ws.subsurf_spectral_field["asymmetry parameter"] = spectral(np.array([0.6, 0.4]), [[0.0], [0.0], [0.0]], ["alt", "lat", "lon"])
ws.subsurf_scat_species = [
    A.HenyeyGreensteinScatterer("scattering extinction", "single scattering albedo", "asymmetry parameter")
]
ws.legendre_polynomial_dimension = 16

# In equilibrium with a sky of its own temperature the radiance within is Planck's, whatever the scattering, so the
# outgoing radiance is Kirchhoff's
ws.subsurf_vdisort_settings_agendaSetup(scattering_setting="ScatteringSpecies", sky_temperature=T0)
for theta in [0.0, 35.0, 60.0]:
    got = radiance(theta)
    assert np.allclose(got, expected(theta, B0), rtol=1e-9, atol=1e-30), (theta, got, expected(theta, B0))

# Away from it, the radiance within along the refracted ray is VDISORT's on its stream when the ray is along one
ws.subsurf_vdisort_settings_agendaSetup(
    scattering_setting="ScatteringSpecies", layer_emission_setting="None", sky_setting="None"
)
radiance(0.0)
ws.ray_pathFromPointAndDepth(depth_profile=DEPTH)
ws.subsurf_vdisort_settings_agendaExecute()
assert np.all(np.asarray(ws.vdisort_settings.single_scattering_albedo) > 0.3)
ws.vdisort_spectral_rad_fieldCalc(azi_grid=[0.0])
zenith = np.asarray(ws.vdisort_spectral_rad_field.zen_grid)
stream = int(np.argmax(zenith))  # the stream most nearly towards the interface, propagating at 180 - zenith from it
streams = np.asarray(ws.vdisort_spectral_rad_field.data)[:, -1, 0, stream]
for iv in range(len(FREQ)):
    theta = np.degrees(np.arcsin(n.real[iv] * np.sin(np.radians(180.0 - zenith[stream]))))
    got = radiance(theta)[iv]
    rv, rh = reflectances(theta)
    mean, diff = 0.5 * (rv[iv] + rh[iv]), 0.5 * (rv[iv] - rh[iv])
    # got = [[1 - mean, -diff], [-diff, 1 - mean]] (I, Q)_within + (mean, diff) B_c
    inside = np.linalg.solve([[1 - mean, -diff], [-diff, 1 - mean]], got[:2] - np.array([mean, diff]) * BC[iv])
    assert np.allclose(inside, streams[iv, :2], rtol=1e-9), (iv, inside, streams[iv])
