"""Polarized emission of a flat Fresnel surface on a perfectly spherical Moon.

There is no absorption, so the observed radiance is the surface emission plus
the reflected cosmic background.  On a sphere of radius R, a ray leaving an
observer at altitude h with nadir angle a reaches the surface at the incidence
angle inc given by sin(inc) = (R + h) / R * sin(a), so the vertical and
horizontal brightness temperatures are computed at exact incidence angles for
several refractive indices.
"""

import os

import numpy as np
import pyarts3 as pyarts

T_SURF = 300.0
FREQ = 89e9
ALT = 100e3
N2S = [1.5, 3.0, 9.0]
INCIDENCE = [0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 89]

ws = pyarts.workspace.Workspace()
ws.freq_grid = [FREQ]

# %% No absorption

ws.abs_speciesSet(species=[])
ws.spectral_propmat_agendaAuto()

ws.surf_fieldMoon(model="Sphere")
ws.surf_field["t"] = T_SURF

ws.atm_fieldInit(toa=200e3)
ws.atm_field["t"] = 0.0
ws.atm_field["p"] = 0.0

ws.spectral_rad_surface_agendaSet(option="SurfaceReflectance")
ws.spectral_surf_refl_agendaSet(option="FlatRealFresnel")
ws.spectral_rad_transform_operatorSet(option="RJBT")

R = ws.surf_field.ellipsoid[0]

# %% Calculations, Stokes vector per incidence angle


def linpol(stokes, psi):
    """Brightness temperature [K] of linear polarization at psi [deg] from vertical

    The Stokes reference plane is rotated by psi, and the vertical component is
    taken in the rotated basis.
    """
    rotation = pyarts.arts.rtepack.stokes_rotation(np.cos(np.deg2rad(2 * psi)), np.sin(np.deg2rad(2 * psi)))
    return np.array([(rotation @ s) @ pyarts.arts.Stokvec("Iv") for s in stokes])


STOKES = {}
for n2 in N2S:
    ws.surf_field["scalar refractive index"] = n2
    STOKES[n2] = np.empty((len(INCIDENCE), 4))
    for j, inc in enumerate(INCIDENCE):
        nadir = np.rad2deg(np.arcsin(R / (R + ALT) * np.sin(np.deg2rad(inc))))
        ws.ray_pathGeometric(pos=[ALT, 0, 0], los=[180.0 - nadir, 0.0], max_stepsize=1000.0)
        ws.spectral_radClearskyEmission()
        ws.spectral_radApplyUnitFromSpectralRadiance()
        STOKES[n2][j] = ws.spectral_rad[0]

TB = {n2: np.column_stack([linpol(STOKES[n2], 0), linpol(STOKES[n2], 90)]) for n2 in N2S}

# %% Plot

if "ARTS_HEADLESS" not in os.environ:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, len(N2S), sharey=True, figsize=(13, 4.5), constrained_layout=True)
    for a, n2 in zip(ax, N2S):
        for psi in range(10, 90, 10):
            a.plot(INCIDENCE, linpol(STOKES[n2], psi), color="0.6", lw=0.8,
                   label="Linear, every 10° from Tv to Th" if psi == 10 else None)
        a.plot(INCIDENCE, TB[n2][:, 0], "o-", color="C0", lw=2, label="Tv (0°)")
        a.plot(INCIDENCE, TB[n2][:, 1], "o-", color="C1", lw=2, label="Th (90°)")
        a.set_title(f"n = {n2:g}")
        a.set_xlabel("Incidence angle [deg]")
        a.grid(alpha=0.3)
    ax[0].set_ylabel("RJBT [K]")
    ax[0].legend(loc="lower left")
    fig.suptitle(f"Flat Fresnel surface on a spherical Moon: {FREQ / 1e9:g} GHz, $T_s$ = {T_SURF:g} K, no atmosphere")
    plt.show()

# %% Test

expected = {
    1.5: [
        [285.999656, 285.999656],
        [286.482908, 285.507213],
        [287.942875, 283.898433],
        [290.376879, 280.718773],
        [293.623136, 274.973330],
        [296.896818, 264.619757],
        [297.334691, 245.472923],
        [285.260647, 208.966563],
        [227.596394, 138.044745],
        [151.499430, 80.550477],
        [40.029451, 19.086934],
    ],
    3.0: [
        [223.683473, 223.683473],
        [225.197651, 222.164125],
        [229.830249, 217.450434],
        [237.836711, 209.059721],
        [249.560917, 196.137196],
        [265.140481, 177.390058],
        [283.403805, 151.021981],
        [297.425876, 114.711170],
        [273.073149, 65.704643],
        [202.829286, 35.531979],
        [60.287721, 8.360046],
    ],
    9.0: [
        [107.953419, 107.953419],
        [109.251450, 106.666626],
        [113.300180, 102.792306],
        [120.607127, 96.292055],
        [132.185716, 87.113449],
        [149.899422, 75.207965],
        [177.178029, 60.556772],
        [220.197910, 43.204116],
        [283.181632, 23.296738],
        [293.752529, 12.469241],
        [141.013989, 3.432732],
    ],
}

for n2 in N2S:
    assert np.allclose(TB[n2], expected[n2], rtol=0, atol=1e-5), \
        f"Values have drifted from expected results for n = {n2}:\n{TB[n2]}"
