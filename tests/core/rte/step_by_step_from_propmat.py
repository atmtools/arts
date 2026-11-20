"""spectral_radStepByStepEmissionFromPropmat against the methods it replaces

The method computes the transmission and the source along the path a block of
frequencies at a time.  Its radiance and derivatives must be those of
spectral_tramat_pathFromPath, spectral_rad_srcvec_pathFromPropmat,
spectral_radStepByStepEmission and spectral_rad_jacFromBackground in turn.  The
path data are random: polarized propagation matrices with layers from optically
thin to thick, non-LTE sources, three Jacobian targets and a frequency grid
shifted differently at each point.  The frequencies are not a whole number of
blocks.  Every transmittance option is tested, with and without the
hypsometric derivatives.

A point without a non-LTE source has nlte_on false, and its source is then
left out.  That must give the same to the bit as sources of zero that are used.
"""

import numpy as np
import pyarts3 as pyarts

rng = np.random.default_rng(42)

NP = 37  # Path points
NF = 45  # Frequencies

ws = pyarts.Workspace()
ws.surf_fieldEarth()
ws.atm_field = pyarts.arts.AtmField(toa=100e3)
ws.atm_field["t"] = 250.0
ws.atm_field["p"] = 1e4
ws.atm_field["O2"] = 0.21
ws.atm_field["H2O"] = 1e-3
ws.abs_bands = pyarts.arts.AbsorptionBands()
ws.measurement_sensors = pyarts.arts.ArrayOfSensor()
ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
ws.jac_targetsAddSpeciesVMR(species="O2")
ws.jac_targetsAddSpeciesVMR(species="H2O")
ws.jac_targetsFinalize()
NQ = ws.jac_targets.target_count()
NX = ws.jac_targets.x_size()

# Layers of 100 m with absorption from 1e-5 to 1e-1 per m, so optical depths from 1e-3 to 10
A = 10 ** rng.uniform(-5, -1, (NP, NF))
K = np.zeros((NP, NF, 7))
K[..., 0] = A
K[..., 1:] = A[..., None] * rng.uniform(-0.4, 0.4, (NP, NF, 6))
dK = rng.uniform(-1, 1, (NP, NQ, NF, 7)) * A[:, None, :, None]
S = rng.uniform(-1, 1, (NP, NF, 4)) * A[..., None] * 1e-12
dS = rng.uniform(-1, 1, (NP, NQ, NF, 4)) * A[:, None, :, None] * 1e-12

f = np.linspace(100e9, 101e9, NF)
T = rng.uniform(150, 300, NP)

ws.spectral_propmat_path = [pyarts.arts.SpectralPropmat(K[i], dK[i], S[i], dS[i]) for i in range(NP)]
ws.freq_grid_path = [pyarts.arts.AscendingGrid(f * (1 + rng.uniform(-1e-5, 1e-5))) for _ in range(NP)]
ws.ray_path = [pyarts.arts.PropagationPathPoint() for _ in range(NP)]
ws.atm_path = [pyarts.arts.AtmPoint() for _ in range(NP)]
for i in range(NP):
    ws.ray_path[i].pos[0] = 100.0 * (NP - i)
    ws.ray_path[i].pos_type = "atm"
    ws.ray_path[i].los_type = "atm"
    ws.atm_path[i]["t"] = T[i]
ws.ray_path[-1].los_type = "surface"
ws.spectral_rad_bkg = pyarts.arts.StokvecVector(rng.uniform(0, 1e-15, (NF, 4)))
ws.spectral_rad_bkg_jac = pyarts.arts.StokvecMatrix(rng.uniform(-1e-17, 1e-17, (NX, NF, 4)))


def worst(new, ref):
    """The largest difference over the largest reference value"""
    new, ref = np.asarray(new), np.asarray(ref)
    assert new.shape == ref.shape, f"{new.shape} != {ref.shape}"
    return np.abs(new - ref).max() / np.abs(ref).max()


for option in ["constant", "linsrc", "linprop", "magop", "magop_linsrc"]:
    for hse in [0, 1]:
        ws.rte_option = option

        ws.spectral_tramat_pathFromPath(hse_derivative=hse)
        ws.spectral_rad_srcvec_pathFromPropmat()
        ws.spectral_radStepByStepEmission()
        ws.spectral_rad_jacFromBackground()
        ref = [np.array(ws.spectral_rad, copy=True), np.array(ws.spectral_rad_jac, copy=True),
               np.array(ws.spectral_rad_jac_path, copy=True)]

        ws.spectral_radStepByStepEmissionFromPropmat(hse_derivative=hse)
        new = [np.array(ws.spectral_rad), np.array(ws.spectral_rad_jac), np.array(ws.spectral_rad_jac_path)]

        errors = [worst(n, r) for n, r in zip(new, ref)]
        assert max(errors) < 1e-14, f"{option}, hse_derivative {hse}: relative errors {errors}"

# Sources of zero, used with nlte_on true and left out with nlte_on false
zero_S, zero_dS = np.zeros((NF, 4)), np.zeros((NQ, NF, 4))

for option in ["constant", "linsrc", "linprop", "magop", "magop_linsrc"]:
    ws.rte_option = option
    results = {}
    for nlte_on in [True, False]:
        ws.spectral_propmat_path = [pyarts.arts.SpectralPropmat(K[i], dK[i], zero_S, zero_dS) for i in range(NP)]
        for point in ws.spectral_propmat_path:
            point.nlte_on = nlte_on

        ws.spectral_rad_srcvec_pathFromPropmat()
        ws.spectral_radStepByStepEmissionFromPropmat()
        results[nlte_on] = [np.array(ws.spectral_rad_srcvec_path.J, copy=True), np.array(ws.spectral_rad_srcvec_path.dJ, copy=True),
                            np.array(ws.spectral_rad, copy=True), np.array(ws.spectral_rad_jac, copy=True),
                            np.array(ws.spectral_rad_jac_path, copy=True)]

    assert all(np.abs(x).max() > 0 for x in results[True]), f"{option}: a result is all zero"
    for name, used, left_out in zip(["source", "source derivative", "radiance", "Jacobian", "path Jacobian"], results[True], results[False]):
        assert np.array_equal(used, left_out), f"{option}: the {name} differs when the non-LTE source is left out"
