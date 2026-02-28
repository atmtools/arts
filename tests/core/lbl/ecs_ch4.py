"""Test ECS (Energy Corrected Sudden) line mixing for CH4 spherical top.

This test exercises the VP_ECS_SPHTOP lineshape for spherical top molecules.
CH4 is a tetrahedral molecule (Td symmetry) whose sub-level structure
(A1, A2, E, F1, F2) is already separated into distinct ARTS bands.
Within each symmetry sub-band, lines couple only through J — identical
to the linear molecule (Hartmann) formalism with l=0.

The test uses the CH4 nu3 fundamental band (~3019 cm-1) A1->A2 symmetry
sub-band, which contains closely-spaced Q-branch lines where line mixing
is clearly visible at high pressure.  Conditions emulate the deep Jovian
atmosphere (T=200 K, P=10 bar, 86% H2 / 14% He).

References
----------
Pieroni, D., Nguyen-Van-Thanh, C. Brodbeck, C. Boulet, J.-M. Hartmann,
    T. Gabard, J.-P. Champion, D. Bermejo, J.-L. Domenech, and
    C. Claveau, "Experimental and theoretical study of line mixing in
    methane spectra. I. The N2-broadened nu3 band at room temperature,"
    J. Chem. Phys., 110(16), 7717-7732, 1999.
    (ECS-EP line mixing parameters for CH4)

Pine, A. S. and T. Gabard, "Multispectrum fits for line mixing in the
    nu3 band Q branch of methane," J. Mol. Spectrosc., 217, 105-114,
    2003.
    (Experimental line mixing data for CH4 nu3 Q-branch)

Hartmann, J.-M., C. Boulet, and D. Robert, "Collisional Effects on
    Molecular Spectra," Elsevier, 2008.
    (General ECS formalism reference, spherical top coupling)
"""

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

print("ECS CH4 spherical top tests")

ws = pyarts.Workspace()
ws.abs_speciesSet(species=["CH4-211"])

print("Reading catalog data...")
ws.ReadCatalogData()

# Select the nu3 fundamental A1->A2 sub-band (15 lines including 5 Q-branch)
# The Q-branch lines are closely spaced (~2-5 cm-1 apart) near band center,
# so at 10 bar (HWHM ~4-5 cm-1) they overlap and line mixing is visible.
bandkey = (
    "CH4-211 alpha 24 1 ElecStateLabel X X n 1 1 "
    "rovibSym A1 A2 v1 0 0 v2 0 0 v3 1 0 v4 0 0 "
    "vibSym F2 A1"
)

ws.abs_bands = {bandkey: ws.abs_bands[bandkey]}
print(
    "Selected band with "
    + str(len(ws.abs_bands[bandkey].lines))
    + " lines"
)

# Deep Jovian atmospheric conditions (~10 bar level)
ws.jac_targets = pyarts.arts.JacobianTargets()
ws.atm_pointInit()
ws.atm_point.temperature = 200  # ~Jupiter at 10 bar level
ws.atm_point.pressure = 10e5  # 10 bar
ws.atm_point[pyarts.arts.SpeciesEnum("CH4")] = 2e-3  # ~0.2%
ws.atm_point[pyarts.arts.SpeciesEnum("H2")] = 0.86  # 86% H2
ws.atm_point[pyarts.arts.SpeciesEnum("He")] = 0.14  # 14% He
ws.atm_point.mag = [0, 0, 0]
ws.ray_point

ws.jac_targetsInit()
ws.WignerInit()

# Set up ECS data for CH4 with H2/He broadeners
ws.abs_ecs_dataInit()
ws.abs_ecs_dataAddPieroni1999()
ws.abs_ecs_dataAddMeanAir(vmrs=[0.86, 0.14], species=["H2", "He"])

# Frequency grid spanning the nu3 P, Q, and R branches
ws.freq_grid = np.linspace(2700 * 29979245800, 3300 * 29979245800, 5001)


def calc(ws, lineshape=None):
    if lineshape is not None:
        ws.abs_bands[bandkey].lineshape = lineshape
    ws.spectral_propmatInit()
    ws.spectral_propmatAddLines()
    return 1.0 * ws.spectral_propmat[:, 0]


freq_cm = ws.freq_grid / 29979245800

plt.clf()

# Standard Voigt (LTE)
pm_lte = calc(ws)
plt.semilogy(freq_cm, pm_lte, label="Online (VP_LTE)", lw=3)

# ECS full relaxation matrix
pm_ecs = calc(ws, "VP_ECS_SPHTOP")
plt.semilogy(freq_cm, pm_ecs, label="ECS (VP_ECS_SPHTOP)")

# Verify ECS produces reasonable absorption
assert np.max(pm_ecs) > 0, "ECS absorption must have positive values"
neg_fraction = np.sum(pm_ecs < 0) / len(pm_ecs)
print(
    "ECS: max="
    + str(np.max(pm_ecs))
    + ", min="
    + str(np.min(pm_ecs))
    + ", negative fraction="
    + str(round(neg_fraction, 3))
)

# Verify that ECS differs meaningfully from LTE in the wings
# At 10 bar the lines overlap substantially, so line mixing should matter
mask = pm_lte > np.max(pm_lte) * 1e-4
rdiff = np.abs((pm_ecs[mask] - pm_lte[mask]) / pm_lte[mask])
print(
    "Relative ECS-LTE diff: max="
    + str(round(np.max(rdiff), 4))
    + ", mean="
    + str(round(np.mean(rdiff), 4))
)
assert np.max(rdiff) > 0.05, (
    "ECS should differ from LTE by at least 5% somewhere"
)

# Reset to no line mixing
ws.abs_bands.clear_linemixing()
pm_nolm = calc(ws, "VP_LTE")
plt.semilogy(freq_cm, pm_nolm, label="No line mixing")

# 1st order Rosenkranz adaptation from ECS
ws.abs_bands[bandkey].lineshape = "VP_ECS_SPHTOP"
ws.abs_bandsLineMixingAdaptation(
    temperatures=np.linspace(100, 300, 11),
    band_key=bandkey,
    rosenkranz_fit_order=1,
)
pm_ros1 = calc(ws, "VP_LTE")
plt.semilogy(
    freq_cm,
    pm_ros1,
    "--",
    label="1st Order Rosenkranz",
)

plt.xlabel("Wavenumber [cm$^{-1}$]")
plt.ylabel("Propagation matrix [1/m]")
plt.title(r"CH$_4$ $\nu_3$ ECS line mixing test (deep Jovian atmosphere)")
plt.legend()

print("ECS CH4 test completed successfully")
