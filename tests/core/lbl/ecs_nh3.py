"""Plot a low-J NH3 nu2 band example using the Hadded C++ core.

Run this file to display the plots; it also saves ecs_nh3.png. ARTS_HEADLESS
suppresses the window, as in the other ECS tests.

The 12 catalogue Q lines (J <= 3) cover both inversion subbranches. Collision
cross sections are the final He/H2 values in Hadded et al. (2004), Table I,
doi:10.1063/1.1630306. This subspace only needs L <= 6. Off-diagonal Q(L,M,M')
with M != M' are omitted as in that paper; signed diagonal channels use Eq. 3.

This is an IOS-with-detailed-balance demonstration (all Omega = 1). Constant
He half-widths of 0.03 cm^-1/amagat use the representative scale in Appendix A
of Hadded et al. (2002), doi:10.1063/1.1463442; H2 widths are 2.7 times larger.
These are illustrative widths, not a calibrated full-band planetary model.
Couplings to lines outside this low-J subset and NH3 self broadening are omitted.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts
from scipy.constants import atomic_mass
from scipy.special import wofz

arts = pyarts.arts
lbl = arts.lbl
c, h, k = arts.constants.c, arts.constants.h, arts.constants.k
T = 296.0
n_stp = 101325.0 / (k * 273.15)  # One amagat, in molecules/m^3.
nh3_vmr = 1e-6  # Trace NH3-4111 abundance relative to the H2/He bath.

ws = pyarts.Workspace()
ws.abs_speciesSet(species=["NH3-4111"])
ws.ReadCatalogData()
ws.WignerInit()

# The catalogue splits nu2 by inversion and symmetry labels. Gather both
# inversion subbranches before constructing the relaxation matrix.
nu2 = {arts.QuantumNumberType(name): value for name, value in
       [("v1", "0 0"), ("v2", "1 0"), ("v3", "0 0"), ("v4", "0 0")]}
selected = []
for qid, band in ws.abs_bands.items():
    if not all(key in qid.state and str(qid.state[key]) == value
               for key, value in nu2.items()):
        continue
    if arts.QuantumNumberType("vibInv") not in qid.state:
        continue
    inversion = str(qid.state["vibInv"])
    if inversion not in ("s a", "a s"):
        continue
    for line in band.lines:
        Ju, Jl = int(line.qn["J"].upper.value), int(line.qn["J"].lower.value)
        Ku, Kl = int(line.qn["K"].upper.value), int(line.qn["K"].lower.value)
        if Ju == Jl and 0 < Kl == Ku <= Jl <= 3:
            selected.append((line, lbl.hadded_rotational_line(Ju, Jl, Kl, inversion == "s a")))
selected.sort(key=lambda item: item[0].f0)
assert len(selected) == 12, "Expected both low-J nu2 Q subbranches"
lines, quantum = zip(*selected)
isotope = arts.SpeciesIsotope("NH3-4111")
f0 = np.array([line.f0 for line in lines])
e0 = np.array([line.e0 for line in lines])
gl = np.array([line.gl for line in lines])
gu = np.array([line.gu for line in lines])
a = np.array([line.a for line in lines])
dipr = np.array([lbl.hadded_reduced_dipole(q) for q in quantum])

# Hadded uses LOWER-state populations. This signed dipole preserves the
# catalogue Einstein-A line strength, including its statistical weights.
population = gl * np.exp(-e0 / (k * T)) / isotope.Q(T)
dipole = np.sign(dipr) * c * np.sqrt(a * gu / (8 * np.pi * f0**3 * gl))
gd_fac = np.sqrt(arts.constants.doppler_broadening_const_squared * T / isotope.mass)
np.testing.assert_allclose(
    population * dipole**2,
    c**2 / (8 * np.pi) * np.array([line.s(T, isotope.Q(T)) for line in lines]),
    rtol=2e-14, atol=0,
)

# L, |M|, Q_He, Q_H2 [angstrom^2]; He values already include the fitted scale.
q_table = [
    (1, 0, 2.83, 20.48), (2, 0, 12.54, 12.51),
    (3, 0, 2.14, 16.80), (3, 3, 3.37, 4.48),
    (4, 0, 1.07, 0.53), (4, 3, 1.23, 0.00),
    (5, 0, 0.61, 0.00), (5, 3, 1.02, 0.20),
    (6, 0, 0.25, 0.00), (6, 3, 0.15, 0.00), (6, 6, 0.52, 0.90),
]
channels, sigma_he, sigma_h2 = [], [], []
for L, M, q_he, q_h2 in q_table:
    for signed_M in ((0,) if M == 0 else (-M, M)):
        channels.append(lbl.hadded_collision_channel(L, signed_M, signed_M))
        sigma_he.append(q_he * 1e-20)
        sigma_h2.append(q_h2 * 1e-20)

# Build each partner at one amagat of that partner, then apply bath fractions.
W_per_amagat = np.zeros((len(lines), len(lines)))
for partner, fraction, cross_sections, half_width in [
    ("H2-11", 0.85, sigma_h2, 0.081),
    ("He-4", 0.15, sigma_he, 0.030),
]:
    mass = arts.SpeciesIsotope(partner).mass
    reduced_mass = atomic_mass * isotope.mass * mass / (isotope.mass + mass)
    mean_speed = np.sqrt(8 * k * T / (np.pi * reduced_mass))
    Q = n_stp * mean_speed * np.array(cross_sections) / (2 * np.pi)
    basis = lbl.hadded_basis_data(channels, Q, np.ones(len(channels)))
    W_per_amagat += fraction * np.asarray(lbl.hadded_relaxation_matrix_offdiagonal(
        quantum, basis, e0, np.ones(len(lines)), T,
        np.full(len(lines), half_width * 100 * c),
    ))

# This population convention must satisfy detailed balance even across the
# two inversion subbranches. W[i, j] is the coupling into line i from line j,
# so it is the population of the source state j that symmetrises the flux.
# Ortho and para lines remain uncoupled.
flux = W_per_amagat * population[None, :]
np.testing.assert_allclose(flux, flux.T, rtol=3e-13, atol=1e-15 * np.max(abs(flux)))
orthopara = np.array([q.K % 3 == 0 for q in quantum])
assert np.all(W_per_amagat[orthopara[:, None] != orthopara[None, :]] == 0)

wavenumber = np.linspace(925, 975, 6001)
frequency = wavenumber * 100 * c
profiles = []
for density in (1.0, 10.0, 30.0):
    W = density * W_per_amagat
    diagonal = np.diag(np.diag(W))
    scale = (density * n_stp * nh3_vmr * frequency
             * (-np.expm1(-h * frequency / (k * T))) / np.sqrt(np.pi))

    def absorption(matrix):
        shape = arts.lbl.relaxation_matrix_profile(
            frequency, f0, matrix, population, dipole, gd_fac,
        )
        return scale * np.asarray(shape).real

    independent = absorption(diagonal)
    mixed = absorption(W)

    # Check the no-mixing limit independently against a sum of Voigt profiles.
    doppler = gd_fac * f0
    z = (frequency[:, None] - f0 + 1j * np.diag(W)) / doppler
    reference = scale * np.sum(population * dipole**2 * wofz(z).real / doppler, axis=1)
    np.testing.assert_allclose(independent, reference, rtol=2e-10, atol=1e-13 * reference.max())
    assert np.all(np.isfinite(mixed)) and np.min(mixed) >= 0
    effect = np.max(abs(mixed - independent)) / independent.max()
    assert effect > 1e-3, "The example must exercise a visible mixing contribution"
    pressure_bar = density * n_stp * k * T / 1e5  # Ideal-gas pressure.
    print(f"NH3 nu2: {pressure_bar:.2f} bar, max mixing change / peak = {effect:.2%}")
    profiles.append((pressure_bar, independent, mixed))

# Solid: mixed; dashed: independent. Columns show the two inversion branches.
fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained", sharex="col")
for col, (limits, title) in enumerate([
    ((928, 937), "Lower inversion a → upper s"),
    ((964, 972), "Lower inversion s → upper a"),
]):
    for pressure_bar, independent, mixed in profiles:
        color = axes[0, col].plot(
            wavenumber, mixed, label=f"{pressure_bar:.1f} bar",
        )[0].get_color()
        axes[0, col].plot(wavenumber, independent, "--", color=color, alpha=0.75)
        axes[1, col].plot(
            wavenumber, 100 * (mixed - independent) / independent.max(), color=color,
        )
    axes[0, col].set_title(title)
    axes[0, col].legend()
    axes[1, col].axhline(0, color="0.5", linewidth=0.7)
    axes[1, col].set_xlabel("Wavenumber [cm⁻¹]")
    for axis in axes[:, col]:
        axis.set_xlim(*limits)
        axis.grid(alpha=0.25)
axes[0, 0].set_ylabel("Absorption [m⁻¹]")
axes[1, 0].set_ylabel("Mixing change / independent peak [%]")
fig.suptitle(
    "NH₃ ν₂ Q branches, J ≤ 3 — 296 K, 85% H₂ / 15% He, 1 ppm NH₃-4111\n"
    "IOS + detailed balance; representative widths; solid: mixed, dashed: independent",
    fontsize=12,
)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()
