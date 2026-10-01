Line-by-line interfaces
=======================

See :doc:`concept.absorption.lbl` for line-shape and line-mixing physics.
Implementation and derivative conventions are in :doc:`dev.lbl`.

CO2 and O2 ECS models
---------------------

Set a band's ``lineshape`` to ``VP_ECS_HARTMANN`` for the supported CO2 band
model or ``VP_ECS_MAKAROV`` for the O2 microwave model. Initialize Wigner tables
with ``ws.WignerInit()``. ECS line mixing currently does not support Zeeman
splitting within a band.

The Hartmann rotational energy model currently covers CO2-626. Its collision
parameter presets also contain CO2-628 and CO2-636, but those entries alone do
not provide the missing rotational energy models. Makarov covers O2-66.

For each collision partner, provide a line-shape model with its diagonal
widths, its ECS parameters in ``abs_ecs_data``, and its atmospheric volume
mixing ratio. ARTS constructs and combines the partner matrices before
calculating the spectrum. A Bath partner takes the atmospheric fraction left
after the explicit partners. Without Bath, the listed partners' VMRs are
normalised to sum to one. If all are absent, the collision contribution is zero.

If the catalogue supplies only Bath/air widths,
``ws.abs_ecs_dataAddMeanAir(vmrs=..., species=...)`` can form a Bath ECS
parameter set by averaging the specified partners' temperature-model
coefficients. This approximation retains the available Bath widths; it does
not infer separate partner widths. Its nonlinear collision matrix generally
differs from a weighted sum of separately constructed partner matrices.

Mean-air weights must be finite, nonnegative, and sum to one within
``1e-4``; accepted weights are normalised. Zero-weight species may lack data.
Each positive-weight species must have ECS data for every processed
isotopologue, with matching temperature-model types and coefficient counts.
The input species must be explicit partners, not Bath. These weights define
the parameter set; the runtime Bath fraction follows the atmospheric VMR rule.

The supported workspace Jacobians include temperature, pressure, composition,
isotopologue ratios, line frequencies, lower-state energies, Einstein
coefficients, and ``G0``/``D0`` model coefficients. Spectroscopic lower-state
energy targets change optical populations, not the rotational collision-energy
model. Perturbations that require derivatives of unresolved eigenmodes or an
ill-conditioned eigenbasis raise an error.

Rosenkranz adaptation
---------------------

``abs_bandsLineMixingAdaptation`` obtains per-partner Rosenkranz coefficients
from equivalent lines and fits their temperature dependence with the ``POLY``
model. ``rosenkranz_fit_order=1`` retains the first-order strength correction;
``rosenkranz_fit_order=2`` also fits the quadratic strength and frequency
corrections. The reference pressure is ``atm_point.pressure``, in pascals.
Evaluation restores the corresponding powers of pressure.

Check the fit over the intended pressure and temperature range. Per-partner
fits cannot generally reproduce the cross terms of a fully mixed second-order
calculation; see :ref:`lbl-ecs-rosenkranz-fitting` for the equations and limits.

NH3 prepared-core interface
----------------------------

The NH3 parallel-band model is available through flat helpers in
``pyarts3.arts.lbl``. These accept prepared collision data and line states;
there is no NH3 workspace line-shape option or calibrated H2/He preset.
The physical domain and angular expression are in :ref:`lbl-ecs-nh3`.

The helpers are available after::

    import numpy as np
    import pyarts3 as pyarts

    arts = pyarts.arts
    lbl = arts.lbl
    c, h, k = arts.constants.c, arts.constants.h, arts.constants.k


.. list-table:: NH3 helpers in ``lbl``
   :header-rows: 1
   :widths: 47 53

   * - Helper
     - Purpose
   * - ``hadded_rotational_line(Ju, Jl, K, lower_antisymmetric=False)``
     - Parallel-band line with equal upper/lower K and opposite inversion symmetry.
   * - ``hadded_collision_channel(L, Mi, Mf)``
     - Signed angular channel; Mi and Mf are multiples of 3.
   * - ``hadded_basis_data(channels, Q, Omega)``
     - Collision rates and adiabatic factors in channel order.
   * - ``hadded_rotational_energy(J, K, B, C)``
     - Rigid symmetric-top energy; B, C, and the result are in joules.
   * - ``hadded_reduced_dipole(line)``
     - Signed angular dipole in the lower-state population convention.
   * - ``hadded_adiabatic_factors(gap, duration)``
     - Factors from energy gaps in joules and collision duration in seconds.
   * - ``hadded_relaxation_matrix_offdiagonal(lines, basis, e0, Omega_line, T, widths)``
     - Relaxation matrix with supplied diagonal widths and calculated couplings.

Call ``ws.WignerInit()`` before evaluating the angular kernel. Arrays describing
lines must all use the same order. ``e0`` and ``Omega_line`` describe original
lower states; temperature is in kelvin. All adiabatic factors use the paper's
convention, greater than or equal to one. Setting them to one gives IOS with
detailed balance.

For spectra, use Hz for ``Q``, ``widths``, and the returned matrix. Convert
collision cross sections in square metres using
``Q = number_density * mean_relative_speed * cross_section / (2*pi)``.
Convert widths from cm^-1 by multiplying by ``100*c``. Supply every signed
collision channel explicitly; omitted channels have zero rates. Energies and
adiabatic gaps must come from one consistent molecular model.

Prepared-matrix spectra
~~~~~~~~~~~~~~~~~~~~~~~~

``lbl.relaxation_matrix_profile(frequency, f0, W, population, dipole, gd_fac)``
evaluates a prepared real relaxation matrix. ``frequency`` and ``f0`` are in Hz;
``W[from, to]`` is in Hz and includes nonnegative diagonal half-widths. Pressure
shifts may be included in ``f0``. ``gd_fac`` is the Gaussian 1/e half-width
divided by frequency. The returned complex shape does not contain absorber
density, isotopic abundance, or the stimulated-emission factor.

For the NH3 lower-state convention, prepare the dimensionless populations and
signed dipoles as follows, where ``gl`` and ``gu`` are the catalogue statistical
weights, ``A`` is the Einstein coefficient, and ``Qpart`` the partition function::

    population = gl * np.exp(-e0 / (k*T)) / Qpart
    dipole = np.sign(dipr) * c * np.sqrt(A * gu / (8*np.pi*f0**3*gl))
    gd_fac = np.sqrt(arts.constants.doppler_broadening_const_squared * T / mass)

Here ``dipr`` comes from ``hadded_reduced_dipole`` and ``mass`` is the ARTS
isotopologue mass in atomic mass units. With the returned ``shape``, absorption
in inverse metres is::

    absorption = (n_abs * frequency * (-np.expm1(-h*frequency/(k*T)))
                  * np.asarray(shape).real / np.sqrt(np.pi))

``n_abs`` is the number density of the absorbing isotopologue in inverse cubic
metres, including abundance and isotopic ratio. Populations, signed dipoles,
and the matrix must use the same normalization. The equivalent-line Doppler
approximation and its limits are described in :ref:`lbl-ecs-lineshape`.

Plotted nu2 example
~~~~~~~~~~~~~~~~~~~

From the repository root, with the current ``pyarts3`` build, its catalogue
data, NumPy, SciPy, and Matplotlib available, run::

    python tests/core/lbl/ecs_nh3.py

For a run without a display window::

    ARTS_HEADLESS=1 python tests/core/lbl/ecs_nh3.py

The script saves ``ecs_nh3.png`` in the current directory. It gathers twelve
catalogue nu2 Q lines with J <= 3 across both inversion subbranches and plots
mixed and independent-line absorption, plus their differences. Conditions are
296 K, an 85% H2 / 15% He bath, 1 ppm NH3-4111, and densities of 1, 10, and
30 amagat; the plot labels show the corresponding ideal-gas pressures.

The example uses the diagonal collision factors in Table I of
:cite:t:`Hadded2004`, including the equality of positive and negative diagonal
projections. It sets all adiabatic factors to one and uses representative
constant widths: 0.03 cm^-1/amagat for He and 2.7 times that for H2. It is a
core-use demonstration, not a calibrated full-band planetary model. Coupling
to lines outside the low-J subset and NH3 self broadening are omitted.

The test checks detailed balance, separation of the ortho and para blocks,
finite nonnegative spectra, and visible mixing. Its no-mixing result must agree
with an independent sum of Voigt profiles.
