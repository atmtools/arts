Line-by-line implementation
===========================

The physics of line mixing is described in :doc:`concept.absorption.lbl`.
Public interfaces and examples are in :doc:`user.lbl`.

ECS preparation and equivalent lines
------------------------------------

The Hartmann and Makarov kernels receive rotational states and energies already
in matrix order. Prepare these once per band and reuse them for all collision
partners; neither kernel needs a catalogue-order permutation. The common
``e0`` vector contains original lower-state energies, including when Hartmann
interchanges its angular labels. Separate reference-rotor arrays remain indexed
by angular momentum and supply the collision-basis energies.

For Makarov, ``rotational_energy(N)`` supplies the reference rotor and
``level_energy(N, J)`` supplies the resolved spin-triplet energy. Both use the
same molecular constants and ground-state reference. Collision energies do not
replace catalogue energies in optical populations or line frequencies.

``ComputeData`` stores the relaxation matrix as the spectral operator itself,
in the orientation the resolvent and the eigendecomposition use, so nothing is
transposed between the kernels and ``core_calc_eqv``. The element ``W[i, j]``
is the coupling into line ``i`` from line ``j``, and the optical sum rule is a
row closure, ``sum_j dipr[j] * W[i, j] == 0``. Papers differ in which index
they print first, so coupling coefficients taken from the literature may need
transposing before they match this storage; a new kernel should be checked
against the row closure above rather than against a printed equation.
``core_calc_eqv`` solves the eigenvector system for the amplitude coefficients
instead of explicitly forming an inverse. It checks solver status, eigenvector
conditioning, and damping eigenvalues before using the equivalent lines. The sum-rule residual remains a separate diagnostic;
a completed rescaling does not establish closure or positive absorption.

``calculate`` combines the partner matrices before diagonalisation.
``equivalent_values`` keeps separate partner matrices when preparing equivalent
lines at several temperatures for adaptation. These operations serve different
purposes; adding pure-partner spectra does not reproduce matrix mixing.

Batched derivatives
--------------------

Hartmann and Makarov propagate the requested targets together in a derivative
tensor with axes ``[target, row, column]``. Angular couplings and the
unperturbed eigendecomposition are shared. Changes to the centered operator
require resolved, distinct modes and a well-conditioned eigenvector basis;
unresolved degeneracies raise an error. Targets that leave that operator fixed
reuse its eigenbasis. Frequency derivatives follow the other LBL models' wind
convention.

The final profile reuses the ordinary Voigt model's Faddeeva derivative routine,
including its numerical approximation. Do not substitute the
cancellation-prone closed-form derivative. Temperature finite-difference
checks should stay within a partition-function interpolation interval. Weak
mixing perturbations need steps large enough to resolve the equivalent-line
shifts against the carrier frequency in floating-point arithmetic.

NH3 prepared core
------------------

``src/core/lbl/lbl_lineshape_voigt_ecs_hadded.h`` defines the C++ namespace
``lbl::voigt::ecs::hadded``. Its four-term angular kernel and corrections are
specified in :ref:`lbl-ecs-nh3`. This core is independent of catalogue
regrouping and workspace line-shape dispatch.

``rotational_line`` contains upper and lower ``rotational_state`` objects,
each with integer ``J``, nonnegative ``K``, and an ``inversion`` symmetry.
Select physical states allowed by nuclear-spin symmetry before calling the
kernel. The core validates parallel-band transitions and separates ortho and
para blocks.

With ``n`` lines, ``c`` collision channels, and ``q`` derivative targets,
prepare the following arrays:

.. list-table:: Prepared inputs and outputs
   :header-rows: 1
   :widths: 25 20 55

   * - Data
     - Shape
     - Convention
   * - ``W``
     - ``[n, n]``
     - Real relaxation matrix, ``W[to, from]`` — coupling into ``to`` from ``from``; Hz for spectral calculations.
   * - ``e0``, ``Omega_line``
     - ``[n]``
     - Original lower-state energies in J and corresponding adiabatic factors.
   * - ``basis_data.channels``
     - ``[c]``
     - Unique signed ``(L, Mi, Mf)`` channels, with projections divisible by 3.
   * - ``basis_data.Q``, ``basis_data.Omega``
     - ``[c]``
     - Rates in the units of ``W`` and paper-I factors greater than or equal to 1.
   * - ``dW``
     - ``[q, n, n]``
     - All requested matrix derivatives, with the target axis first.
   * - ``derivative_data.dT``
     - ``[q]``
     - Temperature perturbations.
   * - ``de0``, ``dOmega_line``
     - ``[q, n]``
     - Lower-state energy and line adiabatic-factor perturbations.
   * - ``dQ``, ``dOmega_basis``
     - ``[q, c]``
     - Collision-rate and basis adiabatic-factor perturbations.

Empty derivative views mean zero derivatives. ``adiabatic_factors`` also
accepts batched gap and duration derivatives. Use the same energy model for
``e0`` and the supplied gaps; the core does not select nearest-lower levels.
Signed channels are explicit: omitted channels have zero rates, and the core
does not infer symmetries between different dynamical factors.

``relaxation_matrix_offdiagonal`` fills every off-diagonal entry of ``W`` and
``dW`` and preserves their supplied diagonals. It selects the downward direction
from ``e0`` and obtains the reverse element from detailed balance, including
lower-state degeneracies. Ties use a deterministic source choice; an energy
perturbation across a tie may change that branch. No sorting occurs inside the
function.

Initialize the shared Wigner tables before calling the kernel. One
``coupling_kernel`` call owns the complete line-pair loop and its derivative
pages, with a single ``arts_wigner_thread_init`` / ``arts_wigner_thread_free``
lifetime. There is no per-pair scratch initialization.

``sum_rule_diagonal`` is an optional diagnostic estimator. It replaces only the
diagonals of ``W`` and ``dW`` for fixed reduced dipoles; the result depends on the
finite line set. It is not the Hartmann/Makarov off-diagonal rescaling.

The flat ``hadded_*`` Python bindings expose the prepared primal core. The
batched derivative views remain C++ interfaces. ``relaxation_matrix_profile``
constructs the complex frequency/relaxation matrix and delegates to
``ComputeData``; it does not duplicate the eigensolver or Voigt evaluation.

Validation
-----------

``src/core/lbl/test/test_lbl_ecs_hadded.cpp`` covers independently evaluated
angular coefficients, inversion selection, detailed balance, permutations,
and finite differences of the batched derivatives.
``tests/core/lbl/ecs_nh3.py`` exercises a catalogue band through Python, checks
the no-mixing limit against an independent Voigt sum, and produces a plot.
