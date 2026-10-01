#pragma once

#include <matpack.h>

#include <span>
#include <vector>

// Paper I: Hadded et al., J. Chem. Phys. 116, 7544 (2002), doi:10.1063/1.1463442.
namespace lbl::voigt::ecs::hadded {
//! Inversion symmetry of a rovibrational level (paper I, Eq. 5).
enum class inversion : bool { symmetric, antisymmetric };

//! Prepared physical level; the caller selects levels allowed by nuclear-spin symmetry.
struct rotational_state {
  Index     J{}, K{};
  inversion symmetry{inversion::symmetric};
};

//! One rank-1 parallel-band transition: equal K and opposite inversion symmetry.
struct rotational_line {
  rotational_state upper{}, lower{};
};

//! Signed body-fixed projections; NH3 collision channels have Mi,Mf multiples of 3.
struct collision_channel {
  Index L{}, Mi{}, Mf{};
};

//! Prepared dynamical factors Q(L,Mi,Mf) and paper-I Omega(L,Mi) (Eqs. 13,18,19).
//! Channels are unique; omitted channels are zero. Supply signed projections
//! explicitly: no symmetry of the off-diagonal Q entries is assumed.
//! Q and the output W have the same units. For W in Hz, convert cross sections
//! [m^2] to Q in Hz with number_density * mean_relative_speed / (2*pi).
//! Omega is the paper's factor >= 1, the reciprocal of the CO2/O2 convention.
struct basis_data {
  std::vector<collision_channel> channels{};
  Vector                         Q{}, Omega{};
};

//! All input derivatives have the same leading target axis as dW.
//! Empty views mean zero derivatives; e0 is in J, and T in K.
struct derivative_data {
  ConstVectorView dT{};
  ConstMatrixView de0{}, dQ{}, dOmega_basis{}, dOmega_line{};
};

//! Rigid symmetric-top energy [J], with B and C supplied in J (paper I, Eq. 6).
Numeric rotational_energy(Index J, Index K, Numeric B, Numeric C);

//! Signed reduced dipole in the paper's population convention (Eq. 9).
//! The corresponding population includes g*(2*J_lower+1), not 2*J_upper+1.
Numeric reduced_dipole(const rotational_line& line);

//! Paper-I Eq. 19 for supplied gaps [J] and collision duration length/speed [s].
//! The nearest-lower-level gap prescription is an input, not inferred here.
//! Optional dOmega and dgap are [target,gap]; dduration is [target].
void adiabatic_factors(VectorView      Omega,
                       ConstVectorView gap,
                       Numeric         duration,
                       MatrixView      dOmega    = {},
                       ConstMatrixView dgap      = {},
                       ConstVectorView dduration = {});

//! Paper-I IOS angular kernel (Eqs. 10-12), with detailed balance (Eq. 16) and
//! ECS energy corrections (Eqs. 17-18), for one parallel band in supplied order.
//! Setting all Omega factors to 1 gives IOS with detailed balance enforced.
//! W[i,j] is the coupling into line i from line j, matching ComputeData's stored
//! matrix convention; this is the transpose of the paper's index order.
//! e0 and Omega_line belong to each line's ORIGINAL LOWER state. The same energy
//! model must supply e0 and all adiabatic gaps; no sorting or catalogue access occurs.
//! Downward elements come from the angular kernel; reverse elements use detailed
//! balance including the lower-state rotational degeneracy. Ortho/para blocks
//! do not couple. Diagonals of W and dW are preserved; every off-diagonal is set.
//! No CO2/O2 truncated-band rescaling is applied: the papers use independent widths.
//! At an exact energy tie, derivatives follow a deterministic choice of source;
//! an energy perturbation across the tie may switch the downward approximation.
void relaxation_matrix_offdiagonal(MatrixView                       W,
                                   std::span<const rotational_line> lines,
                                   const basis_data&                basis,
                                   ConstVectorView                  e0,
                                   ConstVectorView                  Omega_line,
                                   Numeric                          T,
                                   Tensor3View                      dW          = {},
                                   const derivative_data&           derivatives = {});

//! Diagnostic width estimate from the optical sum rule (paper I, Eq. 14).
//! Replaces only the diagonal, also for every dW target. The reduced dipoles are
//! fixed state data; the estimate depends on the supplied finite set of lines.
//! This is explicit and optional, not applied by relaxation_matrix_offdiagonal.
void sum_rule_diagonal(MatrixView W, ConstVectorView dipr, Tensor3View dW = {});
}  // namespace lbl::voigt::ecs::hadded
