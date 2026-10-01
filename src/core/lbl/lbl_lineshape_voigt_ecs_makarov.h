#pragma once

#include <matpack.h>

#include <span>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs {
struct rotational_line;
struct energy_data;
}  // namespace lbl::voigt::ecs

namespace lbl::voigt::ecs::makarov {
/*! Returns the reduced dipole
 * 
 * @param[in] Ju Main rotational number with spin of the upper level
 * @param[in] Jl Main rotational number with spin of the lower level
 * @param[in] N Main rotational number of both levels
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N);

//! Validate the O2-66 spin-triplet, constant-N microwave model domain.
void validate_band(const QuantumIdentifier& bnd_qid, const band_data& bnd);

//! O2-66 reference-rotor energy [J], relative to the N=1, J=0 ground state.
//! Allows the formal angular momenta needed by the ECS basis, including L-2.
Numeric rotational_energy(Rational N);

//! Approximate resolved level energy [J] on the same reference as rotational_energy.
Numeric level_energy(Rational N, Rational J);

//! Prepare resolved lower-state energies and the spin-free ECS reference ladder.
void prepare_energies(energy_data& energies, const QuantumIdentifier& qid, std::span<const rotational_line> lines);

//! All per-line inputs follow the matrix ordering.
//! Optional dW holds all target derivatives, with diagonal derivatives already supplied.
void relaxation_matrix_offdiagonal(MatrixView&                      W,
                                   const QuantumIdentifier&         bnd_qid,
                                   std::span<const rotational_line> lines,
                                   Numeric                          T0,
                                   const SpeciesEnum                broadening_species,
                                   const linemixing::species_data&  rovib_data,
                                   const Vector&                    dipr,
                                   const energy_data&               energies,
                                   const AtmPoint&                  atm,
                                   Tensor3View                      dW     = {},
                                   ConstVectorView                  dT     = {},
                                   ConstMatrixView                  dQ     = {},
                                   ConstMatrixView                  dOmega = {});
}  // namespace lbl::voigt::ecs::makarov
