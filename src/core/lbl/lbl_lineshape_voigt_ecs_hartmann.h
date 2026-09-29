#pragma once

#include <matpack.h>

#include <span>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs {
struct rotational_line;
struct energy_data;
}  // namespace lbl::voigt::ecs

namespace lbl::voigt::ecs::hartmann {
Numeric reduced_dipole(
    const Rational Jf, const Rational Ji, const Rational lf, const Rational li, const Rational k = Rational{1});

//! CO2-626 reference-rotor and resolved-state energies [J].
Numeric rotational_energy(Rational J);
Numeric level_energy(Rational J);

//! Validate the CO2-626 linear-rotor model domain of a whole band.
//! This depends only on the catalogue and is not repeated per atmospheric point.
void validate_band(const QuantumIdentifier& bnd_qid, const band_data& bnd);

//! Prepare lower-state and reference-rotor energies before any angular-kernel swap.
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
}  // namespace lbl::voigt::ecs::hartmann
