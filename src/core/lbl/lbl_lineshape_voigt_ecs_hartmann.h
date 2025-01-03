#pragma once

#include <matpack.h>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs::hartmann {
Numeric reduced_dipole(const Rational Jf,
                       const Rational Ji,
                       const Rational lf,
                       const Rational li,
                       const Rational k = 1);

void relaxation_matrix_offdiagonal(MatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesEnum broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm);
}  // namespace lbl::voigt::ecs::hartmann
