#pragma once

#include <matpack.h>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs::makarov {
/*! Returns the reduced dipole
 * 
 * @param[in] Ju Main rotational number with spin of the upper level
 * @param[in] Jl Main rotational number with spin of the lower level
 * @param[in] N Main rotational number of both levels
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N);

void relaxation_matrix_offdiagonal(MatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesEnum broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm);
}  // namespace lbl::voigt::ecs::makarov
