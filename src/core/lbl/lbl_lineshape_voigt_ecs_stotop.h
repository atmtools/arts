#pragma once

#include <matpack.h>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs::stotop {
/*! Returns the reduced dipole for a symmetric top molecule
 *
 * For symmetric top pure rotational transitions (ΔK=0, Δv=0),
 * the reduced dipole matrix element is:
 *
 *   d(Jf, Ji, K) = (-1)^{K+Jf} sqrt(2*Jf+1) * 3j(Jf, 1, Ji; K, 0, -K)
 *
 * @param[in] Jf  Lower state total angular momentum
 * @param[in] Ji  Upper state total angular momentum
 * @param[in] K   Projection of angular momentum on molecular axis (same for
 *                upper and lower)
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Jf,
                       const Rational Ji,
                       const Rational K);

/*! Compute the off-diagonal elements of the relaxation matrix
 *  for symmetric top molecules (NH3, PH3, etc.)
 *
 * Uses the ECS-EP (Energy Corrected Sudden with Exponential Power law)
 * formalism.  The coupling between lines ℓ and ℓ' within a K-sub-band
 * involves 3j and 6j symbols with K replacing the vibrational angular
 * momentum quantum number l used in the linear molecule (CO2) case.
 * Lines with different K are not coupled (ΔK=0 collisions).
 * K is read from each line's quantum numbers.
 *
 * @param[in,out] W                  Relaxation matrix (diagonal already set)
 * @param[in]     bnd_qid            Band quantum identifier
 * @param[in]     bnd                Band data (lines carry K quantum number)
 * @param[in]     sorting            Index sorting of the lines
 * @param[in]     broadening_species Which broadening species to use
 * @param[in]     rovib_data         ECS-EP parameters for the species
 * @param[in]     dipr               Reduced dipole moments
 * @param[in]     atm                Atmospheric point
 */
void relaxation_matrix_offdiagonal(MatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesEnum broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm);
}  // namespace lbl::voigt::ecs::stotop
