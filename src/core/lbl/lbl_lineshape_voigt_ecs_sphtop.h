#pragma once

#include <matpack.h>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs::sphtop {
/*! Returns the reduced dipole for a spherical top molecule
 *
 * For a spherical top (Td symmetry), the sub-level structure (A1, A2, E,
 * F1, F2) is already resolved at the band level in the ARTS catalog.
 * Within each symmetry sub-band, the lines are characterized only by J
 * (total angular momentum).  The reduced dipole is the l=0 limit of the
 * linear molecule formula:
 *
 *   d(Jf, Ji) = (-1)^{Jf+1} sqrt(2*Jf+1) * 3j(Jf, 1, Ji; 0, 0, 0)
 *
 * This is equivalent to the Hartmann formula with l_i = l_f = 0.
 *
 * @param[in] Jf  Lower state total angular momentum
 * @param[in] Ji  Upper state total angular momentum
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Jf, const Rational Ji);

/*! Compute the off-diagonal elements of the relaxation matrix
 *  for spherical top molecules (CH4, etc.)
 *
 * Uses the ECS-EP (Energy Corrected Sudden with Exponential Power law)
 * formalism.  For spherical tops, the tetrahedral sub-level structure
 * is already separated into distinct ARTS bands (by rovibSym and alpha),
 * so within a single band the coupling is identical to the linear molecule
 * case with l=0.  The 3j symbols simplify to 3j(J, J', L; 0, 0, 0).
 *
 * @param[in,out] W                  Relaxation matrix (diagonal already set)
 * @param[in]     bnd_qid            Band quantum identifier
 * @param[in]     bnd                Band data
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
}  // namespace lbl::voigt::ecs::sphtop
