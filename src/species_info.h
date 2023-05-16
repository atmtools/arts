#include "matpack_concepts.h"
#include "species.h"

/**
 * @file species_info.h
 * @author Richard Larsson
 * @date 2018-03-29
 * 
 * @brief Some molecular constants
 */

/** Get the Lande spin constant
 * 
 * H. Christensen, and L. Veseth, On the High-Precision Zeeman Effect in 02 and SO.
 * Journal of Molecular Spectroscopy 72, 438-444, 1978.
 * 
 * L. Veseth, Relativistic Corrections to the Zeeman Effect in Diatomic Molecules.
 * Journal of Molecular Spectroscopy 66, 259-271, 1977.
 * 
 * @param species Index-mapped specie
 * @return Numeric Lande spin constant
 */
Numeric get_lande_spin_constant(const Species::Species species) noexcept;

/** Get the Lande Lambda constant
 * 
 * @return Numeric Lande Lambda constant
 */
Numeric get_lande_lambda_constant() noexcept;
