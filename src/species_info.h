/* Copyright 2018, Richard Larsson.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

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
