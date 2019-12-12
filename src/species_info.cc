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

/**
 * @file species_info.cc
 * @author Richard Larsson
 * @date 2018-03-29
 * 
 * @brief Some molecular constants
 */

#include "species_info.h"
#include "absorption.h"
#include "wigner_functions.h"

Numeric get_lande_spin_constant(const Index species) noexcept {
  if (species_index_from_species_name("O2") == species)
    return 2.002064;
  else if (species_index_from_species_name("NO") == species)
    return 2.00071;
  else if (species_index_from_species_name("OH") == species)
    return 2.00089;
  else if (species_index_from_species_name("ClO") == species)
    return 2.00072;
  else if (species_index_from_species_name("SO") == species)
    return 2.002106;
  else
    return 2.00231930436182;
}

Numeric get_lande_lambda_constant() noexcept { return 1.0; }
