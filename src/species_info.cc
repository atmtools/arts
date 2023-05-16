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

Numeric get_lande_spin_constant(const Species::Species species) noexcept {
  if (Species::fromShortName("O2") == species)
    return 2.002064;
  if (Species::fromShortName("NO") == species)
    return 2.00071;
  if (Species::fromShortName("OH") == species)
    return 2.00089;
  if (Species::fromShortName("ClO") == species)
    return 2.00072;
  if (Species::fromShortName("SO") == species)
    return 2.002106;
  return 2.00231930436182;
}

Numeric get_lande_lambda_constant() noexcept { return 1.0; }
