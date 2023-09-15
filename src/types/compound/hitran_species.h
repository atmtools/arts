#ifndef hitran_species_h
#define hitran_species_h

#include "quantum_numbers.h"
#include "enums.h"

namespace Hitran {
  /** Finds the ID of the ARTS species from HITRAN */
  QuantumIdentifier id_from_lookup(Index mol, char isochar);
  
  /** Finds the isotopologue ratio of the species from HITRAN
   * 
   * Requires that id_from_lookup is called first for any error handling
   * to happen...
   */
  Numeric ratio_from_lookup(Index mol, char isochar);
  
  SpeciesIsotopologueRatios isotopologue_ratios();
} // namespace Hitran

#endif
