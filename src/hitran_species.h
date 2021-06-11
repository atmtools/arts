#ifndef hitran_species_h
#define hitran_species_h

#include "quantum.h"
#include "enums.h"

namespace Hitran {
  /** Hitran type for selection of data map
   * 
   * Note that this should be appended as few time as possible.
   * 
   * Each type here has to correspond to a data-map in the .cc-file
   * 
   * Only add a new type here if HITRAN changes significantly the
   * interpretation of their data
   */
  ENUMCLASS(Type, char,
            Pre2012CO2Change,
            Newest)
  
  /** Finds the ID of the ARTS species from HITRAN */
  QuantumIdentifier id_from_lookup(Index mol, char isochar, Type type);
  
  /** Finds the isotopologue ratio of the species from HITRAN
   * 
   * Requires that id_from_lookup is called first for any error handling
   * to happen...
   */
  Numeric ratio_from_lookup(Index mol, char isochar, Type type);
  
  SpeciesIsotopologueRatios isotopologue_ratios(Type);
}

#endif
