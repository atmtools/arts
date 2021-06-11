#ifndef jpl_species_h
#define jpl_species_h

#include "quantum.h"

namespace Jpl {
/** Finds the ID of the ARTS species from JPL */
QuantumIdentifier id_from_lookup(Index tag);
}  // Jpl

#endif  // jpl_species_h
