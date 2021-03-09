#ifndef hitran_species_h
#define hitran_species_h

#include "quantum.h"

namespace Hitran {
  QuantumIdentifier from_lookup(Index mol, char isochar);
}

#endif
