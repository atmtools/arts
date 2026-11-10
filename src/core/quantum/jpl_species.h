#pragma once

#include <quantum.h>

#include <span>

#include "jpl_species_info.h"

namespace Jpl {
struct LineDataMod {
  QuantumIdentifier qid;
  Numeric           T0;
  Numeric           QT0;
};

/** The builtin JPL species information, sorted by JPL species tag */
std::span<const JplSpeciesInfo> builtin_data();

/** The JPL species information sorted by JPL species tag
 *
 * Throws if a tag is given twice, or if a partition function or reference temperature is not positive.
 */
ArrayOfJplSpeciesInfo sorted_data(ArrayOfJplSpeciesInfo data);

/** Finds the ARTS species of a JPL species tag in data sorted by tag */
LineDataMod data_lookup(Index mol, std::span<const JplSpeciesInfo> data);
}  // namespace Jpl
