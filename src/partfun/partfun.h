#ifndef partfun_h
#define partfun_h

#include <isotopologues.h>

namespace PartitionFunctions {
Numeric Q(Numeric T, const SpeciesIsotope& ir);

Numeric dQdT(Numeric T, const SpeciesIsotope& ir);

bool has_partfun(const SpeciesIsotope& ir) noexcept;
}  // namespace PartitionFunctions

#endif  // partfun_h
