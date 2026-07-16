#include "partfun.h"

#include <auto_partfun.h>
#include <auto_partition_functions.h>

namespace PartitionFunctions {
Numeric Q(Numeric T, const SpeciesIsotope& ir) { return partfun_impl<Derivatives::No>(T, ir); }

Numeric dQdT(Numeric T, const SpeciesIsotope& ir) { return partfun_impl<Derivatives::Yes>(T, ir); }
}  // namespace PartitionFunctions
