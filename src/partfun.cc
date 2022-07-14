#include "partfun.h"

namespace PartitionFunctions::detail {

template Numeric partfun_impl<Derivatives::Yes>(
    Numeric T, const Species::IsotopeRecord& ir);

template Numeric partfun_impl<Derivatives::No>(
    Numeric T, const Species::IsotopeRecord& ir);

}  // namespace PartitionFunctions::detail

Numeric Q(Numeric T, const Species::IsotopeRecord& ir) {
  return PartitionFunctions::detail::partfun_impl<
      PartitionFunctions::Derivatives::No>(T, ir);
}

Numeric dQdT(Numeric T, const Species::IsotopeRecord& ir) {
  return PartitionFunctions::detail::partfun_impl<
      PartitionFunctions::Derivatives::Yes>(T, ir);
}
