#include "partfun.h"

namespace PartitionFunctions {
namespace detail {

template Numeric partfun_impl<Derivatives::Yes>(
    Numeric T, const Species::IsotopeRecord& ir);

template Numeric partfun_impl<Derivatives::No>(
    Numeric T, const Species::IsotopeRecord& ir);

}  // namespace detail

Numeric Q(Numeric T, const Species::IsotopeRecord& ir) {
  return detail::partfun_impl<Derivatives::No>(T, ir);
}

Numeric dQdT(Numeric T, const Species::IsotopeRecord& ir) {
  return detail::partfun_impl<Derivatives::Yes>(T, ir);
}
}  // namespace PartitionFunctions
