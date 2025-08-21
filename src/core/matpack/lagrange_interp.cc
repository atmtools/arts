#include "lagrange_interp.h"

namespace lagrange_interp {
static_assert(cyclic<loncross>);
static_assert(!cyclic<identity>);
static_assert(constant_lagrange_type_list<std::array<lag_t<-1>, 5>>);
}  // namespace lagrange_interp
