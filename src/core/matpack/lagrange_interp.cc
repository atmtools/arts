#include "lagrange_interp.h"

namespace lagrange_interp {
static_assert(cyclic<loncross>);
static_assert(!cyclic<identity>);
}  // namespace lagrange_interp
