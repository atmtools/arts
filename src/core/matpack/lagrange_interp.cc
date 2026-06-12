#include "lagrange_interp.h"

namespace lagrange_interp {
static_assert(grid_transformer<grid_identity>);
static_assert(value_transformer<value_identity>);
static_assert(value_transformer<log_transform>);
static_assert(cyclic<loncross>);
static_assert(!cyclic<grid_identity>);
static_assert(constant_lagrange_type_list<std::array<lag_t<-1>, 5>>);

Numeric log_transform::forward(Numeric x) noexcept { return std::log(x); }

Numeric log_transform::inverse(Numeric x) noexcept { return std::exp(x); }
}  // namespace lagrange_interp
