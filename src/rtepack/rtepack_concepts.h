#pragma once

#include "matpack_constexpr.h"

namespace rtepack {
using vec7 = matpack::matpack_constant_data<Numeric, 7>;
using vec4 = matpack::matpack_constant_data<Numeric, 4>;
using mat44 = matpack::matpack_constant_data<Numeric, 4, 4>;
}  // namespace rtepack
