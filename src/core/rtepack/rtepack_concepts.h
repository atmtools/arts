#pragma once

#include <matpack.h>

namespace rtepack {
using vec7   = matpack::cdata_t<Numeric, 7>;
using vec4   = matpack::cdata_t<Numeric, 4>;
using mat44  = matpack::cdata_t<Numeric, 4, 4>;
using cmat44 = matpack::cdata_t<Complex, 4, 4>;
}  // namespace rtepack
