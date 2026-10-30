#pragma once

#include <matpack.h>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
muelmat rayleigh_scattering(const Vector2 &los_in, const Vector2 &los_out, const Numeric depolarization_factor);
}  // namespace rtepack
