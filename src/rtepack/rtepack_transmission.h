#pragma once

#include "matpack_view.h"

#include "rtepack_multitype.h"
#include "rtepack_source.h"

#include <concepts>
#include <type_traits>
#include <vector>

namespace rtepack {
void two_level_exp(mueller &t,
                   mueller_view &dt1,
                   mueller_view &dt2,
                   const propmat &k1,
                   const propmat &k2,
                   const const_propmat_view &dk1,
                   const const_propmat_view &dk2,
                   const Numeric r,
                   const ExhaustiveConstVectorView &dr1,
                   const ExhaustiveConstVectorView &dr2);
} // namespace rtepack
