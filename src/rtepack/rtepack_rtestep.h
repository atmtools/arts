#pragma once

#include "matpack_view.h"

#include "rtepack_multitype.h"
#include "rtepack_source.h"

#include <concepts>
#include <type_traits>
#include <vector>

namespace rtepack {
using mueller_view = matpack::matpack_view<mueller, 1, false, false>;
using const_mueller_view = matpack::matpack_view<mueller, 1, true, false>;

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

constexpr stokes linear_step(const mueller &t, const stokes &i, const stokes &j) {
  return t * (i - j) + j;
}

constexpr stokes
two_level_linear_step(stokes_view &di1, stokes_view &di2, const mueller &t,
            const mueller &pit, const stokes &i, const stokes &j1, const stokes &j2,
            const const_mueller_view &dt1, const const_mueller_view &dt2,
            const const_stokes_view &dj1, const const_stokes_view &dj2) {
  ARTS_ASSERT(dt1.size() == dt1.size() and dt1.size() == dj2.size() and
              dt1.size() == dj2.size())

  const stokes j = avg(j1, j2);

  for (Index k = 0; k < dt1.size(); k++) {
    di1[k] += static_cast<stokes>(pit * (dt1[k] * (i - j) + dj1[k] - t * dj1[k]));
    di2[k] += static_cast<stokes>(pit * (dt2[k] * (i - j) + dj2[k] - t * dj2[k]));
  }

  return linear_step(t, i, j);
}
} // namespace rtepack
