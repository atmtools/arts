#pragma once

#include "matpack_view.h"

#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_source.h"

#include <concepts>
#include <type_traits>
#include <vector>

namespace rtepack {
using muelmat_view = matpack::matpack_view<muelmat, 1, false, false>;
using const_muelmat_view = matpack::matpack_view<muelmat, 1, true, false>;

void two_level_exp(muelmat &t,
                   muelmat_view &dt1,
                   muelmat_view &dt2,
                   const propmat &k1,
                   const propmat &k2,
                   const propmat_vector_const_view &dk1,
                   const propmat_vector_const_view &dk2,
                   const Numeric r,
                   const ExhaustiveConstVectorView &dr1,
                   const ExhaustiveConstVectorView &dr2);

constexpr stokvec linear_step(const muelmat &t, const stokvec &i, const stokvec &j) {
  return t * (i - j) + j;
}

constexpr stokvec
two_level_linear_step(stokvec_vector_view &di1, stokvec_vector_view &di2, const muelmat &t,
            const muelmat &pit, const stokvec &i, const stokvec &j1, const stokvec &j2,
            const const_muelmat_view &dt1, const const_muelmat_view &dt2,
            const stokvec_vector_const_view &dj1, const stokvec_vector_const_view &dj2) {
  ARTS_ASSERT(dt1.size() == dt1.size() and dt1.size() == dj2.size() and
              dt1.size() == dj2.size())

  const stokvec j = avg(j1, j2);

  for (Index k = 0; k < dt1.size(); k++) {
    di1[k] += static_cast<stokvec>(pit * (dt1[k] * (i - j) + dj1[k] - t * dj1[k]));
    di2[k] += static_cast<stokvec>(pit * (dt2[k] * (i - j) + dj2[k] - t * dj2[k]));
  }

  return linear_step(t, i, j);
}
} // namespace rtepack
