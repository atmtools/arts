#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_source.h"

#include <concepts>
#include <type_traits>
#include <vector>

namespace rtepack {
constexpr stokvec linear_step(const muelmat &t, const stokvec &i, const stokvec &j) {
  return t * (i - j) + j;
}

constexpr stokvec
two_level_linear_step(stokvec_vector_view &di1, stokvec_vector_view &di2, const muelmat &t,
            const muelmat &pit, const stokvec &i, const stokvec &j1, const stokvec &j2,
            const muelmat_vector_const_view &dt1, const muelmat_vector_const_view &dt2,
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

ENUMCLASS(RadiativeTransferSolver, char, Emission, Transmission)

void two_level_linear_emission_step(
    stokvec_vector_view I, stokvec_matrix_view dI1, stokvec_matrix_view dI2,
    const stokvec_vector_const_view &J1, const stokvec_vector_const_view &J2,
    const stokvec_matrix_const_view &dJ1, const stokvec_matrix_const_view &dJ2,
    const muelmat_vector_const_view &T, const muelmat_vector_const_view &PiT,
    const muelmat_matrix_const_view &dT1, const muelmat_matrix_const_view &dT2);

void two_level_linear_emission_step(stokvec_vector_view I,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const muelmat_vector_const_view &T);

void two_level_linear_transmission_step(
    stokvec_vector_view I, stokvec_matrix_view dI1, stokvec_matrix_view dI2,
    const muelmat_vector_const_view &T, const muelmat_vector_const_view &PiT,
    const muelmat_matrix_const_view &dT1, const muelmat_matrix_const_view &dT2);
} // namespace rtepack
