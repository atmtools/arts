#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_source.h"

namespace rtepack {
/** A single linear step of the point-to-point radiative transfer equation
 * 
 * Returns T * (I - J) + J
 *
 * @param T Transmission matrix of the medium
 * @param I Radiation attenuated by the medium
 * @param J Source function of the medium
 * @return constexpr stokvec Radiation leaving the medium
 */
constexpr stokvec linear_step(const muelmat &T,
                              const stokvec &I,
                              const stokvec &J) {
  return T * (I - J) + J;
}

/** A single linear step of the point-to-point radiative transfer equation
 * 
 * Returns T * (I - J) + J
 *
 * Also returns the derivatives.
 *
 * The assumption is however that T and J are from two points, 1 and 2.
 * The value of T has to be computed as if it consists of theses two points.
 * The derivatives are computed as if it was.
 * 
 * @param dI1 Radiation attenuated by the medium derivative wrt point 1
 * @param dI2 Radiation attenuated by the medium derivative wrt point 2
 * @param T Transmission matrix of the medium
 * @param PiT Cumulated Transmission to the end-of-the-path, to scale derivatives
 * @param I Radiation attenuated by the medium
 * @param J1 Source function of the medium at point 1
 * @param J2 Source function of the medium at point 2
 * @param dT1 Transmission matrix of the medium wrt point 1
 * @param dT2 Transmission matrix of the medium wrt point 2
 * @param dJ1 Source function derivative of the medium wrt point 1
 * @param dJ2 Source function derivative of the medium wrt point 1
 * @return constexpr stokvec as in linear_step with J = avg(J1, J2)
 */
constexpr stokvec two_level_linear_step(stokvec_vector_view &dI1,
                                        stokvec_vector_view &dI2,
                                        const muelmat &T,
                                        const muelmat &PiT,
                                        const stokvec &I,
                                        const stokvec &J1,
                                        const stokvec &J2,
                                        const muelmat_vector_const_view &dT1,
                                        const muelmat_vector_const_view &dT2,
                                        const stokvec_vector_const_view &dJ1,
                                        const stokvec_vector_const_view &dJ2) {
  ARTS_ASSERT(dI1.size() == dI2.size() and dT1.size() == dT2.size() and
              dJ1.size() == dJ2.size() and dI1.size() == dT1.size())

  const auto J = avg(J1, J2);

  for (Index k = 0; k < dT1.size(); k++) {
    dI1[k] += PiT * (dT1[k] * (I - J) + dJ1[k] - T * dJ1[k]);
    dI2[k] += PiT * (dT2[k] * (I - J) + dJ2[k] - T * dJ2[k]);
  }

  return linear_step(T, I, J);
}

void two_level_linear_emission_step(stokvec_vector_view I,
                                    stokvec_matrix_view dI1,
                                    stokvec_matrix_view dI2,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const stokvec_matrix_const_view &dJ1,
                                    const stokvec_matrix_const_view &dJ2,
                                    const muelmat_vector_const_view &T,
                                    const muelmat_vector_const_view &PiT,
                                    const muelmat_matrix_const_view &dT1,
                                    const muelmat_matrix_const_view &dT2);

void two_level_linear_emission_step(stokvec_vector_view I,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const muelmat_vector_const_view &T);

void two_level_linear_transmission_step(stokvec_vector_view I,
                                        stokvec_matrix_view dI1,
                                        stokvec_matrix_view dI2,
                                        const muelmat_vector_const_view &T,
                                        const muelmat_vector_const_view &PiT,
                                        const muelmat_matrix_const_view &dT1,
                                        const muelmat_matrix_const_view &dT2);

void two_level_linear_emission_step(stokvec_vector &I,
                                    std::vector<stokvec_matrix> &dI,
                                    const std::vector<muelmat_vector> &Ts,
                                    const std::vector<muelmat_vector> &Pi,
                                    const std::vector<muelmat_tensor3> &dTs,
                                    const std::vector<stokvec_vector> &Js,
                                    const std::vector<stokvec_matrix> &dJs,
                                    const stokvec_vector &I0);

void two_level_linear_transmission_step(stokvec_vector &I,
                                        std::vector<stokvec_matrix> &dI,
                                        const std::vector<muelmat_vector> &Ts,
                                        const std::vector<muelmat_vector> &Pi,
                                        const std::vector<muelmat_tensor3> &dTs,
                                        const stokvec_vector &I0);
}  // namespace rtepack
