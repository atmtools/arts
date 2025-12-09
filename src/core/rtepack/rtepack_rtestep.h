#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {

/** A single linear step of the point-to-point radiative transfer equation
 * 
 * Sets I := T * (I - J) + J
 * where
 * T = exp(- (K0 + K1) * r / 2)
 * J = (S0 + S1) / 2
 * S0 = planck(T0, f) + inv(K0) * J0
 * S1 = planck(T1, f) + inv(K1) * J1
 *
 * @param[inout] I incoming and then outgoing radiation
 * @param[in] f frequency grid
 * @param[in] K0 absorption matrix at point 0
 * @param[in] K1 absorption matrix at point 1
 * @param[in] J0 source function at point 0
 * @param[in] J1 source function at point 1
 * @param[in] T0 temperature at point 0
 * @param[in] T1 temperature at point 1
 * @param[in] r distance between the two points
 */
void nlte_step(stokvec_vector_view I,
               const Vector &f,
               const propmat_vector_const_view &K0,
               const propmat_vector_const_view &K1,
               const stokvec_vector_const_view &J0,
               const stokvec_vector_const_view &J1,
               const Numeric &T0,
               const Numeric &T1,
               const Numeric &r);

void two_level_linear_emission_step_by_step_full(
    stokvec_vector &I,
    std::vector<stokvec_matrix> &dI,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<muelmat_vector> &Pi,
    const std::vector<muelmat_tensor3> &dTs,
    const std::vector<stokvec_vector> &Js,
    const std::vector<stokvec_matrix> &dJs,
    const stokvec_vector &I0);

void two_level_linear_emission_step_by_step_full(
    std::vector<stokvec_vector> &Is,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<stokvec_vector> &Js);

void two_level_linear_transmission_step(stokvec_vector &I,
                                        std::vector<stokvec_matrix> &dI,
                                        const std::vector<muelmat_vector> &Ts,
                                        const std::vector<muelmat_vector> &Pi,
                                        const std::vector<muelmat_tensor3> &dTs,
                                        const stokvec_vector &I0);
}  // namespace rtepack
