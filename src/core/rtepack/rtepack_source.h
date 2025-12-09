#pragma once

#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack::source {
stokvec level_lte(Numeric B);

stokvec level_lte(stokvec_vector_view dj, Numeric B, const ConstVectorView &dB);

stokvec level_nlte(Numeric B, const propmat &k, const stokvec &n);

stokvec level_nlte(stokvec_vector_view dj,
                   Numeric B,
                   const ConstVectorView &dB,
                   const propmat &k,
                   const propmat_vector_view &dk,
                   const stokvec &n,
                   const stokvec_vector_view &dn);

stokvec unpolarized_basis_vector();

stokvec unpolarized_basis_vector(stokvec_vector_view dj);

void level_nlte_and_scattering_and_sun(stokvec_vector_view J,
                                       stokvec_matrix_view dJ,
                                       const stokvec_vector_const_view &J_add,
                                       const propmat_vector_const_view &K,
                                       const stokvec_vector_const_view &a,
                                       const stokvec_vector_const_view &S,
                                       const propmat_matrix_const_view &dK,
                                       const stokvec_matrix_const_view &da,
                                       const stokvec_matrix_const_view &dS,
                                       const ConstVectorView &B,
                                       const ConstMatrixView &dB);

void level_nlte_and_scattering(stokvec_vector_view J,
                               stokvec_matrix_view dJ,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const propmat_matrix_const_view &dK,
                               const stokvec_matrix_const_view &da,
                               const stokvec_matrix_const_view &dS,
                               const ConstVectorView &B,
                               const ConstMatrixView &dB);

void level_nlte_and_scattering(stokvec_vector_view J,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const ConstVectorView &B);

void level_nlte(stokvec_vector_view J,
                stokvec_matrix_view dJ,
                const propmat_vector_const_view &K,
                const stokvec_vector_const_view &S,
                const propmat_matrix_const_view &dK,
                const stokvec_matrix_const_view &dS,
                const ConstVectorView &f,
                const Numeric &t,
                const Index &it);
}  // namespace rtepack::source
