#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"

namespace rtepack {
void two_level_exp(muelmat &t,
                   muelmat_vector_view &dt1,
                   muelmat_vector_view &dt2,
                   const propmat &k1,
                   const propmat &k2,
                   const propmat_vector_const_view &dk1,
                   const propmat_vector_const_view &dk2,
                   const Numeric r,
                   const ExhaustiveConstVectorView &dr1,
                   const ExhaustiveConstVectorView &dr2);

void two_level_exp(muelmat_vector_view &t,
                   muelmat_matrix_view &dt1,
                   muelmat_matrix_view &dt2,
                   const propmat_vector_const_view &k1,
                   const propmat_vector_const_view &k2,
                   const propmat_matrix_const_view &dk1,
                   const propmat_matrix_const_view &dk2,
                   const Numeric r,
                   const ExhaustiveConstVectorView &dr1,
                   const ExhaustiveConstVectorView &dr2);

muelmat exp (propmat k);

void two_level_exp(muelmat_vector_view &t,
                   const propmat_vector_const_view &k1,
                   const propmat_vector_const_view &k2,
                   const Numeric r);
}  // namespace rtepack
