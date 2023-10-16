#pragma once

#include <matpack.h>

#include "rtepack_mueller_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
Array<muelmat_vector> bulk_backscatter(const ConstTensor5View &Pe,
                                       const ConstMatrixView &pnd);

Array<muelmat_matrix> bulk_backscatter_derivative(const ConstTensor5View &Pe,
                                                  const ArrayOfMatrix &dpnd_dx);

void bulk_backscatter_commutative_transmission_rte(
    Array<stokvec_vector> &I,
    Array<Array<stokvec_matrix>> &dI,
    const stokvec_vector &I_incoming,
    const Array<muelmat_vector> &T,
    const Array<muelmat_vector> &PiTf,
    const Array<muelmat_vector> &PiTr,
    const Array<muelmat_vector> &Z,
    const Array<muelmat_matrix> &dT1,
    const Array<muelmat_matrix> &dT2,
    const Array<muelmat_matrix> &dZ);
}  // namespace rtepack