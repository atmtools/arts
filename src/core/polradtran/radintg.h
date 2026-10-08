#pragma once

#include <matpack.h>

#include "polradtran_workdata.h"

/* The subroutines that 3rdparty/polradtran/radintg3.f and radintg4.f share,
   ported to C++: RT3's RT3_X is RT4's X, line by line.  Each follows its
   Fortran step by step, with matpack in place of Evans' matrix helpers
   (radmat.f): MZERO is "= 0.0", MCOPY "=", MADD and MSCALARMULT "+=" and
   "*=", MINVERT is inv_inplace (LAPACK), and MMULT is mult (BLAS DGEMM, or
   DGEMV for a matrix-vector product), whose alpha and beta absorb the
   MIDENTITY, MSUB and MADD around a product.  Their scratch (X, Y, GAMMA
   and their vectors) is that of a workdata, which must be sized for their
   n streams (workdata::resize).

   A slab's reflection and transmission matrices are [2, n, n], the
   Fortran's (N, N, 2): the column-major n x n matrices of the + (l = 0) and
   - (l = 1) parts, so a row-major matpack matrix holds the transpose of the
   Fortran matrix and the products are computed for the Fortran matrices
   (C = A B is mult(C, B, A), y = A x is mult(y, transpose(A), x)).  Its
   sources and radiances are [2, n], the Fortran's (N, 2).  The counts are
   not passed; they are the extents of the arrays. */
namespace polradtran {
/** COMBINE_LAYERS (RT3_COMBINE_LAYERS): the reflection and transmission
 * matrices and the source vectors of two layers combined into one.  The
 * positive side (down) of the first layer is attached to the negative side
 * (up) of the second layer; thus layer 1 is put on top of layer 2.
 *
 *   reflect1, trans1, reflect2, trans2  [2, n, n]  REFLECT1(N, N, 2), ...
 *   source1, source2                    [2, n]     SOURCE1(N, 2), SOURCE2(N, 2)
 *   out_reflect, out_trans              [2, n, n]  OUT_REFLECT(N, N, 2), OUT_TRANS(N, N, 2), output
 *   out_source                          [2, n]     OUT_SOURCE(N, 2), output
 *
 * The outputs must not share memory with the inputs.
 */
void combine_layers(ConstTensor3View reflect1,
                    ConstTensor3View trans1,
                    ConstMatrixView  source1,
                    ConstTensor3View reflect2,
                    ConstTensor3View trans2,
                    ConstMatrixView  source2,
                    Tensor3View      out_reflect,
                    Tensor3View      out_trans,
                    MatrixView       out_source,
                    workdata&        work);

/** INTERNAL_RADIANCE (RT3_INTERNAL_RADIANCE): the internal radiance at a
 * level.  The reflection and transmission matrices and source vectors are
 * given for the atmosphere above (up) and below (down) the level.  The
 * upwelling and downwelling radiances are computed from the two layer
 * properties and the radiance incident on the top and bottom.
 *
 *   upreflect, uptrans, downreflect, downtrans  [2, n, n]  UPREFLECT(N, N, 2), ...
 *   upsource, downsource                        [2, n]     UPSOURCE(N, 2), DOWNSOURCE(N, 2)
 *   intoprad, inbottomrad                       [n]        INTOPRAD(N), INBOTTOMRAD(N)
 *   uprad, downrad                              [n]        UPRAD(N), DOWNRAD(N), output
 */
void internal_radiance(ConstTensor3View upreflect,
                       ConstTensor3View uptrans,
                       ConstMatrixView  upsource,
                       ConstTensor3View downreflect,
                       ConstTensor3View downtrans,
                       ConstMatrixView  downsource,
                       ConstVectorView  intoprad,
                       ConstVectorView  inbottomrad,
                       VectorView       uprad,
                       VectorView       downrad,
                       workdata&        work);
}  // namespace polradtran
