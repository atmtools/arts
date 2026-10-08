#pragma once

#include <matpack.h>

#include "polradtran_workdata.h"

/* The subroutines that 3rdparty/polradtran/radintg3.f and radintg4.f share,
   ported to C++: RT3's RT3_X is RT4's X, line by line, or RT4's X is RT3's
   in the azimuth mode 0 with the thermal source alone.  Each follows its
   Fortran step by step, with matpack in place of Evans' matrix helpers
   (radmat.f): MZERO is "= 0.0", MCOPY "=", MADD and MSCALARMULT "+=" and
   "*=", MINVERT is inv_inplace (LAPACK), and MMULT is mult (BLAS DGEMM, or
   DGEMV for a matrix-vector product), whose alpha and beta absorb the
   MIDENTITY, MSUB and MADD around a product.  The scratch of
   DOUBLING_INTEGRATION, COMBINE_LAYERS and INTERNAL_RADIANCE (X, Y, GAMMA
   and their vectors) is that of a workdata, which must be sized for their
   n streams (workdata::resize).

   A slab's reflection and transmission matrices are [2, n, n], the
   Fortran's (N, N, 2): the column-major n x n matrices of the + (l = 0) and
   - (l = 1) parts, so a row-major matpack matrix holds the transpose of the
   Fortran matrix and the products are computed for the Fortran matrices
   (C = A B is mult(C, B, A), y = A x is mult(y, transpose(A), x)).  As
   (NSTOKES, NUMMU, NSTOKES, NUMMU, 2) they are [2, nummu, nstokes, nummu,
   nstokes], whose [l, j2, i2, j1, i1] is matrix l at row (i1, j1) and
   column (i2, j2).  Its sources and radiances are [2, n], the Fortran's
   (N, 2), or [2, nummu, nstokes].  The counts are not passed; they are the
   extents of the arrays. */
namespace polradtran {
/** DOUBLING_INTEGRATION (RT3_DOUBLING_INTEGRATION): integrates homogeneous
 * thin layers with the doubling algorithm, num_doubles doubling steps.  The
 * initial reflection and transmission matrices are input.  Depending on
 * src_code (1 solar, 2 thermal, 3 both, 0 none) the exponential (solar)
 * and linear (thermal) sources are doubled: exp_source and lin_source are
 * the source vectors at zero optical depth, expfactor is the single-layer
 * attenuation of the exponential source and linfactor the single-layer
 * slope of the linear one.  With symmetric, the minus parts of the
 * reflection and transmission matrices are taken to equal the plus parts
 * instead of being computed.  RT4's DOUBLING_INTEGRATION is RT3's with the
 * linear source alone (the overload below).
 *
 *   reflect     [2, n, n]  REFLECT(N, N, 2), the + and - n x n matrices; overwritten
 *   trans       [2, n, n]  TRANS(N, N, 2); overwritten
 *   exp_source  [2, n]     EXP_SOURCE(N, 2); overwritten; only with the solar source
 *   lin_source  [2, n]     LIN_SOURCE(N, 2); overwritten
 *   t_reflect   [2, n, n]  T_REFLECT(N, N, 2), output
 *   t_trans     [2, n, n]  T_TRANS(N, N, 2), output
 *   t_source    [2, n]     T_SOURCE(N, 2), output
 *
 * The outputs must not share memory with the inputs.  The scratch (X, Y,
 * GAMMA, T_LIN, CONST, T_CONST and two vectors) is work's; T_EXP is
 * t_source, which is written only after the doubling.
 */
void doubling_integration(Index       num_doubles,
                          Index       src_code,
                          bool        symmetric,
                          Tensor3View reflect,
                          Tensor3View trans,
                          MatrixView  exp_source,
                          Numeric     expfactor,
                          MatrixView  lin_source,
                          Numeric     linfactor,
                          Tensor3View t_reflect,
                          Tensor3View t_trans,
                          MatrixView  t_source,
                          workdata&   work);

/** DOUBLING_INTEGRATION of RT4: the above with the linear (thermal) source
 * alone (src_code 2, no exp_source).
 */
void doubling_integration(Index       num_doubles,
                          bool        symmetric,
                          Tensor3View reflect,
                          Tensor3View trans,
                          MatrixView  lin_source,
                          Numeric     linfactor,
                          Tensor3View t_reflect,
                          Tensor3View t_trans,
                          MatrixView  t_source,
                          workdata&   work);

/** NONSCATTER_LAYER (RT3_NONSCATTER_LAYER): the reflection and transmission
 * matrices and the source vectors of a purely absorbing layer of optical
 * depth deltatau.  The source function, the Planck function (planck0 at the
 * top, planck1 at the bottom), varies linearly with optical depth across
 * the layer; it is set in the azimuth mode 0 only, and for deltatau > 0.
 * RT4's NONSCATTER_LAYER is RT3's in mode 0.
 *
 *   mu_values  [nummu]
 *   reflect    [2, nummu, nstokes, nummu, nstokes]  REFLECT, output (0)
 *   trans      [2, nummu, nstokes, nummu, nstokes]  TRANS, output (diagonal exp(-deltatau / mu))
 *   source     [2, nummu, nstokes]                  SOURCE, output, in the unit of the Planck function
 */
void nonscatter_layer(Index           mode,
                      Numeric         deltatau,
                      ConstVectorView mu_values,
                      Numeric         planck0,
                      Numeric         planck1,
                      Tensor5View     reflect,
                      Tensor5View     trans,
                      Tensor3View     source);

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
