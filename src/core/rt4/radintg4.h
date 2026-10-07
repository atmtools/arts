#pragma once

#include <matpack.h>

/* The subroutines of 3rdparty/polradtran/radintg4.f, ported to C++ one at a
   time.  Each follows its Fortran step by step, with matpack in place of
   Evans' matrix helpers (radmat.f): MZERO is "= 0.0", MCOPY "=", MADD and
   MSCALARMULT "+=" and "*=", MINVERT is inv_inplace (LAPACK), and MMULT is
   mult (BLAS DGEMM, or DGEMV for a matrix-vector product), whose alpha and
   beta absorb the MIDENTITY, MSUB and MADD around a product.

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The counts (NSTOKES, NUMMU) are not passed; they are the
   extents of the arrays. */
namespace rt4 {
/** DOUBLING_INTEGRATION: integrates homogeneous thin layers with the
 * doubling algorithm, num_doubles doubling steps.  The initial reflection
 * and transmission matrices are input.  A linear (thermal) source is
 * assumed: lin_source is the source vector at zero optical depth and
 * linfactor the single-layer slope of the linear source.  With symmetric,
 * the minus parts of the reflection and transmission matrices are taken
 * to equal the plus parts instead of being computed.
 *
 *   reflect     [2, n, n]  REFLECT(N, N, 2), the + and - n x n matrices; overwritten
 *   trans       [2, n, n]  TRANS(N, N, 2); overwritten
 *   lin_source  [2, n]     LIN_SOURCE(N, 2); overwritten
 *   t_reflect   [2, n, n]  T_REFLECT(N, N, 2), output
 *   t_trans     [2, n, n]  T_TRANS(N, N, 2), output
 *   t_source    [2, n]     T_SOURCE(N, 2), output
 *
 * As matpack matrices, the Fortran's column-major n x n matrices hold their
 * transposes; the products are computed for the Fortran matrices (C = A B
 * is mult(C, B, A), y = A x is mult(y, transpose(A), x)).
 */
void doubling_integration(Index       num_doubles,
                          bool        symmetric,
                          Tensor3View reflect,
                          Tensor3View trans,
                          MatrixView  lin_source,
                          Numeric     linfactor,
                          Tensor3View t_reflect,
                          Tensor3View t_trans,
                          MatrixView  t_source);

/** COMBINE_LAYERS: the reflection and transmission matrices and the source
 * vectors of two layers combined into one.  The positive side (down) of
 * the first layer is attached to the negative side (up) of the second
 * layer; thus layer 1 is put on top of layer 2.
 *
 *   reflect1, trans1, reflect2, trans2  [2, n, n]  REFLECT1(N, N, 2), ...
 *   source1, source2                    [2, n]     SOURCE1(N, 2), SOURCE2(N, 2)
 *   out_reflect, out_trans              [2, n, n]  OUT_REFLECT(N, N, 2), OUT_TRANS(N, N, 2), output
 *   out_source                          [2, n]     OUT_SOURCE(N, 2), output
 *
 * The outputs must not share memory with the inputs.  As matpack matrices,
 * the Fortran's column-major n x n matrices hold their transposes, as for
 * doubling_integration.
 */
void combine_layers(ConstTensor3View reflect1,
                    ConstTensor3View trans1,
                    ConstMatrixView  source1,
                    ConstTensor3View reflect2,
                    ConstTensor3View trans2,
                    ConstMatrixView  source2,
                    Tensor3View      out_reflect,
                    Tensor3View      out_trans,
                    MatrixView       out_source);

/** INTERNAL_RADIANCE: the internal radiance at a level.  The reflection
 * and transmission matrices and source vectors are given for the
 * atmosphere above (up) and below (down) the level.  The upwelling and
 * downwelling radiances are computed from the two layer properties and the
 * radiance incident on the top and bottom.
 *
 *   upreflect, uptrans, downreflect, downtrans  [2, n, n]  UPREFLECT(N, N, 2), ...
 *   upsource, downsource                        [2, n]     UPSOURCE(N, 2), DOWNSOURCE(N, 2)
 *   intoprad, inbottomrad                       [n]        INTOPRAD(N), INBOTTOMRAD(N)
 *   uprad, downrad                              [n]        UPRAD(N), DOWNRAD(N), output
 *
 * As matpack matrices, the Fortran's column-major n x n matrices hold their
 * transposes, as for doubling_integration.
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
                       VectorView       downrad);

/** INITIALIZE: the reflection and transmission matrices of the initial,
 * thin sublayer of a scattering layer for the doubling, to first order in
 * its thickness.
 *
 *   delta_z         the thickness of the sublayer
 *   mu_values       [nummu]                              MU_VALUES
 *   quad_weights    [nummu]                              QUAD_WEIGHTS
 *   gas_extinct     the gas extinction, per unit length (on the Stokes diagonal)
 *   extinct_matrix  [2, nummu, nstokes, nstokes]         EXTINCT_MATRIX(NSTOKES, NSTOKES, NUMMU, 2),
 *                                                        per unit length
 *   scatter_matrix  [4, nummu, nstokes, nummu, nstokes]  SCATTER_MATRIX(NSTOKES, NUMMU, NSTOKES, NUMMU, 4),
 *                                                        per unit length and steradian
 *   reflect         [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans           [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 */
void initialize(Numeric          delta_z,
                ConstVectorView  mu_values,
                ConstVectorView  quad_weights,
                Numeric          gas_extinct,
                ConstTensor4View extinct_matrix,
                ConstTensor5View scatter_matrix,
                Tensor5View      reflect,
                Tensor5View      trans);

/** INITIAL_SOURCE: the source vectors of the initial, thin sublayer of a
 * scattering layer for the doubling, from the emission (absorption) vector
 * of the particles and the gas extinction (which emits in I only).
 *
 *   delta_z      the thickness of the sublayer
 *   mu_values    [nummu]               MU_VALUES
 *   planck       the Planck function
 *   emis_vector  [2, nummu, nstokes]   EMIS_VECTOR(NSTOKES, NUMMU, 2), per unit length
 *   gas_extinct  the gas extinction, per unit length
 *   source       [2, nummu, nstokes]   SOURCE(NSTOKES, NUMMU, 2), output, in the unit of
 *                                      the Planck function
 */
void initial_source(Numeric          delta_z,
                    ConstVectorView  mu_values,
                    Numeric          planck,
                    ConstTensor3View emis_vector,
                    Numeric          gas_extinct,
                    Tensor3View      source);

/** NONSCATTER_LAYER: the reflection and transmission matrices and the source
 * vectors of a purely absorbing layer.  The source function, the Planck
 * function, varies linearly with optical depth across the layer.
 *
 *   deltatau   the vertical optical thickness of the layer
 *   mu_values  [nummu]                           MU_VALUES
 *   planck0    the Planck function at the top of the layer
 *   planck1    the Planck function at the bottom of the layer
 *   reflect    [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans      [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source     [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output, in the unit of
 *                                                   the Planck function
 */
void nonscatter_layer(Numeric         deltatau,
                      ConstVectorView mu_values,
                      Numeric         planck0,
                      Numeric         planck1,
                      Tensor5View     reflect,
                      Tensor5View     trans,
                      Tensor3View     source);
}  // namespace rt4
