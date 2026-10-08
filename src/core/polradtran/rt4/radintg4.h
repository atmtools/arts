#pragma once

#include <matpack.h>
#include <polradtran_workdata.h>

/* The subroutines of 3rdparty/polradtran/radintg4.f, ported to C++ one at a
   time; COMBINE_LAYERS and INTERNAL_RADIANCE, which RT3 shares, are
   polradtran's (radintg.h).  Each follows its Fortran step by step, with matpack in place of
   Evans' matrix helpers (radmat.f): MZERO is "= 0.0", MCOPY "=", MADD and
   MSCALARMULT "+=" and "*=", MINVERT is inv_inplace (LAPACK), and MMULT is
   mult (BLAS DGEMM, or DGEMV for a matrix-vector product), whose alpha and
   beta absorb the MIDENTITY, MSUB and MADD around a product.  The scratch of
   DOUBLING_INTEGRATION (X, Y, GAMMA and their vectors) is that of a
   workdata, which must be sized for its n streams (workdata::resize).

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The counts (NSTOKES, NUMMU) are not passed; they are the
   extents of the arrays. */
namespace polradtran::rt4 {
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
                          MatrixView  t_source,
                          workdata&   work);

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
}  // namespace polradtran::rt4
