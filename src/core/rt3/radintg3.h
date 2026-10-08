#pragma once

#include <matpack.h>

/* The subroutines of 3rdparty/polradtran/radintg3.f, ported to C++ one at a
   time.  The reflection and transmission matrices of a slab are the
   Fortran's (NSTOKES, NUMMU, NSTOKES, NUMMU, 2): [2, nummu, nstokes, nummu,
   nstokes], whose [l, j2, i2, j1, i1] is the column-major n x n matrix l
   (n = nstokes nummu) at row (i1, j1) and column (i2, j2), and its sources
   and radiances are (NSTOKES, NUMMU, 2): [2, nummu, nstokes].  l = 0 is
   the + (downwelling for reflection, see RADTRAN) and l = 1 the - part.
   The counts are not passed; they are the extents of the arrays. */
namespace rt3 {
/** DOUBLING_INTEGRATION (RT3_DOUBLING_INTEGRATION): integrates homogeneous
 * thin layers with the doubling algorithm, num_doubles doubling steps.  The
 * initial reflection and transmission matrices are input.  Depending on
 * src_code (1 solar, 2 thermal, 3 both, 0 none) the exponential (solar)
 * and linear (thermal) sources are doubled: exp_source and lin_source are
 * the source vectors at zero optical depth, expfactor is the single-layer
 * attenuation of the exponential source and linfactor the single-layer
 * slope of the linear one.  With symmetric, the minus parts of the
 * reflection and transmission matrices are taken to equal the plus parts
 * instead of being computed.
 *
 *   reflect     [2, n, n]  REFLECT(N, N, 2), the + and - n x n matrices; overwritten
 *   trans       [2, n, n]  TRANS(N, N, 2); overwritten
 *   exp_source  [2, n]     EXP_SOURCE(N, 2); overwritten
 *   lin_source  [2, n]     LIN_SOURCE(N, 2); overwritten
 *   t_reflect   [2, n, n]  T_REFLECT(N, N, 2), output
 *   t_trans     [2, n, n]  T_TRANS(N, N, 2), output
 *   t_source    [2, n]     T_SOURCE(N, 2), output
 *
 * As matpack matrices, the Fortran's column-major n x n matrices hold their
 * transposes; the products are computed for the Fortran matrices (C = A B
 * is mult(C, B, A), y = A x is mult(y, transpose(A), x)).  MINVERT is
 * LAPACK's inv_inplace.
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
                          MatrixView  t_source);

/** COMBINE_LAYERS (RT3_COMBINE_LAYERS): combines the reflection and
 * transmission matrices and source vectors of two layers into those of the
 * combined layer.  The positive side (down) of the first layer is attached
 * to the negative side (up) of the second layer; thus layer 1 is put on top
 * of layer 2.
 *
 *   reflect1, trans1, reflect2, trans2  [2, n, n]  REFLECT1(N, N, 2), ..., the + and - n x n matrices
 *   source1, source2                    [2, n]     SOURCE1(N, 2), SOURCE2(N, 2)
 *   out_reflect, out_trans              [2, n, n]  OUT_REFLECT(N, N, 2), OUT_TRANS(N, N, 2), output
 *   out_source                          [2, n]     OUT_SOURCE(N, 2), output
 *
 * The outputs must not overlap the inputs.  As matpack matrices, the
 * Fortran's column-major n x n matrices hold their transposes, as for
 * doubling_integration; MINVERT is LAPACK's inv_inplace.
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
 *
 * As matpack matrices, the Fortran's column-major n x n matrices hold their
 * transposes, as for doubling_integration; MINVERT is LAPACK's
 * inv_inplace.
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

/** INITIALIZE (RT3_INITIALIZE): infinitesimal generator initialization of
 * the local reflection and transmission matrices of a layer of thickness
 * delta_z from the phase function matrix, extinction and albedo.  With
 * f = delta_z / mu extinction for the angle of the row, the reflection is
 * f albedo times P+- (l = 0) and P-+ (l = 1), the transmission
 * 1 - f (1 - albedo P++) and 1 - f (1 - albedo P--).
 *
 *   mu_values       [nummu]
 *   phase_function  [4, nummu, nstokes, nummu, nstokes]  PHASE_FUNCTION(N, N, 4), the
 *                                                        SCATTER_MATRIX of get_scattering
 *   reflect         [2, nummu, nstokes, nummu, nstokes]  REFLECT(N, N, 2), output
 *   trans           [2, nummu, nstokes, nummu, nstokes]  TRANS(N, N, 2), output
 */
void initialize(Numeric          delta_z,
                ConstVectorView  mu_values,
                Numeric          extinction,
                Numeric          albedo,
                ConstTensor5View phase_function,
                Tensor5View      reflect,
                Tensor5View      trans);

/** INITIAL_SOURCE (RT3_INITIAL_SOURCE): infinitesimal generator
 * initialization of a source vector, for the thin layer of thickness
 * delta_z with which the doubling starts: delta_z / mu times extinction
 * times source_vector, for each angle.
 *
 *   mu_values      [nummu]
 *   source_vector  [2, nummu, nstokes]  SOURCE_VECTOR(N, 2)
 *   source         [2, nummu, nstokes]  SOURCE(N, 2), output
 */
void initial_source(
    Numeric delta_z, ConstVectorView mu_values, Numeric extinction, ConstTensor3View source_vector, Tensor3View source);

/** NONSCATTER_LAYER (RT3_NONSCATTER_LAYER): the reflection and transmission
 * matrices and the source vectors of a purely absorbing layer of optical
 * depth deltatau.  The source function, the Planck function (planck0 at the
 * top, planck1 at the bottom), varies linearly with optical depth across
 * the layer; it is set for mode 0 only, and for deltatau > 0.
 *
 *   mu_values  [nummu]
 *   reflect    [2, nummu, nstokes, nummu, nstokes]  REFLECT, output (0)
 *   trans      [2, nummu, nstokes, nummu, nstokes]  TRANS, output (diagonal exp(-deltatau / mu))
 *   source     [2, nummu, nstokes]                  SOURCE, output
 */
void nonscatter_layer(Index           mode,
                      Numeric         deltatau,
                      ConstVectorView mu_values,
                      Numeric         planck0,
                      Numeric         planck1,
                      Tensor5View     reflect,
                      Tensor5View     trans,
                      Tensor3View     source);
}  // namespace rt3
