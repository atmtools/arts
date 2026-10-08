#pragma once

#include <matpack.h>

#include "rt3_workdata.h"

/* The subroutines of 3rdparty/polradtran/radintg3.f, ported to C++ one at a
   time; RT3_COMBINE_LAYERS and RT3_INTERNAL_RADIANCE, which RT4 shares, are
   polradtran's (radintg.h).  The reflection and transmission matrices of a slab are the
   Fortran's (NSTOKES, NUMMU, NSTOKES, NUMMU, 2): [2, nummu, nstokes, nummu,
   nstokes], whose [l, j2, i2, j1, i1] is the column-major n x n matrix l
   (n = nstokes nummu) at row (i1, j1) and column (i2, j2), and its sources
   and radiances are (NSTOKES, NUMMU, 2): [2, nummu, nstokes].  l = 0 is
   the + (downwelling for reflection, see RADTRAN) and l = 1 the - part.
   The counts are not passed; they are the extents of the arrays. */
namespace polradtran::rt3 {
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
 * LAPACK's inv_inplace.  The scratch (X, Y, GAMMA, T_EXP, T_LIN, CONST,
 * T_CONST and two vectors) is work's, which must be sized for the n
 * streams (rt3_workdata::resize).
 */
void doubling_integration(Index         num_doubles,
                          Index         src_code,
                          bool          symmetric,
                          Tensor3View   reflect,
                          Tensor3View   trans,
                          MatrixView    exp_source,
                          Numeric       expfactor,
                          MatrixView    lin_source,
                          Numeric       linfactor,
                          Tensor3View   t_reflect,
                          Tensor3View   t_trans,
                          MatrixView    t_source,
                          rt3_workdata& work);

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
}  // namespace polradtran::rt3
