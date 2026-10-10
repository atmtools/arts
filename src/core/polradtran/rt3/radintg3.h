#pragma once

#include <matpack.h>

/* The subroutines of 3rdparty/polradtran/radintg3.f that RT4 does not
   share, ported to C++ one at a time: INITIALIZE and INITIAL_SOURCE, which
   make the thin initial sublayer of the doubling from RT3's scattering
   sets.  The others (RT3_DOUBLING_INTEGRATION, RT3_COMBINE_LAYERS,
   RT3_INTERNAL_RADIANCE and RT3_NONSCATTER_LAYER) are polradtran's
   (radintg.h).  The reflection and transmission matrices of a slab are the
   Fortran's (NSTOKES, NUMMU, NSTOKES, NUMMU, 2): [2, nummu, nstokes, nummu,
   nstokes], whose [l, j2, i2, j1, i1] is the column-major n x n matrix l
   (n = nstokes nummu) at row (i1, j1) and column (i2, j2), and its sources
   and radiances are (NSTOKES, NUMMU, 2): [2, nummu, nstokes].  l = 0 is
   the + (downwelling for reflection, see RADTRAN) and l = 1 the - part.
   The counts are not passed; they are the extents of the arrays. */
namespace polradtran::rt3 {
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
}  // namespace polradtran::rt3
