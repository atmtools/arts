#pragma once

#include <matpack.h>

/* The subroutines of 3rdparty/polradtran/radintg4.f that RT3 does not
   share, ported to C++ one at a time: INITIALIZE and INITIAL_SOURCE, which
   make the thin initial sublayer of the doubling from RT4's optics sets.
   The others (DOUBLING_INTEGRATION, COMBINE_LAYERS, INTERNAL_RADIANCE and
   NONSCATTER_LAYER) are polradtran's (radintg.h).  Each follows its Fortran
   step by step.  A Fortran array A(d1, ..., dk) is the row-major matpack
   array [dk, ..., d1].  The counts (NSTOKES, NUMMU) are not passed; they
   are the extents of the arrays. */
namespace polradtran::rt4 {
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
}  // namespace polradtran::rt4
