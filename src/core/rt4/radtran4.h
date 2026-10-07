#pragma once

#include <matpack.h>

#include "rt4.h"

namespace rt4 {
/** RADTRANO of 3rdparty/polradtran/radtran4.f, ported to C++.
 *
 * The method is RADTRANO's, step by step, in SI units, and every
 * subroutine it calls is ported too: NONSCATTER_LAYER, INITIAL_SOURCE,
 * INITIALIZE, DOUBLING_INTEGRATION, COMBINE_LAYERS and INTERNAL_RADIANCE
 * (radintg4.h), EXTERNAL_SURFACE and THERMAL_RADIANCE (radutil4.h).  The
 * quadratures are ARTS's (rt4::get_quadrature), the Planck function is
 * ARTS's planck(), MZERO is "= 0.0", MIDENTITY matpack::identity and MCOPY
 * "=".  It calls no Fortran and keeps no state between calls.  The STOPs
 * of RADTRANO throw instead.
 *
 * The arguments are RADTRANO's, except that the ground is external data
 * (GROUND_TYPE 'A', made for every kind of ground by rt4::ground_surface,
 * so GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX and GROUND_REFLEC
 * are gone) and that the counts are not passed:
 * NSTOKES, NUMMU, NUUMMU, NUM_LAYERS and NSL are the extents of the arrays
 * (nstokes of up_rad, nummu of mu_values, nuummu of extra_mu, num_layers of
 * height, nsl of extinct_matrix), and every other extent must agree.  A
 * Fortran array A(d1, ..., dk) is the row-major matpack array [dk, ..., d1]:
 *
 *   surf_reflect   [nummu, nstokes, nummu, nstokes]  SURF_REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU)
 *   gnd_radiance   [nummu, nstokes]                  GND_RADIANCE(NSTOKES, NUMMU)
 *   height         [num_layers + 1]                  HEIGHT
 *   temperatures   [num_layers + 1]                  TEMPERATURES
 *   gas_extinct    [num_layers]                      GAS_EXTINCT; clipped at 0 in place
 *   scatlayers     [num_layers]                      SCATLAYERS, the 1-based optics set of each
 *                                                    layer, < 1 for a non-scattering layer
 *   extinct_matrix [nsl, 2, nummu, nstokes, nstokes] EXTINCT_MATRIX(NSTOKES, NSTOKES, NUMMU, 2, NSL)
 *   emis_vector    [nsl, 2, nummu, nstokes]          EMIS_VECTOR(NSTOKES, NUMMU, 2, NSL)
 *   scatter_matrix [nsl, 4, nummu, nstokes,          SCATTER_MATRIX(NSTOKES, NUMMU, NSTOKES, NUMMU,
 *                   nummu, nstokes]                    4, NSL)
 *   extra_mu       [nuummu]                          the extra, zero-weight angles, input
 *   mu_values      [nummu]                           MU_VALUES, output: the nummu - nuummu
 *                                                    quadrature nodes followed by extra_mu
 *   up_rad         [num_layers + 1, nummu, nstokes]  UP_RAD(NSTOKES, NUMMU, NUM_LAYERS + 1)
 *   down_rad       [num_layers + 1, nummu, nstokes]  DOWN_RAD(NSTOKES, NUMMU, NUM_LAYERS + 1)
 *
 * quad_type is the rule of rt4::get_quadrature (QUAD_TYPE 'D', 'G' or
 * 'L').
 * frequency is in Hz, and the radiances (gnd_radiance, up_rad, down_rad)
 * are in W m-2 Hz-1 sr-1 (RADTRANO: WAVELENGTH in um and W m-2 sr-1 um-1).
 * The temperatures, sky_temp too, must be >= 0 K.
 */
void radtrano(Numeric          max_delta_tau,
              quadrature_type  quad_type,
              ConstTensor4View surf_reflect,
              ConstMatrixView  gnd_radiance,
              Numeric          sky_temp,
              Numeric          frequency,
              ConstVectorView  height,
              ConstVectorView  temperatures,
              VectorView       gas_extinct,
              ConstVectorView  scatlayers,
              ConstTensor5View extinct_matrix,
              ConstTensor4View emis_vector,
              ConstTensor6View scatter_matrix,
              ConstVectorView  extra_mu,
              VectorView       mu_values,
              Tensor3View      up_rad,
              Tensor3View      down_rad);
}  // namespace rt4
