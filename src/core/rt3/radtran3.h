#pragma once

#include <matpack.h>

#include "rt3.h"

namespace rt3 {
/** RADTRAN of 3rdparty/polradtran/radtran3.f, ported to C++.
 *
 * The method is RADTRAN's, step by step, calling the same subroutines,
 * all ported: those of radscat3.f (radscat3.h), THERMAL_RADIANCE
 * (radutil3.h), NONSCATTER_LAYER, INITIAL_SOURCE, INITIALIZE,
 * DOUBLING_INTEGRATION, COMBINE_LAYERS and INTERNAL_RADIANCE (radintg3.h),
 * and the grounds of radutil3.f behind rt3::ground_surface.  It calls no
 * Fortran and keeps no state between calls.  The quadratures are ARTS's
 * (rt3::get_quadrature), the Planck function is ARTS's planck(), and
 * Evans' matrix helpers are matpack: MZERO is "= 0.0", MIDENTITY
 * matpack::identity, MCOPY "=" and MSCALARMULT "*=".  The STOPs of RADTRAN
 * throw instead.
 *
 * The arguments are RADTRAN's, except that the ground is external data,
 * that the extra, zero-weight angles are an input of their own (RADTRAN's
 * QUAD_TYPE 'E' marked them by non-zero entries of MU_VALUES), and that
 * the counts are not passed: NSTOKES, NUMMU, AZIORDER,
 * NUM_LAYERS, NSL, LDCOEF and NOUTLEVELS are the extents of the arrays
 * (nstokes and aziorder + 1 of up_rad, nummu of mu_values, num_layers of
 * height, nsl and ldcoef of scat_coef, noutlevels of outlevels).
 *
 * The ground replaces GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO and
 * GROUND_INDEX.  Every ground of RT3 is a surface layer that reflects back
 * up only, with a radiance of its own and, with the solar source, the
 * direct beam it reflects; rt3::ground_surface makes the three arrays for
 * the Lambertian and the Fresnel ground.  Each azimuth mode has its own
 * reflection and radiance, and RADTRAN's GND_RADIANCE of a mode is
 * gnd_radiance + F direct_reflect, F the direct flux that reaches the
 * ground (on the horizontal).
 *
 * A Fortran array A(d1, ..., dk) is the row-major matpack array
 * [dk, ..., d1]:
 *
 *   surf_reflect    [aziorder + 1, nummu, nstokes, nummu, nstokes]
 *                                          REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2) of the ground in each mode
 *   gnd_radiance    [aziorder + 1, nummu, nstokes]
 *                                          the ground's radiance in each mode, without the reflected direct beam
 *   direct_reflect  [aziorder + 1, nummu, nstokes]
 *                                          the radiance the ground reflects from the direct beam per unit of direct
 *                                          flux (sr-1) in each mode; used only with the solar source
 *   height          [num_layers + 1]       HEIGHT
 *   temperatures    [num_layers + 1]       TEMPERATURES
 *   gas_extinct     [num_layers]           GAS_EXTINCT (clipped at 0 where used)
 *   scat_extinct    [nsl]                  SCAT_EXTINCT
 *   scat_scatter    [nsl]                  SCAT_SCATTER
 *   scat_nlegen     [nsl]                  SCAT_NLEGEN
 *   scat_coef       [nsl, ldcoef, 6]       SCAT_COEF(6, LDCOEF, NSL)
 *   scatlayers      [num_layers]           SCATLAYERS, the 1-based set of each layer, 0 for none
 *   outlevels       [noutlevels]           OUTLEVELS, 1 (top) to num_layers + 1 (bottom)
 *   extra_mu        [nuummu]               the extra angles, input; only with the gauss quadrature
 *   mu_values       [nummu]                MU_VALUES, output: the quadrature nodes, then extra_mu
 *   up_flux         [noutlevels, nstokes]  UP_FLUX(NSTOKES, NOUTLEVELS)
 *   down_flux       [noutlevels, nstokes]  DOWN_FLUX(NSTOKES, NOUTLEVELS)
 *   up_rad          [noutlevels, aziorder + 1, nummu, nstokes]
 *                                          UP_RAD(NSTOKES, NUMMU, AZIORDER+1, NOUTLEVELS)
 *   down_rad        [noutlevels, aziorder + 1, nummu, nstokes]
 *                                          DOWN_RAD(NSTOKES, NUMMU, AZIORDER+1, NOUTLEVELS)
 *
 * src_code: 0 none, 1 solar, 2 thermal, 3 both.  quad_type is the rule of
 * the quadrature nodes (QUAD_TYPE 'G', 'D' or 'L'; 'G' with extra_mu is
 * 'E').  delta_m is DELTAM = 'Y'.  The radiances are in SI
 * (W m-2 Hz-1 sr-1) at the frequency in Hz, as is direct_flux
 * (W m-2 Hz-1); RADTRAN took the wavelength in micrometres and worked per
 * micrometre.
 */
void radtran(Numeric             max_delta_tau,
             Index               src_code,
             quadrature_type     quad_type,
             bool                delta_m,
             Numeric             direct_flux,
             Numeric             direct_mu,
             ConstTensor5View    surf_reflect,
             ConstTensor3View    gnd_radiance,
             ConstTensor3View    direct_reflect,
             Numeric             sky_temp,
             Numeric             frequency,
             ConstVectorView     height,
             ConstVectorView     temperatures,
             ConstVectorView     gas_extinct,
             ConstVectorView     scat_extinct,
             ConstVectorView     scat_scatter,
             const ArrayOfIndex& scat_nlegen,
             ConstTensor3View    scat_coef,
             const ArrayOfIndex& scatlayers,
             const ArrayOfIndex& outlevels,
             ConstVectorView     extra_mu,
             VectorView          mu_values,
             MatrixView          up_flux,
             MatrixView          down_flux,
             Tensor4View         up_rad,
             Tensor4View         down_rad);
}  // namespace rt3
