#pragma once

#include <matpack.h>

#include "rt3.h"

/* The ground of 3rdparty/polradtran/radutil3.f as rt3::radtran's input.
   Its routines (RT3_LAMBERT_SURFACE, RT3_LAMBERT_RADIANCE,
   RT3_FRESNEL_SURFACE, RT3_FRESNEL_RADIANCE and RT3_THERMAL_RADIANCE) are
   RT4's too, and polradtran's (radutil.h).  Radiances are in SI,
   W m-2 Hz-1 sr-1, at the frequency in Hz (the Fortran's are per
   micrometre at the wavelength in micrometres).

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The counts are not passed; they are the extents of the
   arrays. */
namespace polradtran::rt3 {
/** The ground as rt3::radtran's input, for every azimuth mode.
 *
 * Both grounds of RT3 (LAMBERT_SURFACE, FRESNEL_SURFACE) make the same
 * surface layer, polradtran::external_surface_layer's: only the reflection back up,
 * REFLECT(..., 2), depends on the ground and the mode.  With the ground's
 * radiance (LAMBERT_RADIANCE, FRESNEL_RADIANCE), that reflection is all a
 * ground is, and rt3::radtran takes the two as data.  LAMBERT_RADIANCE also
 * adds the direct beam that the ground reflects, DIRECT_SFC_FLUX A / pi.
 * The direct flux that reaches the ground is radtran's, so that part is
 * given per unit of direct flux, as direct_reflect, and radtran scales it.
 * The ground routines are ported (LAMBERT_SURFACE and FRESNEL_SURFACE as
 * polradtran::lambert_surface_layer and fresnel_surface_layer, since
 * polradtran::lambertian_surface and polradtran::fresnel_surface are the types of
 * rt3.h), so this calls no Fortran.
 *
 *   mu_values       [nummu]           MU_VALUES, the streams of rt3::radtran
 *   quad_weights    [nummu]           their weights, 0 for the extra angles
 *   surf_reflect    [aziorder + 1, nummu, nstokes, nummu, nstokes]
 *                                     REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2) of the ground in each mode,
 *                                     output
 *   gnd_radiance    [aziorder + 1, nummu, nstokes]
 *                                     GND_RADIANCE(NSTOKES, NUMMU) of each mode without the reflected direct
 *                                     beam, in W m-2 Hz-1 sr-1, output
 *   direct_reflect  [aziorder + 1, nummu, nstokes]
 *                                     the radiance reflected from the direct beam per unit of direct flux on
 *                                     the ground (on the horizontal), in sr-1, output
 *
 * src_code is RADTRAN's (0 none, 1 solar, 2 thermal, 3 both): a Lambertian
 * ground emits only with the thermal source, a Fresnel ground always (as
 * in RT3), and a Fresnel ground throws with the solar source, because RT3
 * cannot reflect the direct beam specularly.  ground_temp is in K and
 * frequency in Hz.
 */
void ground_surface(const surface&  ground,
                    Index           src_code,
                    ConstVectorView mu_values,
                    ConstVectorView quad_weights,
                    Numeric         frequency,
                    Numeric         ground_temp,
                    Tensor5View     surf_reflect,
                    Tensor3View     gnd_radiance,
                    Tensor3View     direct_reflect);
}  // namespace polradtran::rt3
