#pragma once

#include <matpack.h>

#include "rt3.h"

/* The subroutines of 3rdparty/polradtran/radutil3.f, ported to C++ one at a
   time; RT3_LAMBERT_SURFACE and RT3_FRESNEL_SURFACE, which RT4 shares, are
   polradtran's (radutil.h).  Radiances are in SI, W m-2 Hz-1 sr-1, at the frequency in Hz (the
   Fortran's are per micrometre at the wavelength in micrometres).

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The counts are not passed; they are the extents of the
   arrays. */
namespace polradtran::rt3 {
/** THERMAL_RADIANCE (RT3_THERMAL_RADIANCE): the polarized radiance vector of
 * thermal emission at frequency for a body with albedo and temperature (K):
 * (1 - albedo) times ARTS's planck() in I, in the azimuth mode 0 only (the
 * emission is isotropic and unpolarized), the same for all mu.  Throws for
 * a negative temperature (the Fortran gave 0).
 *
 *   radiance  [2, nummu, nstokes]  RADIANCE(NSTOKES, NUMMU, 2), output
 */
void thermal_radiance(Index mode, Numeric temperature, Numeric albedo, Numeric frequency, Tensor3View radiance);

/** LAMBERT_RADIANCE (RT3_LAMBERT_RADIANCE): the ground radiance of a
 * Lambertian ground of albedo ground_albedo, in the azimuth mode 0 only:
 * with the thermal source (src_code 2 or 3) its emission,
 * (1 - ground_albedo) times ARTS's planck() at ground_temp [K] and the
 * frequency [Hz], and with the solar source (src_code 1 or 3) the direct
 * beam it reflects, direct_sfc_flux ground_albedo / pi, both in I and
 * unpolarized.  direct_sfc_flux is the direct flux on the ground, on the
 * horizontal, in W m-2 Hz-1; the radiance is in W m-2 Hz-1 sr-1.  With the
 * thermal source, it throws for a negative ground_temp (the Fortran gave
 * 0) and needs a positive frequency.
 *
 *   radiance  [nummu, nstokes]  RADIANCE(NSTOKES, NUMMU), output
 */
void lambert_radiance(Index      mode,
                      Index      src_code,
                      Numeric    ground_albedo,
                      Numeric    ground_temp,
                      Numeric    frequency,
                      Numeric    direct_sfc_flux,
                      MatrixView radiance);

/** FRESNEL_RADIANCE (RT3_FRESNEL_RADIANCE): the thermal radiance of a plane
 * ground of complex refractive index index under a medium of index 1, in
 * the azimuth mode 0 only: (1 - R) B for the unpolarized B, ARTS's
 * planck() at ground_temp [K] and the frequency [Hz], and the Fresnel
 * reflection R of fresnel_surface_layer, i.e. [(1 - R1) B, -R2 B, 0, 0],
 * in W m-2 Hz-1 sr-1.  It cannot reflect the direct beam (specularly).  In
 * mode 0 it throws for a negative ground_temp (the Fortran gave 0) and
 * needs a positive frequency.
 *
 *   mu_values  [nummu]           MU_VALUES
 *   radiance   [nummu, nstokes]  RADIANCE(NSTOKES, NUMMU), output
 */
void fresnel_radiance(
    Index mode, ConstVectorView mu_values, Complex index, Numeric ground_temp, Numeric frequency, MatrixView radiance);

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
 * rt3::lambertian_surface and rt3::fresnel_surface are the types of
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
