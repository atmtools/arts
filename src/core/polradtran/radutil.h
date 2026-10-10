#pragma once

#include <matpack.h>

/* The ground routines that 3rdparty/polradtran/radutil3.f and radutil4.f
   share, ported to C++: RT3's RT3_X is RT4's X, line by line, or RT4's X
   is RT3's in the azimuth mode 0 with the thermal source alone.  The
   X_SURFACE routines make the ground as a surface layer for the adding:
   only its reflection back up, REFLECT(..., 2), depends on the ground; the
   reflection from above is 0, the transmission the identity, and the
   source 0.  They are X_surface_layer, since the grounds of rt3.h and
   rt4.h are the X_surface types.  The X_RADIANCE routines give the
   radiance of a ground or of the sky, in SI, W m-2 Hz-1 sr-1, at the
   frequency in Hz (the Fortran's are per micrometre at the wavelength in
   micrometres), with ARTS's planck().  MZERO is "= 0.0" and MIDENTITY
   identity.

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The counts (NSTOKES, NUMMU) are not passed; they are the
   extents of the arrays. */
namespace polradtran {
/** LAMBERT_SURFACE (RT3_LAMBERT_SURFACE): the surface layer of a
 * Lambertian ground of albedo ground_albedo, which reflects the flux
 * equally into all directions and completely unpolarizes the radiation, in
 * the azimuth mode 0 (2 ground_albedo mu_j w_j from stream j into every
 * stream, I to I only; any other mode reflects nothing).
 *
 *   mu_values     [nummu]                              MU_VALUES
 *   quad_weights  [nummu]                              QUAD_WEIGHTS
 *   reflect       [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans         [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source        [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output
 */
void lambert_surface_layer(Index           mode,
                           ConstVectorView mu_values,
                           ConstVectorView quad_weights,
                           Numeric         ground_albedo,
                           Tensor5View     reflect,
                           Tensor5View     trans,
                           Tensor3View     source);

/** FRESNEL_SURFACE (RT3_FRESNEL_SURFACE): the surface layer of a plane
 * ground of complex refractive index index under a medium of index 1: the
 * Fresnel reflection of each stream into itself, with ARTS's fresnel()
 * amplitudes and rtepack::fresnel_reflectance as the Mueller matrix
 * R = [[R1, R2, 0, 0], [R2, R1, 0, 0], [0, 0, R3, -R4], [0, 0, R4, R3]].
 * It is the same in every azimuth mode.
 *
 *   mu_values  [nummu]                              MU_VALUES
 *   reflect    [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans      [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source     [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output
 */
void fresnel_surface_layer(
    ConstVectorView mu_values, Complex index, Tensor5View reflect, Tensor5View trans, Tensor3View source);

/** EXTERNAL_SURFACE: the surface layer of a ground given by its reflection
 * back up, surf_reflect.  RT4 also passes the ground's RADIANCE, which
 * EXTERNAL_SURFACE does not use (RADTRANO gives it to INTERNAL_RADIANCE as
 * the radiance from below), so it is not an argument here.  RT3 has no
 * EXTERNAL_SURFACE; rt3::radtran takes its ground as this layer, too.
 *
 *   surf_reflect  [nummu, nstokes, nummu, nstokes]     SURF_REFL(NSTOKES, NUMMU, NSTOKES, NUMMU)
 *   reflect       [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans         [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source        [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output
 */
void external_surface_layer(ConstTensor4View surf_reflect, Tensor5View reflect, Tensor5View trans, Tensor3View source);

/** THERMAL_RADIANCE (RT3_THERMAL_RADIANCE): the polarized radiance vector of
 * thermal emission at frequency for a body with albedo and temperature (K):
 * (1 - albedo) times ARTS's planck() in I, in the azimuth mode 0 only (the
 * emission is isotropic and unpolarized), the same for all mu.  Throws for
 * a negative temperature (the Fortran gave 0).  RT4's THERMAL_RADIANCE is
 * RT3's in mode 0.
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
 * 0) and needs a positive frequency.  RT4's LAMBERT_RADIANCE is RT3's in
 * mode 0 with the thermal source alone.
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
 * needs a positive frequency.  RT4's FRESNEL_RADIANCE is RT3's in mode 0.
 *
 *   mu_values  [nummu]           MU_VALUES
 *   radiance   [nummu, nstokes]  RADIANCE(NSTOKES, NUMMU), output
 */
void fresnel_radiance(
    Index mode, ConstVectorView mu_values, Complex index, Numeric ground_temp, Numeric frequency, MatrixView radiance);
}  // namespace polradtran
