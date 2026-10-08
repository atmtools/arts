#pragma once

#include <matpack.h>

#include "rt4.h"

/* The ground of 3rdparty/polradtran/radutil4.f, as the external surface of
   RADTRANO, and its routines ported to C++ one at a time.  Each follows its
   Fortran step by step, with matpack in place of Evans' matrix helpers
   (MZERO is "= 0.0", MIDENTITY a unit diagonal) and ARTS's planck() in SI
   in place of PLANCK_FUNCTION.  A Fortran array A(d1, ..., dk) is the
   row-major matpack array [dk, ..., d1]; the counts (NSTOKES, NUMMU) are
   the extents of the arrays. */
namespace polradtran::rt4 {
/** LAMBERT_SURFACE: the reflection matrix of a Lambertian ground of albedo
 * ground_albedo, which reflects the flux equally into all directions and
 * completely unpolarizes the radiation (for mode 0; any other mode
 * reflects nothing), the identity as transmission, and no source.
 *
 *   mu_values     [nummu]  MU_VALUES
 *   quad_weights  [nummu]  QUAD_WEIGHTS
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

/** LAMBERT_RADIANCE: the thermal radiance of a Lambertian ground,
 * (1 - ground_albedo) B(ground_temp) in I, in W m-2 Hz-1 sr-1 at the
 * frequency [Hz].  ground_temp [K] must be >= 0 (planck() is negative
 * below 0 K, where PLANCK_FUNCTION gave 0).
 *
 *   radiance  [nummu, nstokes]  RADIANCE(NSTOKES, NUMMU), output
 */
void lambert_radiance(Numeric ground_albedo, Numeric ground_temp, Numeric frequency, MatrixView radiance);

/** FRESNEL_SURFACE: the reflection matrix of a plane surface of complex
 * refractive index index under a medium of index 1, from the Fresnel
 * formulae (ARTS's fresnel() amplitudes and rtepack::fresnel_reflectance),
 * the identity as transmission, and no source.
 *
 *   mu_values  [nummu]                              MU_VALUES
 *   reflect    [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans      [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source     [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output
 */
void fresnel_surface_layer(
    ConstVectorView mu_values, Complex index, Tensor5View reflect, Tensor5View trans, Tensor3View source);

/** FRESNEL_RADIANCE: the thermal radiance of a plane surface of complex
 * refractive index index under a medium of index 1, (1 - R) B(ground_temp)
 * for the unpolarized B and the Fresnel reflection matrix R, in
 * W m-2 Hz-1 sr-1 at the frequency [Hz].  ground_temp [K] must be >= 0.
 *
 *   mu_values  [nummu]           MU_VALUES
 *   radiance   [nummu, nstokes]  RADIANCE(NSTOKES, NUMMU), output
 */
void fresnel_radiance(
    ConstVectorView mu_values, Complex index, Numeric ground_temp, Numeric frequency, MatrixView radiance);

/** SPECULAR_SURFACE: the reflection matrix of a plane surface with the
 * fixed reflectivity ground_reflec, applied specularly to every stream, the
 * identity as transmission, and no source.
 *
 *   ground_reflec  [nstokes, nstokes]                   R(out, in): GROUND_REFLEC(NSTOKES, NSTOKES)
 *                                                       read row-major, as SPECULAR_SURFACE reads it
 *   reflect        [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans          [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source         [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output
 */
void specular_surface_layer(ConstMatrixView ground_reflec, Tensor5View reflect, Tensor5View trans, Tensor3View source);

/** SPECULAR_RADIANCE: the thermal radiance of a plane surface with the
 * fixed reflectivity ground_reflec, (1 - R) B(ground_temp) for the
 * unpolarized B, in W m-2 Hz-1 sr-1 at the frequency [Hz].  This is RT4's
 * [(1 - R(I, I)) B, -R(Q, I) B]; RT4 gave 0 for U and V, this
 * -R(U, I) B and -R(V, I) B.  ground_temp [K] must be >= 0.
 *
 *   ground_reflec  [nstokes, nstokes]  R(out, in), as for specular_surface_layer
 *   radiance       [nummu, nstokes]    RADIANCE(NSTOKES, NUMMU), output
 */
void specular_radiance(ConstMatrixView ground_reflec, Numeric ground_temp, Numeric frequency, MatrixView radiance);

/** EXTERNAL_SURFACE: the surface layer of a ground given by its reflection
 * back up, SURF_REFL: no reflection from above, SURF_REFL as the reflection
 * back up, the identity as transmission, and no source.  RT4 also passes
 * the ground's RADIANCE, which EXTERNAL_SURFACE does not use (RADTRANO
 * gives it to INTERNAL_RADIANCE as the radiance from below), so it is not
 * an argument here.
 *
 *   surf_reflect  [nummu, nstokes, nummu, nstokes]     SURF_REFL(NSTOKES, NUMMU, NSTOKES, NUMMU)
 *   reflect       [2, nummu, nstokes, nummu, nstokes]  REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   trans         [2, nummu, nstokes, nummu, nstokes]  TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2), output
 *   source        [2, nummu, nstokes]                  SOURCE(NSTOKES, NUMMU, 2), output
 */
void external_surface_layer(ConstTensor4View surf_reflect, Tensor5View reflect, Tensor5View trans, Tensor3View source);

/** THERMAL_RADIANCE: the polarized radiance vector of the thermal emission
 * of a body of albedo albedo at temperature [K], (1 - albedo) B(temperature)
 * in I in both hemispheres, in W m-2 Hz-1 sr-1 at the frequency [Hz].  The
 * emission is isotropic and unpolarized.  temperature must be >= 0
 * (planck() is negative below 0 K, where PLANCK_FUNCTION gave 0).
 *
 *   radiance  [2, nummu, nstokes]  RADIANCE(NSTOKES, NUMMU, 2), output
 */
void thermal_radiance(Numeric temperature, Numeric albedo, Numeric frequency, Tensor3View radiance);

/** The ground as RADTRANO's external surface (ground type 'A').
 *
 * Every ground routine of RT4 (LAMBERT_SURFACE, FRESNEL_SURFACE,
 * SPECULAR_SURFACE, EXTERNAL_SURFACE) makes the same surface layer: no
 * reflection from above, the identity as transmission, no source, and only
 * the reflection back up, REFLECT(..., 2), depending on the ground.  With
 * the ground's own radiance (LAMBERT_, FRESNEL_ or SPECULAR_RADIANCE), that
 * reflection is all a ground is, and RADTRANO takes the two as SURF_REFLECT
 * and GND_RADIANCE.  The ground routines are ported (the X_SURFACE
 * routines as X_surface_layer, which make the ground as a layer for the
 * adding; rt4::fresnel_surface and rt4::specular_surface are the types of
 * rt4.h), so this calls no Fortran.
 */
void ground_surface(const surface&  ground,
                    ConstVectorView mu_values,
                    ConstVectorView quad_weights,
                    Numeric         frequency,
                    Numeric         ground_temp,
                    Tensor4View     surf_reflect,
                    MatrixView      gnd_radiance);
}  // namespace polradtran::rt4
