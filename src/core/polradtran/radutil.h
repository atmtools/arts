#pragma once

#include <matpack.h>

/* The ground routines that 3rdparty/polradtran/radutil3.f and radutil4.f
   share, ported to C++: RT3's RT3_X is RT4's X, line by line.  Each makes
   the ground as a surface layer for the adding: only its reflection back
   up, REFLECT(..., 2), depends on the ground; the reflection from above is
   0, the transmission the identity, and the source 0.  The X_SURFACE
   routines are X_surface_layer, since the grounds of rt3.h and rt4.h are
   the X_surface types.  MZERO is "= 0.0" and MIDENTITY identity.

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
}  // namespace polradtran
