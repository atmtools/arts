#pragma once

#include <matpack.h>

#include "rt4.h"

/* The ground of 3rdparty/polradtran/radutil4.f, as the external surface of
   RADTRANO, and its routines ported to C++ one at a time.  Those RT3
   shares (LAMBERT_SURFACE, LAMBERT_RADIANCE, FRESNEL_SURFACE,
   FRESNEL_RADIANCE, EXTERNAL_SURFACE and THERMAL_RADIANCE) are
   polradtran's (radutil.h); SPECULAR_SURFACE and SPECULAR_RADIANCE are
   RT4's own.  Each follows its Fortran step by step, with matpack in place
   of Evans' matrix helpers (MZERO is "= 0.0", MIDENTITY a unit diagonal)
   and ARTS's planck() in SI in place of PLANCK_FUNCTION.  A Fortran array
   A(d1, ..., dk) is the row-major matpack array [dk, ..., d1]; the counts
   (NSTOKES, NUMMU) are the extents of the arrays. */
namespace polradtran::rt4 {
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
 * adding; polradtran::fresnel_surface and rt4::specular_surface are the types of
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
