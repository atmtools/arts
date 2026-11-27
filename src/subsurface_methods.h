#pragma once

#include <matpack.h>
#include <path_point.h>
#include <subsurf.h>
#include <surf.h>

/*===========================================================================
  === Functions shared by the subsurface methods
  ===========================================================================*/

/** The spectral point of the subsurface at a point of a path
 *
 * The subsurface fields are functions of the altitude relative to the local
 * surface, so the point's altitude is taken relative to the elevation of
 * surf_field at its latitude and longitude.
 *
 * @param subsurf_spectral_field The spectral subsurface field
 * @param surf_field The surface field, for the elevation
 * @param freq_grid The frequencies of the spectral point [Hz]
 * @param point The point of the path
 * @return The spectral point of the subsurface at the point
 */
[[nodiscard]] SubsurfaceSpectralPoint subsurface_spectral_point(const SubsurfaceSpectralField& subsurf_spectral_field,
                                                                const SurfaceField&            surf_field,
                                                                const AscendingGrid&           freq_grid,
                                                                const PropagationPathPoint&    point);

/** The complex refractive index of the subsurface at each frequency of a spectral point
 *
 * From its "refractive index" n' = sqrt(eps') and "loss tangent" tan(delta),
 * see dielectric_refractive_index.
 *
 * @param point The spectral point, with nf values of each property
 * @param nf The number of frequencies
 * @return The complex refractive index at each frequency, Im(n) >= 0
 */
[[nodiscard]] ComplexVector subsurface_refractive_index(const SubsurfaceSpectralPoint& point, Size nf);
