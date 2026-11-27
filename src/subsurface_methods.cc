#include "subsurface_methods.h"

#include <debug.h>
#include <physics_funcs.h>

SubsurfaceSpectralPoint subsurface_spectral_point(const SubsurfaceSpectralField& subsurf_spectral_field,
                                                  const SurfaceField&            surf_field,
                                                  const AscendingGrid&           freq_grid,
                                                  const PropagationPathPoint&    point) {
  const Numeric lat = point.latitude(), lon = point.longitude();
  const Numeric h = surf_field.single_value(SurfaceKey::h, lat, lon);
  return subsurf_spectral_field.at(freq_grid, point.altitude() - h, lat, lon);
}

ComplexVector subsurface_refractive_index(const SubsurfaceSpectralPoint& point, const Size nf) {
  const Vector& refractive_index = point[SubsurfaceSpectralPropertyTag{"refractive index"}];
  const Vector& loss_tangent     = point[SubsurfaceSpectralPropertyTag{"loss tangent"}];

  if (refractive_index.size() != nf or loss_tangent.size() != nf)
    arts::user_error("The subsurface has {} refractive indices and {} loss tangents, but there are {} frequencies",
                     refractive_index.size(),
                     loss_tangent.size(),
                     nf);

  ComplexVector out(nf);
  for (Size iv = 0; iv < nf; iv++) out[iv] = dielectric_refractive_index(refractive_index[iv], loss_tangent[iv]);
  return out;
}
