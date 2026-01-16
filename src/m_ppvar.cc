#include <array_algo.h>
#include <arts_conversions.h>
#include <arts_omp.h>
#include <atm_path.h>
#include <debug.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <stdexcept>

void ray_path_zeeman_magnetic_fieldFromPath(
    ArrayOfVector3 &ray_path_zeeman_magnetic_field,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &atm_path) try {
  ARTS_TIME_REPORT

  const Size np = atm_path.size();
  ARTS_USER_ERROR_IF(np != ray_path.size(),
                     "ray_path and atm_path must have the same size")

  ray_path_zeeman_magnetic_field.resize(np);

  std::transform(ray_path.begin(),
                 ray_path.end(),
                 atm_path.begin(),
                 ray_path_zeeman_magnetic_field.begin(),
                 [](const auto &p, const auto &a) -> Vector3 {
                   using Conversion::rad2deg;
                   const auto zz = lbl::zeeman::magnetic_angles(a.mag, p.los);
                   return {zz.H, rad2deg(zz.theta()), rad2deg(zz.eta())};
                 });
}
ARTS_METHOD_ERROR_CATCH

void atm_pathFromPath(ArrayOfAtmPoint &atm_path,
                      const ArrayOfPropagationPathPoint &ray_path,
                      const AtmField &atm_field) try {
  ARTS_TIME_REPORT

  forward_atm_path(atm_path_resize(atm_path, ray_path), ray_path, atm_field);
}
ARTS_METHOD_ERROR_CATCH

void freq_grid_pathFromPath(ArrayOfAscendingGrid &freq_grid_path,
                            ArrayOfVector3 &freq_wind_shift_jac_path,
                            const AscendingGrid &freq_grid,
                            const ArrayOfPropagationPathPoint &ray_path,
                            const ArrayOfAtmPoint &atm_path) try {
  ARTS_TIME_REPORT

  std::string error;

  freq_grid_path.resize(ray_path.size());
  freq_wind_shift_jac_path.resize(ray_path.size());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size ip = 0; ip < ray_path.size(); ip++) {
    try {
      freq_gridWindShift(freq_grid_path[ip] = freq_grid,
                         freq_wind_shift_jac_path[ip],
                         atm_path[ip],
                         ray_path[ip]);
    } catch (std::exception &e) {
#pragma omp critical
      error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH
