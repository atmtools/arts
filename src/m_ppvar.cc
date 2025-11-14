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
    const ArrayOfAtmPoint &atm_point_path) try {
  ARTS_TIME_REPORT

  const Size np = atm_point_path.size();
  ARTS_USER_ERROR_IF(np != ray_path.size(),
                     "ray_path and atm_point_path must have the same size")

  ray_path_zeeman_magnetic_field.resize(np);

  std::transform(ray_path.begin(),
                 ray_path.end(),
                 atm_point_path.begin(),
                 ray_path_zeeman_magnetic_field.begin(),
                 [](const auto &p, const auto &a) -> Vector3 {
                   using Conversion::rad2deg;
                   const auto zz = lbl::zeeman::magnetic_angles(a.mag, p.los);
                   return {zz.H, rad2deg(zz.theta()), rad2deg(zz.eta())};
                 });
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_srcvec_pathFromPropmat(
    ArrayOfStokvecVector &spectral_rad_srcvec_path,
    ArrayOfStokvecMatrix &spectral_rad_srcvec_jac_path,
    const ArrayOfPropmatVector &spectral_propmat_path,
    const ArrayOfStokvecVector &ray_path_source_vector_nonlte,
    const ArrayOfPropmatMatrix &spectral_propmat_jac_path,
    const ArrayOfStokvecMatrix &ray_path_source_vector_nonlte_jacobian,
    const ArrayOfAscendingGrid &freq_grid_path,
    const ArrayOfAtmPoint &atm_point_path,
    const JacobianTargets &jac_targets) try {
  ARTS_TIME_REPORT

  String error{};

  const Index np = atm_point_path.size();
  if (np == 0) {
    spectral_rad_srcvec_path.resize(0);
    spectral_rad_srcvec_jac_path.resize(0);
    return;
  }

  const Index nf = spectral_propmat_path.front().size();
  const Index nq = jac_targets.target_count();

  const Index it = jac_targets.target_position(AtmKey::t);

  spectral_rad_srcvec_path.resize(np);
  for (auto &t : spectral_rad_srcvec_path) {
    t.resize(nf);
    t = 0;
  }
  spectral_rad_srcvec_jac_path.resize(np);
  for (auto &t : spectral_rad_srcvec_jac_path) {
    t.resize(nq, nf);
    t = 0;
  }

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Index ip = 0; ip < np; ip++) {
    try {
      rtepack::source::level_nlte(spectral_rad_srcvec_path[ip],
                                  spectral_rad_srcvec_jac_path[ip],
                                  spectral_propmat_path[ip],
                                  ray_path_source_vector_nonlte[ip],
                                  spectral_propmat_jac_path[ip],
                                  ray_path_source_vector_nonlte_jacobian[ip],
                                  freq_grid_path[ip],
                                  atm_point_path[ip].temperature,
                                  it);
    } catch (const std::runtime_error &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void atm_point_pathFromPath(ArrayOfAtmPoint &atm_point_path,
                            const ArrayOfPropagationPathPoint &ray_path,
                            const AtmField &atm_field) try {
  ARTS_TIME_REPORT

  forward_atm_path(
      atm_path_resize(atm_point_path, ray_path), ray_path, atm_field);
}
ARTS_METHOD_ERROR_CATCH

void freq_grid_pathFromPath(ArrayOfAscendingGrid &freq_grid_path,
                            ArrayOfVector3 &freq_wind_shift_jac_path,
                            const AscendingGrid &freq_grid,
                            const ArrayOfPropagationPathPoint &ray_path,
                            const ArrayOfAtmPoint &atm_point_path) try {
  ARTS_TIME_REPORT

  std::string error;

  freq_grid_path.resize(ray_path.size());
  freq_wind_shift_jac_path.resize(ray_path.size());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size ip = 0; ip < ray_path.size(); ip++) {
    try {
      freq_gridWindShift(freq_grid_path[ip] = freq_grid,
                         freq_wind_shift_jac_path[ip],
                         atm_point_path[ip],
                         ray_path[ip]);
    } catch (std::exception &e) {
#pragma omp critical
      error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void spectral_tramat_cumulative_pathFromPath(
    ArrayOfMuelmatVector &spectral_tramat_cumulative_path,
    const ArrayOfMuelmatVector &spectral_tramat_path) try {
  ARTS_TIME_REPORT

  forward_cumulative_transmission(spectral_tramat_cumulative_path,
                                  spectral_tramat_path);
}
ARTS_METHOD_ERROR_CATCH
