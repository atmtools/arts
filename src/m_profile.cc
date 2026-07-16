#include <array_algo.h>
#include <workspace.h>

#include <algorithm>
#include <ranges>

void Profile2Path(ArrayOfPropmatVector&              spectral_propmat_path,
                  ArrayOfStokvecVector&              spectral_nlte_srcvec_path,
                  ArrayOfPropmatMatrix&              spectral_propmat_jac_path,
                  ArrayOfStokvecMatrix&              spectral_nlte_srcvec_jac_path,
                  ArrayOfAtmPoint&                   atm_path,
                  ArrayOfAscendingGrid&              freq_grid_path,
                  ArrayOfVector3&                    freq_wind_shift_jac_path,
                  const AscendingGrid&               freq_grid,
                  const AscendingGrid&               alt_grid,
                  const ArrayOfPropagationPathPoint& ray_path,
                  const ArrayOfAtmPoint&             atm_profile,
                  const ArrayOfPropmatVector&        spectral_propmat_profile,
                  const ArrayOfPropmatMatrix&        spectral_propmat_jac_profile,
                  const ArrayOfStokvecVector&        spectral_nlte_srcvec_profile,
                  const ArrayOfStokvecMatrix&        spectral_nlte_srcvec_jac_profile) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not arr::same_size(alt_grid,
                                        atm_profile,
                                        spectral_propmat_profile,
                                        spectral_propmat_jac_profile,
                                        spectral_nlte_srcvec_profile,
                                        spectral_nlte_srcvec_jac_profile),
                     R"(Profile input must have the same size:

alt_grid size                         {}
atm_profile size                      {}
spectral_propmat_profile size         {}
spectral_propmat_jac_profile size     {}
spectral_nlte_srcvec_profile size     {}
spectral_nlte_srcvec_jac_profile size {}
)",
                     alt_grid.size(),
                     atm_profile.size(),
                     spectral_propmat_profile.size(),
                     spectral_propmat_jac_profile.size(),
                     spectral_nlte_srcvec_profile.size(),
                     spectral_nlte_srcvec_jac_profile.size())

  const std::vector<Size> ips{
      std::from_range, ray_path | stdv::transform([&](const PropagationPathPoint& ray_point) -> Size {
                         const Size ip =
                             stdr::distance(alt_grid.begin(), stdr::lower_bound(alt_grid, ray_point.altitude()));

                         ARTS_USER_ERROR_IF(ip == alt_grid.size(),
                                            "Ray path point altitude {} is out of bounds of the altitude grid [{}, {}]",
                                            ray_point.altitude(),
                                            alt_grid.front(),
                                            alt_grid.back())

                         ARTS_USER_ERROR_IF(alt_grid[ip] != ray_point.altitude(),
                                            R"(Ray path point altitude {} does not match any point in the altitude grid.

Grid points: {:B,}
)",
                                            ray_point.altitude(),
                                            alt_grid)

                         return ip;
                       })};

  const Size n = ips.size();

  spectral_propmat_path.resize(n);
  spectral_propmat_jac_path.resize(n);
  spectral_nlte_srcvec_path.resize(n);
  spectral_nlte_srcvec_jac_path.resize(n);
  atm_path.resize(n);
  freq_grid_path.resize(n);
  freq_wind_shift_jac_path.resize(n);

  for (Size i = 0; i < n; i++) {
    const Size ip                    = ips[i];
    atm_path[i]                      = atm_profile[ip];
    spectral_propmat_path[i]         = spectral_propmat_profile[ip];
    spectral_propmat_jac_path[i]     = spectral_propmat_jac_profile[ip];
    spectral_nlte_srcvec_path[i]     = spectral_nlte_srcvec_profile[ip];
    spectral_nlte_srcvec_jac_path[i] = spectral_nlte_srcvec_jac_profile[ip];
    freq_grid_path[i]                = freq_grid;
    freq_wind_shift_jac_path[i]      = Vector3{0, 0, 0};
  }
}
ARTS_METHOD_ERROR_CATCH

void ProfileFromAltitude(const Workspace&       ws,
                         ArrayOfAtmPoint&       atm_profile,
                         ArrayOfPropmatVector&  spectral_propmat_profile,
                         ArrayOfPropmatMatrix&  spectral_propmat_jac_profile,
                         ArrayOfStokvecVector&  spectral_nlte_srcvec_profile,
                         ArrayOfStokvecMatrix&  spectral_nlte_srcvec_jac_profile,
                         const Agenda&          spectral_propmat_agenda,
                         const JacobianTargets& jac_targets,
                         const AscendingGrid&   alt_grid,
                         const AtmField&        atm_field,
                         const AscendingGrid&   freq_grid,
                         const Numeric&         lat,
                         const Numeric&         lon) {
  const ArrayOfPropagationPathPoint ray_path{std::from_range,
                                             alt_grid | stdv::transform([lat, lon](const Numeric& alt) {
                                               PropagationPathPoint p{.pos_type = PathPositionType::atm,
                                                                      .los_type = PathPositionType::atm,
                                                                      .pos      = {alt, lat, lon},
                                                                      .los      = {0, 0}};
                                               return p;
                                             })};

  const ArrayOfAscendingGrid freq_grid_path(ray_path.size(), freq_grid);
  const ArrayOfVector3       freq_wind_shift_jac_path(ray_path.size(), Vector3{0, 0, 0});
  atm_pathFromPath(atm_profile, ray_path, atm_field);
  spectral_propmat_pathFromPath(ws,
                                spectral_propmat_profile,
                                spectral_nlte_srcvec_profile,
                                spectral_propmat_jac_profile,
                                spectral_nlte_srcvec_jac_profile,
                                spectral_propmat_agenda,
                                freq_grid_path,
                                freq_wind_shift_jac_path,
                                jac_targets,
                                ray_path,
                                atm_profile);
}
