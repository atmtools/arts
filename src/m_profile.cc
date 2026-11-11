#include <array_algo.h>
#include <arts_omp.h>
#include <lagrange_interp.h>
#include <workspace.h>

#include <algorithm>
#include <ranges>

namespace {
constexpr std::array wind_keys{AtmKey::wind_u, AtmKey::wind_v, AtmKey::wind_w};

//! The wind of atm_field at pos, zero for the components it does not contain
Vector3 wind_at(const AtmField& atm_field, const Vector3& pos) {
  Vector3 wind{0.0, 0.0, 0.0};
  for (Size i = 0; i < wind_keys.size(); i++) {
    if (atm_field.contains(wind_keys[i])) wind[i] = atm_field[wind_keys[i]].at(pos);
  }
  return wind;
}

/*! freq_grid with points added beyond both ends, so that it covers the shift by a wind of speed max_wind_speed

The added points continue the spacing at either end until they reach f0 (1 - v / c) below the first
frequency f0 and fn (1 + v / c) above the last frequency fn.  An added point that would not be
positive is f0 (1 - v / c) instead.
*/
AscendingGrid wind_extended(const AscendingGrid& freq_grid, const Numeric max_wind_speed) {
  const Size    n    = freq_grid.size();
  const Numeric beta = max_wind_speed / Constant::c;

  if (not(beta < 1.0)) arts::user_error("The wind speed of the profile, {} m/s, is not below the speed of light", max_wind_speed);
  if (n == 0 or beta == 0.0) return freq_grid;
  if (freq_grid.front() <= 0.0)
    arts::user_error("The frequencies must be positive to be shifted by the wind, but the first is {} Hz", freq_grid.front());

  const Numeric f_low   = freq_grid.front() * (1.0 - beta);
  const Numeric f_high  = freq_grid.back() * (1.0 + beta);
  const Numeric df_low  = n > 1 ? freq_grid[1] - freq_grid[0] : freq_grid.front() - f_low;
  const Numeric df_high = n > 1 ? freq_grid[n - 1] - freq_grid[n - 2] : f_high - freq_grid.back();
  const Size    n_low   = static_cast<Size>(std::ceil((freq_grid.front() - f_low) / df_low));
  const Size    n_high  = static_cast<Size>(std::ceil((f_high - freq_grid.back()) / df_high));

  Vector out(n_low + n + n_high);
  for (Size i = 0; i < n_low; i++) out[i] = freq_grid.front() - static_cast<Numeric>(n_low - i) * df_low;
  if (n_low > 0 and out[0] <= 0.0) out[0] = f_low;
  out[Range{n_low, n}] = freq_grid.vec();
  for (Size i = 0; i < n_high; i++) out[n_low + n + i] = freq_grid.back() + static_cast<Numeric>(i + 1) * df_high;

  return AscendingGrid{std::move(out)};
}

//! Linear interpolation weights from freq_grid_profile to the frequencies of a path point
auto frequency_lags(const AscendingGrid& freq_grid_profile,
                    const AscendingGrid& freq_grid,
                    const Numeric        f_extrapolfac) {
  const Size n = freq_grid_profile.size();
  if (n < 2)
    arts::user_error(
        "The wind shifts the frequencies of the path, but *freq_grid_profile* has {} frequencies, too few to interpolate",
        n);

  const Numeric f_min = freq_grid_profile[0] - f_extrapolfac * (freq_grid_profile[1] - freq_grid_profile[0]);
  const Numeric f_max = freq_grid_profile[n - 1] + f_extrapolfac * (freq_grid_profile[n - 1] - freq_grid_profile[n - 2]);
  if (freq_grid.front() < f_min or freq_grid.back() > f_max)
    arts::user_error(R"(The wind shifts the frequencies of a path point to [{}, {}] Hz.
This is beyond [{}, {}] Hz, the extrapolation that f_extrapolfac {} allows from *freq_grid_profile*, [{}, {}] Hz.

Compute the profile with at least the wind speed of the path, or allow more extrapolation.)",
                     freq_grid.front(),
                     freq_grid.back(),
                     f_min,
                     f_max,
                     f_extrapolfac,
                     freq_grid_profile.front(),
                     freq_grid_profile.back());

  // The frequencies are inside the limits, as checked above, so the lags need no limit of their own
  return lagrange_interp::make_lags<1, lagrange_interp::grid_identity>(
      freq_grid_profile.vec(), freq_grid.vec(), 0.0, "Frequency");
}
}  // namespace

void Profile2Path(ArrayOfPropmatVector&              spectral_propmat_path,
                  ArrayOfStokvecVector&              spectral_nlte_srcvec_path,
                  ArrayOfPropmatMatrix&              spectral_propmat_jac_path,
                  ArrayOfStokvecMatrix&              spectral_nlte_srcvec_jac_path,
                  ArrayOfAtmPoint&                   atm_path,
                  ArrayOfAscendingGrid&              freq_grid_path,
                  ArrayOfVector3&                    freq_wind_shift_jac_path,
                  const AscendingGrid&               freq_grid,
                  const AscendingGrid&               freq_grid_profile,
                  const AscendingGrid&               alt_grid,
                  const ArrayOfPropagationPathPoint& ray_path,
                  const AtmField&                    atm_field,
                  const JacobianTargets&             jac_targets,
                  const ArrayOfAtmPoint&             atm_profile,
                  const ArrayOfPropmatVector&        spectral_propmat_profile,
                  const ArrayOfPropmatMatrix&        spectral_propmat_jac_profile,
                  const ArrayOfStokvecVector&        spectral_nlte_srcvec_profile,
                  const ArrayOfStokvecMatrix&        spectral_nlte_srcvec_jac_profile,
                  const Numeric&                     f_extrapolfac) try {
  ARTS_TIME_REPORT

  if (not(f_extrapolfac >= 0.0)) arts::user_error("f_extrapolfac must not be negative, it is {}", f_extrapolfac);

  const std::vector<Size> ips{
      std::from_range, ray_path | stdv::transform([&](const PropagationPathPoint& ray_point) -> Size {
                         const Size ip =
                             stdr::distance(alt_grid.begin(), stdr::lower_bound(alt_grid, ray_point.altitude()));

                         if (ip == alt_grid.size())
                           arts::user_error("Ray path point altitude {} is out of bounds of the altitude grid [{}, {}]",
                                            ray_point.altitude(),
                                            alt_grid.front(),
                                            alt_grid.back());

                         if (alt_grid[ip] != ray_point.altitude())
                           arts::user_error(R"(Ray path point altitude {} does not match any point in the altitude grid.

Grid points: {:B,}
)",
                                            ray_point.altitude(),
                                            alt_grid);

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

  String error{};

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n; i++) {
    try {
      const Size ip = ips[i];

      atm_path[i]      = atm_profile[ip];
      atm_path[i].wind = wind_at(atm_field, ray_path[i].pos);

      freq_gridWindShift(freq_grid_path[i] = freq_grid, freq_wind_shift_jac_path[i], atm_path[i], ray_path[i]);

      if (stdr::equal(freq_grid_path[i], freq_grid_profile)) {
        spectral_propmat_path[i]         = spectral_propmat_profile[ip];
        spectral_propmat_jac_path[i]     = spectral_propmat_jac_profile[ip];
        spectral_nlte_srcvec_path[i]     = spectral_nlte_srcvec_profile[ip];
        spectral_nlte_srcvec_jac_path[i] = spectral_nlte_srcvec_jac_profile[ip];
      } else {
        const auto  lags = frequency_lags(freq_grid_profile, freq_grid_path[i], f_extrapolfac);
        const Size  nf   = freq_grid.size();
        const auto& dk   = spectral_propmat_jac_profile[ip];
        const auto& ds   = spectral_nlte_srcvec_jac_profile[ip];

        spectral_propmat_path[i].resize(nf);
        spectral_nlte_srcvec_path[i].resize(nf);
        spectral_propmat_jac_path[i].resize(dk.nrows(), nf);
        spectral_nlte_srcvec_jac_path[i].resize(ds.nrows(), nf);

        lagrange_interp::reinterp(spectral_propmat_path[i], spectral_propmat_profile[ip], lags);
        lagrange_interp::reinterp(spectral_nlte_srcvec_path[i], spectral_nlte_srcvec_profile[ip], lags);
        for (Index it = 0; it < dk.nrows(); it++)
          lagrange_interp::reinterp(spectral_propmat_jac_path[i][it], dk[it], lags);
        for (Index it = 0; it < ds.nrows(); it++)
          lagrange_interp::reinterp(spectral_nlte_srcvec_jac_path[i][it], ds[it], lags);
      }

      // The profile keeps the frequency derivatives of the wind targets, so turn them into wind derivatives here
      spectral_propmat_jacWindFix(spectral_propmat_jac_path[i],
                                  spectral_nlte_srcvec_jac_path[i],
                                  freq_grid_path[i],
                                  jac_targets,
                                  freq_wind_shift_jac_path[i]);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void ProfileFromAltitude(const Workspace&       ws,
                         ArrayOfAtmPoint&       atm_profile,
                         AscendingGrid&         freq_grid_profile,
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
                         const Numeric&         lon) try {
  ARTS_TIME_REPORT

  const ArrayOfPropagationPathPoint ray_path{std::from_range,
                                             alt_grid | stdv::transform([lat, lon](const Numeric& alt) {
                                               PropagationPathPoint p{.pos_type = PathPositionType::atm,
                                                                      .los_type = PathPositionType::atm,
                                                                      .pos      = {alt, lat, lon},
                                                                      .los      = {0, 0}};
                                               return p;
                                             })};

  atm_pathFromPath(atm_profile, ray_path, atm_field);

  // The largest shift is by the strongest wind of the profile, as if it blew along the line of sight
  const Numeric max_wind_speed =
      stdr::fold_left(atm_profile | stdv::transform([](const AtmPoint& p) {
                        return std::hypot(p.wind[0], p.wind[1], p.wind[2]);
                      }),
                      0.0,
                      [](Numeric a, Numeric b) { return std::max(a, b); });
  freq_grid_profile = wind_extended(freq_grid, max_wind_speed);

  // The profile is of a still atmosphere, Profile2Path takes the wind of each path point from atm_field
  for (AtmPoint& atm_point : atm_profile) atm_point.wind = {0.0, 0.0, 0.0};

  // The propagation matrices are computed without wind.  With a wind shift Jacobian of one, spectral_propmat_jacWindFix
  // multiplies the frequency derivatives of the wind targets by the frequency, which is undone below to keep them
  const ArrayOfAscendingGrid freq_grid_path(ray_path.size(), freq_grid_profile);
  const ArrayOfVector3       freq_wind_shift_jac_path(ray_path.size(), Vector3{1.0, 1.0, 1.0});
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

  for (const AtmKey key : wind_keys) {
    const auto ptr = stdr::find_if(jac_targets.atm, Cmp::eq(key), &Jacobian::AtmTarget::type);
    if (ptr == jac_targets.atm.end()) continue;

    if (freq_grid_profile.front() <= 0.0)
      arts::user_error("The frequency derivatives of the wind targets need positive frequencies, but the first is {} Hz",
                       freq_grid_profile.front());

    const Size it = ptr->target_pos;
    for (Size ip = 0; ip < ray_path.size(); ip++) {
      for (Size iv = 0; iv < freq_grid_profile.size(); iv++) {
        spectral_propmat_jac_profile[ip][it, iv]     /= freq_grid_profile[iv];
        spectral_nlte_srcvec_jac_profile[ip][it, iv] /= freq_grid_profile[iv];
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
