#include <array_algo.h>
#include <arts_omp.h>
#include <lagrange_interp.h>
#include <workspace.h>

#include <algorithm>
#include <limits>
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

void Profile2Path(const Workspace&                   ws,
                  ArrayOfSpectralPropmat&            spectral_propmat_path,
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
                  const ArrayOfSpectralPropmat&      spectral_propmat_profile,
                  const Agenda&                      spectral_propmat_agenda,
                  const Numeric&                     f_extrapolfac) try {
  ARTS_TIME_REPORT

  if (not(f_extrapolfac >= 0.0)) arts::user_error("f_extrapolfac must not be negative, it is {}", f_extrapolfac);

  // A path of fewer than two points has no layer to propagate through, so nothing uses its propagation matrices.
  // The profile does not cover points outside alt_grid.  Neither has propagation.  Points inside alt_grid but between
  // its altitudes, such as the tangent point of a limb path, are not in the profile either, so they are computed.
  constexpr Size outside  = std::numeric_limits<Size>::max();
  constexpr Size between  = outside - 1;
  const bool     no_layer = ray_path.size() < 2;

  const std::vector<Size> ips{std::from_range,
                              ray_path | stdv::transform([&](const PropagationPathPoint& ray_point) -> Size {
                                const Numeric alt = ray_point.altitude();
                                if (no_layer or alt_grid.empty() or alt < alt_grid.front() or alt > alt_grid.back())
                                  return outside;

                                const Size ip = stdr::distance(alt_grid.begin(), stdr::lower_bound(alt_grid, alt));
                                return alt_grid[ip] == alt ? ip : between;
                              })};

  const Size n  = ips.size();
  const Size nf = freq_grid.size();
  const Size nq = jac_targets.target_count();

  spectral_propmat_path.resize(n);
  atm_path.resize(n);
  freq_grid_path.resize(n);
  freq_wind_shift_jac_path.resize(n);

  // Points that are not in the profile have the atmosphere and frequencies of atm_pathFromPath and
  // freq_grid_pathFromPath.  Those between the altitudes of alt_grid have the propagation of spectral_propmat_pathFromPath,
  // as for the default path agenda, and those without an altitude of the profile have none.
  for (const Size kind : {outside, between}) {
    const std::vector<Size> points{
        std::from_range, stdv::iota(Size{0}, n) | stdv::filter([&ips, kind](Size i) { return ips[i] == kind; })};
    if (points.empty()) continue;

    const ArrayOfPropagationPathPoint path{std::from_range,
                                           points | stdv::transform([&ray_path](Size i) { return ray_path[i]; })};
    ArrayOfAtmPoint                   part_atm_path{};
    ArrayOfAscendingGrid              part_freq_grid_path{};
    ArrayOfVector3                    part_freq_wind_shift_jac_path{};
    ArrayOfSpectralPropmat            part_spectral_propmat_path{};

    atm_pathFromPath(part_atm_path, path, atm_field);
    freq_grid_pathFromPath(part_freq_grid_path, part_freq_wind_shift_jac_path, freq_grid, path, part_atm_path);
    if (kind == between) {
      spectral_propmat_pathFromPath(ws,
                                    part_spectral_propmat_path,
                                    spectral_propmat_agenda,
                                    part_freq_grid_path,
                                    part_freq_wind_shift_jac_path,
                                    jac_targets,
                                    path,
                                    part_atm_path);
    } else {
      part_spectral_propmat_path.resize(points.size());
      for (auto& part : part_spectral_propmat_path) part.init(nf, nq);
    }

    for (Size k = 0; k < points.size(); k++) {
      const Size i                     = points[k];
      atm_path[i]                      = std::move(part_atm_path[k]);
      freq_grid_path[i]                = std::move(part_freq_grid_path[k]);
      freq_wind_shift_jac_path[i]      = part_freq_wind_shift_jac_path[k];
      spectral_propmat_path[i]         = std::move(part_spectral_propmat_path[k]);
    }
  }

  String error{};

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n; i++) {
    try {
      const Size ip = ips[i];
      if (ip == outside or ip == between) continue;

      atm_path[i]      = atm_profile[ip];
      atm_path[i].wind = wind_at(atm_field, ray_path[i].pos);

      freq_gridWindShift(freq_grid_path[i] = freq_grid, freq_wind_shift_jac_path[i], atm_path[i], ray_path[i]);

      const SpectralPropmat& profile = spectral_propmat_profile[ip];
      SpectralPropmat&       point   = spectral_propmat_path[i];

      if (stdr::equal(freq_grid_path[i], freq_grid_profile)) {
        point = profile;
      } else {
        const auto lags = frequency_lags(freq_grid_profile, freq_grid_path[i], f_extrapolfac);

        point.K.resize(nf);
        point.nlte.resize(nf);
        point.dK.resize(profile.dK.nrows(), nf);
        point.dnlte.resize(profile.dnlte.nrows(), nf);
        point.nlte_on = profile.nlte_on;

        lagrange_interp::reinterp(point.K, profile.K, lags);
        for (Index it = 0; it < profile.dK.nrows(); it++) lagrange_interp::reinterp(point.dK[it], profile.dK[it], lags);

        // Without a non-LTE source, there are only zeros to interpolate
        if (profile.nlte_on) {
          lagrange_interp::reinterp(point.nlte, profile.nlte, lags);
          for (Index it = 0; it < profile.dnlte.nrows(); it++)
            lagrange_interp::reinterp(point.dnlte[it], profile.dnlte[it], lags);
        } else {
          point.nlte  = Stokvec{0.0};
          point.dnlte = Stokvec{0.0};
        }
      }

      // The profile keeps the frequency derivatives of the wind targets, so turn them into wind derivatives here
      spectral_propmatWindFix(point, freq_grid_path[i], jac_targets, freq_wind_shift_jac_path[i]);
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
                         ArrayOfSpectralPropmat& spectral_propmat_profile,
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

  // The propagation matrices are computed without wind.  With a wind shift Jacobian of one, spectral_propmatWindFix
  // multiplies the frequency derivatives of the wind targets by the frequency, which is undone below to keep them
  const ArrayOfAscendingGrid freq_grid_path(ray_path.size(), freq_grid_profile);
  const ArrayOfVector3       freq_wind_shift_jac_path(ray_path.size(), Vector3{1.0, 1.0, 1.0});
  spectral_propmat_pathFromPath(ws,
                                spectral_propmat_profile,
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
        spectral_propmat_profile[ip].dK[it, iv]    /= freq_grid_profile[iv];
        spectral_propmat_profile[ip].dnlte[it, iv] /= freq_grid_profile[iv];
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
