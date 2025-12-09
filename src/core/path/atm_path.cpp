#include "atm_path.h"

#include <arts_conversions.h>
#include <arts_omp.h>
#include <debug.h>
#include <nonstd.h>

#include <algorithm>
#include <vector>

#include "path_point.h"

ArrayOfAtmPoint &atm_path_resize(ArrayOfAtmPoint &atm_path,
                                 const ArrayOfPropagationPathPoint &ppath) {
  atm_path.resize(ppath.size());
  return atm_path;
}

void forward_atm_path(ArrayOfAtmPoint &atm_path,
                      const ArrayOfPropagationPathPoint &rad_path,
                      const AtmField &atm) {
  stdr::transform(
      rad_path,
      atm_path.begin(),
      [&atm](const auto &p) { return atm.at(p); },
      &PropagationPathPoint::pos);
}

ArrayOfAtmPoint forward_atm_path(const ArrayOfPropagationPathPoint &rad_path,
                                 const AtmField &atm) {
  ArrayOfAtmPoint atm_path(rad_path.size());
  forward_atm_path(atm_path, rad_path, atm);
  return atm_path;
}

void forward_path_freq(AscendingGrid &path_freq,
                       const AscendingGrid &main_freq,
                       const PropagationPathPoint &rad_path,
                       const AtmPoint &atm_path) {
  auto dot_prod = [&]() {
    const auto &[u, v, w] = atm_path.wind;
    const auto &[za, aa]  = rad_path.los;
    const auto f          = std::hypot(u, v, w);
    const auto za_f       = std::acos(w / f);
    const auto aa_f       = std::atan2(u, v);
    const auto za_p       = Conversion::deg2rad(180 - za);
    const auto aa_p       = Conversion::deg2rad(aa + 180);
    return f * ((f == 0) ? 1.0
                         : (std::cos(za_f) * std::cos(za_p) +
                            std::sin(za_f) * std::sin(za_p) *
                                std::cos(aa_f - aa_p)));
  };

  const Numeric fac = 1.0 - dot_prod() / Constant::speed_of_light;

  ARTS_USER_ERROR_IF(
      fac < 0 or nonstd::isnan(fac), "Bad frequency scaling factor: {}", fac)

  Vector tmp = std::move(path_freq).rvec();
  stdr::transform(
      main_freq, tmp.begin(), [fac](const auto &f) { return fac * f; });
  path_freq = std::move(tmp);
}

void forward_path_freq(ArrayOfAscendingGrid &path_freq,
                       const AscendingGrid &main_freq,
                       const ArrayOfPropagationPathPoint &rad_path,
                       const ArrayOfAtmPoint &atm_path) {
  if (arts_omp_in_parallel()) {
    for (Size ip = 0; ip < atm_path.size(); ip++) {
      forward_path_freq(path_freq[ip], main_freq, rad_path[ip], atm_path[ip]);
    }
  } else {
    String error{};
#pragma omp parallel for
    for (Size ip = 0; ip < atm_path.size(); ip++) {
      try {
        forward_path_freq(path_freq[ip], main_freq, rad_path[ip], atm_path[ip]);
      } catch (const std::exception &e) {
#pragma omp critical
        error = std::format("{}\n", std::string_view(e.what()));
      }
    }

    if (not error.empty()) throw std::runtime_error(error);
  }
}

namespace {
void extract1D(ArrayOfAtmPoint &atm_path,
               const AtmField &atm_field,
               const Vector &z_grid,
               const Vector &lat_grid,
               const Vector &lon_grid) {
  const Size n = atm_path.size();

  ARTS_USER_ERROR_IF(n < 1, "Empty path")
  ARTS_USER_ERROR_IF(z_grid.size() != 1 and z_grid.size() != n,
                     "Bad altitude grid")
  ARTS_USER_ERROR_IF(lat_grid.size() != 1 and lat_grid.size() != n,
                     "Bad latitude grid")
  ARTS_USER_ERROR_IF(lon_grid.size() != 1 and lon_grid.size() != n,
                     "Bad longitude grid")

  const auto at = [](const Vector &x, const Size i) {
    return x.size() > 0 ? x[i] : x.front();
  };

  for (Size i = 0; i < n; i++) {
    atm_path[i] = atm_field.at(at(z_grid, i), at(lat_grid, i), at(lon_grid, i));
  }
}
}  // namespace

ArrayOfAtmPoint extract1D(const AtmField &atm_field,
                          const Vector &z_grid,
                          const Vector &lat_grid,
                          const Vector &lon_grid) {
  const Index n = z_grid.size();
  const Index m = lat_grid.size();
  const Index l = lon_grid.size();

  ArrayOfAtmPoint atm_path(std::max({n, m, l}));

  extract1D(atm_path, atm_field, z_grid, lat_grid, lon_grid);

  return atm_path;
}
