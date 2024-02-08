#include "fwd_path.h"

#include <__algorithm/ranges_binary_search.h>
#include <__algorithm/ranges_lower_bound.h>

#include <algorithm>
#include <string_view>

#include "path_point.h"
#include "sorted_grid.h"

namespace fwd {
constexpr Size xpos(const Numeric x, const AscendingGrid& x_grid) {
  if (x_grid.size() == 1) return 0;

  return std::min<Size>(
      std::max<Index>(
          0,
          std::distance(x_grid.begin(), std::ranges::lower_bound(x_grid, x)) -
              1),
      static_cast<Size>(x_grid.size()) - 2);
}

constexpr Numeric xweight(const Numeric x,
                          const AscendingGrid& x_grid,
                          const Size i0) {
  const Numeric x0 = x_grid[i0];
  if (x0 == x) return 1.0;

  const Numeric x1 = x_grid[i0 + 1];
  if (x1 == x) return 0.0;

  return (x1 - x) / (x1 - x0);
}

void check_grid(const Numeric x [[maybe_unused]],
                const AscendingGrid& x_grid [[maybe_unused]],
                const char* const name [[maybe_unused]]) {
  ARTS_USER_ERROR_IF(x_grid.size() == 0, "No grid of ", name);

  if (x_grid.size() == 1) return;

  ARTS_USER_ERROR_IF(x_grid.front() > x,
                     "Below lower ",
                     name,
                     " point; ",
                     x,
                     " < ",
                     x_grid.front())
  ARTS_USER_ERROR_IF(x_grid.back() < x,
                     "Above upper ",
                     name,
                     " point; ",
                     x,
                     " > ",
                     x_grid.back());
}

std::vector<path> geometric_planar(const Vector3 pos,
                                   const Vector2 los,
                                   const AscendingGrid& alt,
                                   const AscendingGrid& lat,
                                   const AscendingGrid& lon) {
  const Numeric zenith = los[0];
  ARTS_USER_ERROR_IF(zenith == 90.0,
                     "Not a valid zenith angle for geometric planar radiation")

  const Numeric altitude = pos[0];
  const Numeric latitude = pos[1];
  const Numeric longitude = pos[2];
  check_grid(altitude, alt, "altitude");
  check_grid(latitude, lat, "latitude");
  check_grid(longitude, lon, "longitude");

  const Size ialt = xpos(altitude, alt);
  const Size ilat = xpos(latitude, lat);
  const Size ilon = xpos(longitude, lon);
  const auto walt = xweight(altitude, alt, ialt);
  const auto wlat = xweight(latitude, lat, ilat);
  const auto wlon = xweight(longitude, lon, ilon);

  const bool up_looking = zenith > 90.0;
  const Numeric csc = std::abs(1.0 / Conversion::cosd(zenith));

  std::vector<path> path;
  path.reserve(alt.size());
  path.push_back({
      .point =
          {
              .pos_type = ::path::PositionType::atm,
              .los_type = ::path::PositionType::atm,
              .pos = pos,
              .los = los,
          },

      .alt_index = ialt,
      .lat_index = ilat,
      .lon_index = ilon,

      .alt_weight = walt,
      .lat_weight = wlat,
      .lon_weight = wlon,

      .distance = 0.0,
  });

  if (path.back().alt_weight == 0.0 and
      not std::ranges::binary_search(alt, path.back().point.altitude())) {
    path.back().alt_index -= 1;
    path.back().alt_weight = 1.0;
  }

  if (up_looking) {
    while (path.back().point.altitude() != alt.back()) {
      const auto oldpos = path.back();
      auto& newpos = path.emplace_back(oldpos);

      newpos.alt_index += oldpos.alt_weight == 1.0;
      newpos.alt_weight = 1.0;
      newpos.point.altitude() = alt[newpos.alt_index];
      newpos.distance =
          std::abs(newpos.point.altitude() - oldpos.point.altitude()) * csc;
    }

    path.back().point.los_type = ::path::PositionType::space;
  } else {
    while (path.back().point.altitude() != alt.front()) {
      const auto oldpos = path.back();
      auto& newpos = path.emplace_back(oldpos);

      newpos.alt_index -= oldpos.alt_weight == 1.0;
      newpos.alt_weight = 1.0;
      newpos.point.altitude() = alt[newpos.alt_index];
      newpos.distance =
          std::abs(newpos.point.altitude() - oldpos.point.altitude()) * csc;
    }

    path.back().point.los_type = ::path::PositionType::surface;
  }

  while (path.size() > 1 and path[1].distance == 0.0) {
    path.erase(path.begin());
  }

  return path;
}

std::ostream& operator<<(std::ostream& os, const path& pp) {
  return os << pp.point << ' ' << pp.alt_index << ' ' << pp.lat_index << ' '
            << pp.lon_index << ' ' << pp.alt_weight << ' ' << pp.lat_weight
            << ' ' << pp.lon_weight << ' ' << pp.distance;
}

std::ostream& operator<<(std::ostream& os, const std::vector<path>& pps) {
  for (auto& pp : pps) {
    os << pp << '\n';
  }
  return os;
}
}  // namespace fwd
