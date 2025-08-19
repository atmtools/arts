#include "matpack_mdspan_helpers_gridded_data_t.h"

#include <algorithm>

#include "matpack_mdspan_helpers_grid_t.h"

namespace matpack {
GeodeticField3 make_geodetic(const GriddedField3& gf) {
  GeodeticField3 out{
      .data_name  = gf.data_name,
      .data       = {},  // Will be filled in below
      .grid_names = {"altitude", "latitude", "longitude"},
      .grids      = {},  // Will be filled in below
  };

  std::vector<std::pair<Numeric, Size>> alts;
  std::vector<std::pair<Numeric, Size>> lats;
  std::vector<std::pair<Numeric, Size>> lons;
  alts.reserve(gf.shape()[0]);
  lats.reserve(gf.shape()[1]);
  lons.reserve(gf.shape()[2]);

  for (Size i = 0; i < gf.grid<0>().size(); ++i) {
    alts.emplace_back(gf.grid<0>()[i], i);
  }

  for (Size i = 0; i < gf.grid<1>().size(); ++i) {
    lats.emplace_back(gf.grid<1>()[i], i);
  }

  for (Size i = 0; i < gf.grid<2>().size(); ++i) {
    Numeric lon = gf.grid<2>()[i];
    while (lagrange_interp::loncross::cycle(lon) != lon) {
      lon = lagrange_interp::loncross::cycle(lon);
    }
    lons.emplace_back(lon, i);
  }

  stdr::sort(alts, {}, &std::pair<Numeric, Size>::first);
  auto [e0, x] = stdr::unique(alts, {}, &std::pair<Numeric, Size>::first);
  alts.erase(e0, alts.end());

  stdr::sort(lats, {}, &std::pair<Numeric, Size>::first);
  auto [e1, y] = stdr::unique(lats, {}, &std::pair<Numeric, Size>::first);
  lats.erase(e1, lats.end());
  while (lats.size() > 1 and lats.front().first < -90) lats.erase(lats.begin());
  while (lats.size() > 1 and lats.back().first > 90) lats.erase(lats.end() - 1);

  stdr::sort(lons, {}, &std::pair<Numeric, Size>::first);
  auto [e2, z] = stdr::unique(lons, {}, &std::pair<Numeric, Size>::first);
  lons.erase(e2, lons.end());

  const auto first = [](auto v) { return v.first; };
  out.grid<0>()    = AscendingGrid{alts.begin(), alts.end(), first};
  out.grid<1>()    = Vector{AscendingGrid{lats.begin(), lats.end(), first}};
  out.grid<2>()    = Vector{AscendingGrid{lons.begin(), lons.end(), first}};

  out.data.resize(alts.size(), lats.size(), lons.size());

  for (Size i = 0; i < alts.size(); ++i) {
    for (Size j = 0; j < lats.size(); ++j) {
      for (Size k = 0; k < lons.size(); ++k) {
        out.data[i, j, k] =
            gf.data[alts[i].second, lats[j].second, lons[k].second];
      }
    }
  }

  return out;
}

GeodeticField2 make_geodetic(const GriddedField2& gf) {
  GeodeticField2 out{
      .data_name  = gf.data_name,
      .data       = {},  // Will be filled in below
      .grid_names = {"latitude", "longitude"},
      .grids      = {},  // Will be filled in below
  };

  std::vector<std::pair<Numeric, Size>> lats;
  std::vector<std::pair<Numeric, Size>> lons;
  lats.reserve(gf.shape()[0]);
  lons.reserve(gf.shape()[1]);

  for (Size i = 0; i < gf.grid<0>().size(); ++i) {
    lats.emplace_back(gf.grid<0>()[i], i);
  }

  for (Size i = 0; i < gf.grid<1>().size(); ++i) {
    Numeric lon = gf.grid<1>()[i];
    while (lagrange_interp::loncross::cycle(lon) != lon) {
      lon = lagrange_interp::loncross::cycle(lon);
    }
    lons.emplace_back(lon, i);
  }

  stdr::sort(lats, {}, &std::pair<Numeric, Size>::first);
  auto [e1, y] = stdr::unique(lats, {}, &std::pair<Numeric, Size>::first);
  lats.erase(e1, lats.end());
  while (lats.size() > 1 and lats.front().first < -90) lats.erase(lats.begin());
  while (lats.size() > 1 and lats.back().first > 90) lats.erase(lats.end() - 1);

  stdr::sort(lons, {}, &std::pair<Numeric, Size>::first);
  auto [e2, z] = stdr::unique(lons, {}, &std::pair<Numeric, Size>::first);
  lons.erase(e2, lons.end());

  const auto first = [](auto v) { return v.first; };
  out.grid<0>()    = Vector{AscendingGrid{lats.begin(), lats.end(), first}};
  out.grid<1>()    = Vector{AscendingGrid{lons.begin(), lons.end(), first}};

  out.data.resize(lats.size(), lons.size());

  for (Size j = 0; j < lats.size(); ++j) {
    for (Size k = 0; k < lons.size(); ++k) {
      out.data[j, k] = gf.data[lats[j].second, lons[k].second];
    }
  }

  return out;
}
}  // namespace matpack
