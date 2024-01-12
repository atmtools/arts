#include "atm_path.h"

#include <arts_conversions.h>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "path_point.h"

ArrayOfAtmPoint &atm_path_resize(ArrayOfAtmPoint &atm_path,
                                 const Ppath &ppath) {
  atm_path.resize(ppath.np);
  return atm_path;
}

ArrayOfAtmPoint &atm_path_resize(ArrayOfAtmPoint &atm_path,
                                 const ArrayOfPropagationPathPoint &ppath) {
  atm_path.resize(ppath.size());
  return atm_path;
}

void forward_atm_path(ArrayOfAtmPoint &atm_path,
                      const Ppath &ppath,
                      const AtmField &atm) {
  Vector alt{ppath.pos(joker, 0)};
  Vector lat{ppath.pos(joker, 1)};
  Vector lon{ppath.pos(joker, 2)};
  atm_path = atm.at(alt, lat, lon);
}

ArrayOfAtmPoint forward_atm_path(const Ppath &ppath, const AtmField &atm) {
  ArrayOfAtmPoint atm_path(ppath.np);
  forward_atm_path(atm_path, ppath, atm);
  return atm_path;
}

void forward_atm_path(ArrayOfAtmPoint &atm_path,
                      const ArrayOfPropagationPathPoint &rad_path,
                      const AtmField &atm) {
  //FIXME: THIS IS AN UGLY ALLOCATION, atm.at should take a view somehow
  Vector alt(rad_path.size()), lat(rad_path.size()), lon(rad_path.size());
  for (Size ip = 0; ip < rad_path.size(); ip++) {
    alt[ip] = rad_path[ip].pos(0);
    lat[ip] = rad_path[ip].pos(1);
    lon[ip] = rad_path[ip].pos(2);
  }
  atm.at(atm_path, alt, lat, lon);
}

ArrayOfAtmPoint forward_atm_path(const ArrayOfPropagationPathPoint &rad_path,
                                 const AtmField &atm) {
  ArrayOfAtmPoint atm_path(rad_path.size());
  forward_atm_path(atm_path, rad_path, atm);
  return atm_path;
}

ArrayOfVector &path_freq_resize(ArrayOfVector &path_freq,
                                const Vector &main_freq,
                                const ArrayOfAtmPoint &atm_path) {
  path_freq.resize(atm_path.size());
  for (auto &v : path_freq) v.resize(main_freq.nelem());
  return path_freq;
}

void forward_path_freq(ArrayOfVector &path_freq,
                       const Vector &main_freq,
                       const Ppath &ppath,
                       const ArrayOfAtmPoint &atm_path,
                       const Numeric along_path_atm_speed) {
  auto dot_prod = [&](Index ip) {
    const auto &[u, v, w] = atm_path[ip].wind;
    const auto &za = ppath.los(ip, 0);
    const auto &aa = ppath.los(ip, 1);
    const auto f = std::hypot(u, v, w);
    const auto za_f = std::acos(w / f);
    const auto aa_f = std::atan2(u, v);
    const auto za_p = Conversion::deg2rad(180 - za);
    const auto aa_p = Conversion::deg2rad(aa + 180);
    return f * ((f == 0) ? 1.0
                         : (std::cos(za_f) * std::cos(za_p) +
                            std::sin(za_f) * std::sin(za_p) *
                                std::cos(aa_f - aa_p)));
  };

  for (Size ip = 0; ip < atm_path.size(); ip++) {
    std::transform(main_freq.begin(),
                   main_freq.end(),
                   path_freq[ip].begin(),
                   [fac = 1.0 - (along_path_atm_speed + dot_prod(ip)) /
                                    Constant::speed_of_light](const auto &f) {
                     return fac * f;
                   });
  }
}

void forward_path_freq(ArrayOfVector &path_freq,
                       const Vector &main_freq,
                       const ArrayOfPropagationPathPoint &rad_path,
                       const ArrayOfAtmPoint &atm_path,
                       const Numeric along_path_atm_speed) {
  auto dot_prod = [&](Index ip) {
    const auto &[u, v, w] = atm_path[ip].wind;
    const auto &za = rad_path[ip].los[0];
    const auto &aa = rad_path[ip].los[1];
    const auto f = std::hypot(u, v, w);
    const auto za_f = std::acos(w / f);
    const auto aa_f = std::atan2(u, v);
    const auto za_p = Conversion::deg2rad(180 - za);
    const auto aa_p = Conversion::deg2rad(aa + 180);
    return f * ((f == 0) ? 1.0
                         : (std::cos(za_f) * std::cos(za_p) +
                            std::sin(za_f) * std::sin(za_p) *
                                std::cos(aa_f - aa_p)));
  };

  for (Size ip = 0; ip < atm_path.size(); ip++) {
    std::transform(main_freq.begin(),
                   main_freq.end(),
                   path_freq[ip].begin(),
                   [fac = 1.0 - (along_path_atm_speed + dot_prod(ip)) /
                                    Constant::speed_of_light](const auto &f) {
                     return fac * f;
                   });
  }
}

ArrayOfVector forward_path_freq(const Vector &main_freq,
                                const Ppath &ppath,
                                const ArrayOfAtmPoint &atm_path,
                                const Numeric along_path_atm_speed) {
  ArrayOfVector path_freq(atm_path.size(), Vector(main_freq.nelem()));
  forward_path_freq(
      path_freq, main_freq, ppath, atm_path, along_path_atm_speed);
  return path_freq;
}

void extract1D(ArrayOfAtmPoint &atm_path,
               const AtmField &atm_field,
               const Vector &z_grid,
               const Vector &lat_grid,
               const Vector &lon_grid) {
  const Index n = atm_path.size();

  const auto correct_size = [n](const auto &x) {
    const Index m = x->nelem();
    return m == 1 or m == n;
  };
  const std::array<const Vector *, 3> d{&z_grid, &lat_grid, &lon_grid};
  ARTS_USER_ERROR_IF(not std::ranges::all_of(d, correct_size),
                     "All grids must be 1- or n-long"
                     "\nz_grid: ",
                     z_grid,
                     "\nlat_grid: ",
                     lat_grid,
                     "\nlon_grid: ",
                     lon_grid)

  std::array<std::shared_ptr<Vector>, 3> grids;
  for (Index i = 0; i < 3; ++i) {
    if (d[i]->size() == 1) {
      grids[i] =
          std::shared_ptr<Vector>(const_cast<Vector *>(d[i]), [](void *) {});
    } else {
      grids[i] = std::make_shared<Vector>(n, (*d[i])[0]);
    }
  }

  atm_field.at(atm_path, *grids[0].get(), *grids[1].get(), *grids[2].get());
}

ArrayOfAtmPoint extract1D(const AtmField &atm_field,
                          const Vector &z_grid,
                          const Vector &lat_grid,
                          const Vector &lon_grid) {
  const Index n = z_grid.size();
  const Index m = lat_grid.size();
  const Index l = lon_grid.size();

  ArrayOfAtmPoint atm_point(std::max({n, m, l}));
  if (m == l and n == m) {
    atm_field.at(atm_point, z_grid, lat_grid, lon_grid);
  } else if (m == 1 and l == 1) {
    atm_field.at(atm_point,
                 z_grid,
                 Vector(n, lat_grid.front()),
                 Vector(n, lon_grid.front()));
  } else if (m == n and l == 1) {
    atm_field.at(atm_point, z_grid, lat_grid, Vector(n, lon_grid.front()));
  } else if (m == 1 and l == n) {
    atm_field.at(atm_point, z_grid, Vector(n, lat_grid.front()), lon_grid);
  } else {
    ARTS_USER_ERROR(
        "lat_grid and lon_grid must either be of size 1 or of size z_grid.size()\nz_grid: ",
        z_grid,
        "\nlat_grid: ",
        lat_grid,
        "\nlon_grid: ",
        lon_grid)
  }
  return atm_point;
}
