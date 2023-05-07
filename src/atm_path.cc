#include "atm_path.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "auto_md.h"
#include "gridded_fields.h"
#include "matpack_concepts.h"
#include "matpack_view.h"
#include "propagationmatrix.h"
#include "species_tags.h"
#include "transmissionmatrix.h"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

ArrayOfAtmPoint &atm_path_resize(ArrayOfAtmPoint &atm_path,
                                 const Ppath &ppath) {
  atm_path.resize(ppath.np);
  return atm_path;
}

void forward_atm_path(ArrayOfAtmPoint &atm_path, const Ppath &ppath,
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

ArrayOfVector &path_freq_resize(ArrayOfVector &path_freq,
                                const Vector &main_freq,
                                const ArrayOfAtmPoint &atm_path) {
  path_freq.resize(atm_path.nelem());
  for (auto &v : path_freq)
    v.resize(main_freq.nelem());
  return path_freq;
}

void forward_path_freq(ArrayOfVector &path_freq, const Vector &main_freq,
                       const Ppath &ppath, const ArrayOfAtmPoint &atm_path,
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

  for (Index ip = 0; ip < atm_path.nelem(); ip++) {
    std::transform(main_freq.begin(), main_freq.end(), path_freq[ip].begin(),
                   [fac = 1.0 - (along_path_atm_speed + dot_prod(ip)) /
                                    Constant::speed_of_light](const auto &f) {
                     return fac * f;
                   });
  }
}

ArrayOfVector forward_path_freq(const Vector &main_freq, const Ppath &ppath,
                                const ArrayOfAtmPoint &atm_path,
                                const Numeric along_path_atm_speed) {
  ArrayOfVector path_freq(atm_path.nelem(), Vector(main_freq.nelem()));
  forward_path_freq(path_freq, main_freq, ppath, atm_path,
                    along_path_atm_speed);
  return path_freq;
}

/** Returns a unique z_grid and the positions of these values in the original path
 * 
 * @param ppath As WSV
 * @return std::pair<Vector, ArrayOfArrayOfIndex> z_grid, list of positions
 */
std::pair<Vector, ArrayOfArrayOfIndex> unique_z_grid(const Ppath& ppath) {
  const Vector z{ppath.pos(joker, 0)};
  
  Vector z_grid{z};
  std::sort(z_grid.begin(), z_grid.end());
  z_grid = Vector{std::vector<Numeric>{z_grid.begin(), std::unique(z_grid.begin(), z_grid.end())}};

  ArrayOfArrayOfIndex pos(z_grid.nelem());
  if (z.nelem() == z_grid.nelem()) {
    // The grid is completely sorted with unique altitudes

    if (z_grid.front() == z.front()) {
      for (Index i=0; i<pos.nelem(); i++) {
        pos[i] = ArrayOfIndex{i};
      }
    } else {
      for (Index i=0; i<pos.nelem(); i++) {
        pos[i] = ArrayOfIndex{pos.nelem() - i - 1};
      }
    }
  } else {
    // There are multiple values for some altitudes

    for (Index i=0; i<pos.nelem(); i++) {
      for (Index j=0; j<z_grid.nelem(); j++) {
        if (z[j] == z_grid[i]) pos[i].push_back(j);
      }
    }
  }

  return {z_grid, pos};
}

AtmField forward_1d_atm_field(const ArrayOfAtmPoint& atm_path, const Ppath& ppath) {
  ARTS_ASSERT(atm_path.nelem() == ppath.np);
  ARTS_ASSERT(atm_path.nelem() > 0);

  const auto [z, pos] = unique_z_grid(ppath);
  AtmField atm;
  atm.grid[0] = z;
  atm.grid[1] = {ppath.pos(0, 1)};
  atm.grid[2] = {ppath.pos(0, 2)};
  atm.regularized = true;
  atm.top_of_atmosphere = z.back();

  const auto keys = atm_path.front().keys();

  Tensor3 data(atm.regularized_shape());
  for (auto& key: keys) {
    for (Index i=0; i<z.nelem(); i++) {
      data(i, 0, 0) = atm_path[pos[i][0]][key] / static_cast<Numeric>(pos[i].nelem());
      for (Index j=1; j<pos[i].nelem(); j++) {
        data(i, 0, 0) += atm_path[pos[i][j]][key] / static_cast<Numeric>(pos[i].nelem());
      }
    }

    atm[key] = data;
  }

  return atm;
}
