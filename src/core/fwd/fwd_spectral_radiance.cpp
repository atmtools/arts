#include "fwd_spectral_radiance.h"

#include <path_point.h>

#include <algorithm>
#include <memory>
#include <ostream>
#include <ranges>

#include "arts_constants.h"
#include "arts_omp.h"
#include "configtypes.h"
#include "debug.h"
#include "fwd_path.h"
#include "rtepack.h"

namespace fwd {
Stokvec spectral_radiance::B(
    const Numeric f,
    const std::array<spectral_radiance::weighted_position, 8>& pos) const {
  Numeric out = 0.0;

  for (const auto& p : pos) {
    if (p.w == 0.0) continue;
    out += p.w * planck(f, atm(p.i, p.j, p.k)->temperature);
  }

  return {out, 0.0, 0.0, 0.0};
}

Stokvec spectral_radiance::Iback(
    const Numeric f,
    const std::array<spectral_radiance::weighted_position, 8>& pos,
    const path& pp) const {
  Stokvec out{0.0, 0.0, 0.0, 0.0};

  if (pp.point.los_type == PathPositionType::space) {
    for (const auto& p : pos) {
      if (p.w == 0.0) continue;
      out += p.w * spectral_radiance_space(p.j, p.k)(f, pp.point.los);
    }
  } else if (pp.point.los_type == PathPositionType::surface) {
    for (const auto& p : pos) {
      if (p.w == 0.0) continue;
      out += p.w * spectral_radiance_surface(p.j, p.k)(f, pp.point.los);
    }
  }

  return out;
}

std::pair<Propmat, Stokvec> spectral_radiance::PM(
    const Numeric f,
    const std::array<spectral_radiance::weighted_position, 8>& pos,
    const path& pp) const {
  std::pair<Propmat, Stokvec> out{Propmat{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                  Stokvec{0.0, 0.0, 0.0, 0.0}};

  for (const auto& p : pos) {
    if (p.w == 0.0) continue;
    const auto [propmat, stokvec]  = pm(p.i, p.j, p.k)(f, pp.point.los);
    out.first                     += p.w * propmat;
    out.second                    += p.w * stokvec;
  }

  return out;
}

std::array<spectral_radiance::weighted_position, 8>
spectral_radiance::pos_weights(const path& pp) const {
  std::array<weighted_position, 8> out;
  auto ptr = out.begin();
  for (Size i : {0, 1}) {
    for (Size j : {0, 1}) {
      for (Size k : {0, 1}) {
        ptr->i = pp.alt_index + i;
        ptr->j = pp.lat_index + j;
        ptr->k = pp.lon_index + k;
        ptr->w = (i == 0 ? pp.alt_weight : 1.0 - pp.alt_weight) *
                 (j == 0 ? pp.lat_weight : 1.0 - pp.lat_weight) *
                 (k == 0 ? pp.lon_weight : 1.0 - pp.lon_weight);
        ptr++;
      }
    }
  }
  return out;
}

spectral_radiance::spectral_radiance(
    AscendingGrid alt_,
    AscendingGrid lat_,
    AscendingGrid lon_,
    const AtmField& atm_,
    const SurfaceField& surf,
    const std::shared_ptr<ArrayOfAbsorptionBand>& lines,
    const std::shared_ptr<ArrayOfCIARecord>& cia,
    const std::shared_ptr<ArrayOfXsecRecord>& xsec,
    const std::shared_ptr<PredefinedModelData>& predef,
    Numeric ciaextrap,
    Index ciarobust)
    : alt(std::move(alt_)),
      lat(std::move(lat_)),
      lon(std::move(lon_)),
      atm(alt.size(), lat.size(), lon.size()),
      pm(atm.shape()),
      spectral_radiance_surface(lat.size(), lon.size()),
      spectral_radiance_space(
          spectral_radiance_surface.shape(),
          [](Numeric f, Vector2) -> Stokvec {
            return planck(f, Constant::cosmic_microwave_background_temperature);
          }),
      ellipsoid(surf.ellipsoid) {
  ARTS_USER_ERROR_IF(alt.size() == 0, "Must have a sized atmosphere")

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Index j = 0; j < lat.size(); j++) {
      for (Index k = 0; k < lon.size(); k++) {
        spectral_radiance_surface(j, j) = [surf = surf.at(lat[j], lon[k])](
                                              Numeric f, Vector2) -> Stokvec {
          return planck(f, surf.temperature);
        };
      }
    }

    for (Index i = 0; i < alt.size(); i++) {
      for (Index j = 0; j < lat.size(); j++) {
        for (Index k = 0; k < lon.size(); k++) {
          atm(i, j, k) =
              std::make_shared<AtmPoint>(atm_.at(alt[i], lat[j], lon[k]));
          pm(i, j, k) = propmat(
              atm(i, j, k), lines, cia, xsec, predef, ciaextrap, ciarobust);
        }
      }
    }
  } else {
    String errors{};

#pragma omp parallel for collapse(2)
    for (Index j = 0; j < lat.size(); j++) {
      for (Index k = 0; k < lon.size(); k++) {
        try {
          spectral_radiance_surface(j, j) = [surf = surf.at(lat[j], lon[k])](
                                                Numeric f, Vector2) -> Stokvec {
            return planck(f, surf.temperature);
          };
        } catch (const std::exception& e) {
#pragma omp critical
          errors += e.what();
        }
      }
    }

#pragma omp parallel for collapse(3)
    for (Index i = 0; i < alt.size(); i++) {
      for (Index j = 0; j < lat.size(); j++) {
        for (Index k = 0; k < lon.size(); k++) {
          try {
            atm(i, j, k) =
                std::make_shared<AtmPoint>(atm_.at(alt[i], lat[j], lon[k]));
            pm(i, j, k) = propmat(
                atm(i, j, k), lines, cia, xsec, predef, ciaextrap, ciarobust);
          } catch (const std::exception& e) {
#pragma omp critical
            errors += e.what();
          }
        }
      }
    }

    ARTS_USER_ERROR_IF(not errors.empty(), "{}", errors)
  }
}

Stokvec spectral_radiance::operator()(const Numeric f,
                                      const std::vector<path>& path_points,
                                      const Numeric cutoff_transmission) const {
  using std::views::drop;

  ARTS_ASSERT(path_points.size() > 0, "No path points")
  ARTS_ASSERT(path_points.front().distance == 0.0, "Bad path point")

  auto pos = pos_weights(path_points.front());

  if (path_points.size() == 1) {
    return Iback(f, pos, path_points.front());
  }

  auto [K, N] = PM(f, pos, path_points.front());
  Stokvec J   = inv(K) * N + B(f, pos);
  Muelmat T{1.0};
  Stokvec I{0.0, 0.0, 0.0, 0.0};

  for (auto& pp : path_points | drop(1)) {
    pos = pos_weights(pp);

    if (pp.point.los_type != PathPositionType::atm) {
      return I += T * Iback(f, pos, pp);
    }

    auto [Ki, Ni]    = PM(f, pos, pp);
    const Stokvec Ji = inv(Ki) * Ni + B(f, pos);
    const Muelmat Ti = T * exp(avg(Ki, K), pp.distance);

    if (Ti(0, 0) < cutoff_transmission) {
      return I += Ti * avg(Ji, J);
    }

    I += (T - Ti) * avg(Ji, J);

    J = Ji;
    K = Ki;
    T = Ti;
  }

  return I;
}

StokvecVector spectral_radiance::operator()(
    const Numeric f,
    const std::vector<path>& path_points,
    spectral_radiance::as_vector) const {
  using std::ranges::reverse_view;
  using std::views::drop;

  ARTS_ASSERT(path_points.size() > 0, "No path points")
  ARTS_ASSERT(path_points.front().distance == 0.0, "Bad path point")

  std::vector<Stokvec> out;
  out.reserve(path_points.size());

  auto pos = pos_weights(path_points.back());
  out.emplace_back(Iback(f, pos, path_points.back()));

  auto [K, N] = PM(f, pos, path_points.back());
  Stokvec J   = inv(K) * N + B(f, pos);
  Numeric r   = path_points.back().distance;

  for (auto& pp : reverse_view(path_points) | drop(1)) {
    pos = pos_weights(pp);

    auto [Ki, Ni]    = PM(f, pos, pp);
    const Stokvec Ji = inv(Ki) * Ni + B(f, pos);
    const Muelmat T  = exp(avg(Ki, K), r);

    out.emplace_back(T * (out.back() - avg(J, Ji)) + avg(J, Ji));

    J = Ji;
    K = Ki;
    r = pp.distance;
  }

  std::ranges::reverse(out);
  return out;
}

std::ostream& operator<<(std::ostream& os, const spectral_radiance& sr) {
  return os << "Spectral radiance operator:\n"
            << "  Altitude grid: " << sr.alt << "\n";
}

std::vector<path> spectral_radiance::geometric_planar(const Vector3 pos,
                                                      const Vector2 los) const {
  return fwd::geometric_planar(pos, los, alt, lat, lon);
}

void spectral_radiance::from_path(
    std::vector<path>& out,
    const ArrayOfPropagationPathPoint& propagation_path) const {
  return fwd::path_from_propagation_path(
      out, propagation_path, alt, lat, lon, ellipsoid);
}

std::vector<path> spectral_radiance::from_path(
    const ArrayOfPropagationPathPoint& propagation_path) const {
  return fwd::path_from_propagation_path(
      propagation_path, alt, lat, lon, ellipsoid);
}
}  // namespace fwd