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
          spectral_radiance_surface.shape(), [](Numeric f, Vector2) -> Stokvec {
            return planck(f, Constant::cosmic_microwave_background_temperature);
          }) {
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

    ARTS_USER_ERROR_IF(not errors.empty(), errors)
  }
}

Stokvec spectral_radiance::operator()(
    const Numeric f, const std::vector<path>& path_points) const {
  using std::views::drop;

  ARTS_ASSERT(path_points.size() > 0, "No path points")
  ARTS_ASSERT(path_points.front().distance == 0.0, "Bad path point")

  auto pos = pos_weights(path_points.front());

  if (path_points.size() == 1) {
    return Iback(f, pos, path_points.front());
  }

  auto [K, N] = PM(f, pos, path_points.front());
  Stokvec J = inv(K) * N + B(f, pos);
  Muelmat T{1.0};
  Stokvec I{0.0, 0.0, 0.0, 0.0};

  for (auto& pp : path_points | drop(1)) {
    pos = pos_weights(pp);

    if (pp.point.los_type != PathPositionType::atm) {
      return I += T * Iback(f, pos, pp);
    }

    auto [Ki, Ni] = PM(f, pos, pp);
    const Stokvec Ji = inv(Ki) * Ni + B(f, pos);
    const Muelmat Ti = T * exp(avg(Ki, K), pp.distance);

    if (Ti(0, 0) < 1e-6) {
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
  Stokvec J = inv(K) * N + B(f, pos);
  Numeric r = path_points.back().distance;

  for (auto& pp : reverse_view(path_points) | drop(1)) {
    pos = pos_weights(pp);

    auto [Ki, Ni] = PM(f, pos, pp);
    const Stokvec Ji = inv(Ki) * Ni + B(f, pos);
    const Muelmat T = exp(avg(Ki, K), r);

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

std::vector<path> spectral_radiance::from_path(
    const ArrayOfPropagationPathPoint& propagation_path) const {
  return fwd::path_from_propagation_path(propagation_path, alt, lat, lon);
}
}  // namespace fwd