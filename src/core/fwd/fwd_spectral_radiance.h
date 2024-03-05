#pragma once

#include <path_point.h>
#include <physics_funcs.h>

#include <memory>
#include <ostream>

#include "atm.h"
#include "fwd_path.h"
#include "fwd_propmat.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "surf.h"

namespace fwd {
class spectral_radiance {
  AscendingGrid alt;
  AscendingGrid lat;
  AscendingGrid lon;
  matpack::matpack_data<std::shared_ptr<AtmPoint>, 3> atm;
  matpack::matpack_data<propmat, 3> pm;

  matpack::matpack_data<std::function<Stokvec(Numeric, Vector2)>, 2>
      spectral_radiance_surface;
  matpack::matpack_data<std::function<Stokvec(Numeric, Vector2)>, 2>
      spectral_radiance_space;

  Vector2 ellipsoid;

 public:
  struct as_vector {};

  struct weighted_position {
    Numeric w{0.};
    Index i{0}, j{0}, k{0};
  };

  spectral_radiance() = default;
  spectral_radiance(const spectral_radiance&) = default;
  spectral_radiance(spectral_radiance&&) = default;
  spectral_radiance& operator=(const spectral_radiance&) = default;
  spectral_radiance& operator=(spectral_radiance&&) = default;

  spectral_radiance(AscendingGrid alt,
                    AscendingGrid lat,
                    AscendingGrid lon,
                    const AtmField& atm,
                    const SurfaceField& surf,
                    const std::shared_ptr<ArrayOfAbsorptionBand>& lines,
                    const std::shared_ptr<ArrayOfCIARecord>& cia,
                    const std::shared_ptr<ArrayOfXsecRecord>& xsec,
                    const std::shared_ptr<PredefinedModelData>& predef,
                    Numeric ciaextrap = {},
                    Index ciarobust = {});

  Stokvec operator()(const Numeric f,
                     const std::vector<path>& path_points,
                     const Numeric cutoff_transmission = 1e-6) const;

  StokvecVector operator()(const Numeric f,
                           const std::vector<path>& path_points,
                           spectral_radiance::as_vector) const;

  [[nodiscard]] const AscendingGrid& altitude() const { return alt; }
  [[nodiscard]] const AscendingGrid& latitude() const { return lat; }
  [[nodiscard]] const AscendingGrid& longitude() const { return lon; }

  friend std::ostream& operator<<(std::ostream&, const spectral_radiance&);

  [[nodiscard]] std::vector<path> geometric_planar(const Vector3 pos,
                                                   const Vector2 los) const;
  [[nodiscard]] std::vector<path> from_path(
      const ArrayOfPropagationPathPoint& propagation_path) const;
  void from_path(std::vector<path>& out,
                 const ArrayOfPropagationPathPoint& propagation_path) const;

  [[nodiscard]] constexpr std::array<weighted_position, 8> pos_weights(
      const path& pp) const {
    std::array<weighted_position, 8> out;
    auto* ptr = out.begin();
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

  template <Size N>
  [[nodiscard]] constexpr Stokvec B(
      const Numeric f, const std::array<weighted_position, N>& pos) const {
    Numeric out = 0.0;

    for (const auto& p : pos) {
      if (p.w == 0.0) continue;
      out += p.w * planck(f, atm(p.i, p.j, p.k)->temperature);
    }

    return {out, 0.0, 0.0, 0.0};
  }

  template <Size N>
  [[nodiscard]] constexpr Stokvec Iback(
      const Numeric f,
      const std::array<weighted_position, N>& pos,
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

  template <Size N>
  constexpr std::pair<Propmat, Stokvec> PM(
      const Numeric f,
      const std::array<weighted_position, N>& pos,
      const path& pp) const {
    std::pair<Propmat, Stokvec> out{Propmat{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                    Stokvec{0.0, 0.0, 0.0, 0.0}};

    for (const auto& p : pos) {
      if (p.w == 0.0) continue;
      const auto [propmat, stokvec] = pm(p.i, p.j, p.k)(f, pp.point.los);
      out.first += p.w * propmat;
      out.second += p.w * stokvec;
    }

    return out;
  }
};
}  // namespace fwd

using SpectralRadianceOperator = fwd::spectral_radiance;
