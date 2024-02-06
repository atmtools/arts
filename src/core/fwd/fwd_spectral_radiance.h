#pragma once

#include <path_point.h>

#include <limits>
#include <memory>
#include <ostream>

#include "atm.h"
#include "fwd_propmat.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "surf.h"

namespace fwd {
class spectral_radiance_1d {
  AscendingGrid alt;
  Numeric latitude;
  Numeric longitude;
  std::vector<std::shared_ptr<AtmPoint>> atm;
  std::vector<propmat> pm;

  std::function<Stokvec(Numeric, Vector2)> spectral_radiance_surface;
  std::function<Stokvec(Numeric)> spectral_radiance_space;

 public:
  struct PosDistance {
    Size i{std::numeric_limits<Size>::max()};
    Numeric r{0.};
  };

  using PathVector = matpack::matpack_data<PosDistance, 1>;
  using PathMatrix = matpack::matpack_data<PosDistance, 2>;
  using PathTensor3 = matpack::matpack_data<PosDistance, 3>;
  using PathVectorView = matpack::matpack_view<PosDistance, 1, true, false>;

  spectral_radiance_1d() = default;
  spectral_radiance_1d(const spectral_radiance_1d&) = default;
  spectral_radiance_1d(spectral_radiance_1d&&) = default;
  spectral_radiance_1d& operator=(const spectral_radiance_1d&) = default;
  spectral_radiance_1d& operator=(spectral_radiance_1d&&) = default;

  spectral_radiance_1d(AscendingGrid alt,
                       Numeric lat,
                       Numeric lon,
                       const AtmField& atm,
                       const SurfaceField& surf,
                       std::shared_ptr<AbsorptionBands> lines,
                       std::shared_ptr<ArrayOfCIARecord> cia,
                       std::shared_ptr<ArrayOfXsecRecord> xsec,
                       std::shared_ptr<PredefinedModelData> predef,
                       Numeric ciaextrap = {},
                       Index ciarobust = {});

  Stokvec operator()(const Numeric frequency,
                     const Vector2 los,
                     const PathVectorView) const;

  void operator()(StokvecVectorView out,
                  const Numeric frequency,
                  const Vector2 los,
                  const PathVectorView) const;

  [[nodiscard]] const Vector& altitude() const { return alt; }

  friend std::ostream& operator<<(std::ostream&, const spectral_radiance_1d&);

  [[nodiscard]] PathVector geometric_planar(const Numeric altitude,
                                            const Numeric zenith) const;

  [[nodiscard]] ArrayOfAtmPoint get_atm(const PathVectorView) const;
  [[nodiscard]] ExhaustiveConstVectorView altitude_grid() const;
  [[nodiscard]] ExhaustiveConstVectorView latitude_grid() const;
  [[nodiscard]] ExhaustiveConstVectorView longitude_grid() const;
};
}  // namespace fwd

using SpectralRadianceOperator = fwd::spectral_radiance_1d;
