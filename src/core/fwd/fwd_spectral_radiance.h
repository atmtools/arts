#pragma once

#include <path_point.h>

#include <limits>
#include <memory>
#include <ostream>

#include "atm.h"
#include "fwd_propmat.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "surf.h"

namespace fwd {
class spectral_radiance {
  AscendingGrid alt;
  std::vector<std::shared_ptr<AtmPoint>> atm;
  std::vector<propmat> pm;
  SurfacePoint surf;

 public:
  struct PosDistance {
    Size i{std::numeric_limits<Size>::max()};
    Numeric r{0.};
  };

  using PathVector = matpack::matpack_data<PosDistance, 1>;
  using PathMatrix = matpack::matpack_data<PosDistance, 2>;
  using PathVectorView = matpack::matpack_view<PosDistance, 1, true, false>;

  spectral_radiance() = default;
  spectral_radiance(const spectral_radiance&) = default;
  spectral_radiance(spectral_radiance&&) = default;
  spectral_radiance& operator=(const spectral_radiance&) = default;
  spectral_radiance& operator=(spectral_radiance&&) = default;

  spectral_radiance(AscendingGrid alt,
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

  [[nodiscard]] const Vector& altitude() const { return alt; }

  friend std::ostream& operator<<(std::ostream&, const spectral_radiance&);

  [[nodiscard]] PathVector geometric_planar(const Numeric altitude, const Numeric zenith) const;
};
}  // namespace fwd

using SpectralRadianceOperator = fwd::spectral_radiance;
