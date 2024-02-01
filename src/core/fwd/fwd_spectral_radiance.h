#pragma once

#include "fwd_propmat.h"

namespace fwd {
ENUMCLASS(pathing, char, geometric)
class spectral_radiance {
  std::vector<Numeric> alt;
  std::vector<propmat> pm;
  Vector2 ellipsoid;
  pathing pathing{pathing::geometric};

 public:
 [[nodiscard]] Stokvec geometric(const Numeric frequency, const Vector3 pos, const Vector2 los) const;
 Stokvec operator()(const Numeric frequency, const Vector3 pos, const Vector2 los) const;
};
}  // namespace fwd