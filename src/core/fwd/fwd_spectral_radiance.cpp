#include "fwd_spectral_radiance.h"

namespace fwd {
  Stokvec spectral_radiance::geometric(
    const Numeric frequency, const Vector3 pos, const Vector2 los) const {
      return Stokvec{};
    }
    
Stokvec spectral_radiance::operator()(
    const Numeric frequency, const Vector3 pos, const Vector2 los) const {
      switch(pathing) {
        case pathing::geometric:
          return geometric(frequency, pos, los);
        default:
          throw std::runtime_error("Invalid pathing type");
      }
    }
}  // namespace fwd