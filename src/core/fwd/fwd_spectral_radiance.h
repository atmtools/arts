#pragma once

#include "fwd_propmat.h"

namespace fwd {
class geometric_spectral_radiance_operator {
  std::vector<Numeric> alt;
  std::vector<propmat_operator> pm;

 public:
 
 Stokvec operator()(Numeric frequency, Numeric zenith, Numeric azimuth=0.0) const;
};
}  // namespace fwd