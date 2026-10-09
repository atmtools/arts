#include "polradtran.h"

#include <debug.h>
#include <integration.h>

namespace polradtran {
quadrature get_quadrature(Index nmu, quadrature_type type) {
  ARTS_USER_ERROR_IF(nmu < 1, "RT3 and RT4 need at least one quadrature node per hemisphere, got nmu = {}", nmu);

  // The positive half of ARTS's 2 nmu-point rule on [-1, 1]
  const auto positive_half = [nmu](const auto& rule) {
    return quadrature{.mu      = Vector{rule.get_nodes()[Range{nmu, nmu}]},
                      .weights = Vector{rule.get_weights()[Range{nmu, nmu}]}};
  };
  switch (type) {
    case quadrature_type::double_gauss: return positive_half(scattering::DoubleGaussQuadrature(2 * nmu));
    case quadrature_type::gauss:        return positive_half(scattering::GaussLegendreQuadrature(2 * nmu));
    case quadrature_type::lobatto:      return positive_half(scattering::LobattoQuadrature(2 * nmu));
  }
  ARTS_USER_ERROR("Unknown quadrature type {}", static_cast<int>(type));
}
}  // namespace polradtran
