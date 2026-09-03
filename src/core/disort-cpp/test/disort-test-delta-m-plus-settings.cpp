#include <disort.h>

#include <stdexcept>

int main() {
  constexpr Index retained_moments = 2;
  constexpr Index full_moments     = 4;

  DisortSettings settings;
  settings.resize(2, retained_moments, 1, AscendingGrid{1.0}, DescendingGrid{1000.0, 0.0});
  settings.legendre_coefficients.resize(1, 1, full_moments);
  settings.legendre_coefficients[0, 0, joker] = Vector{1.0, 0.9, 0.8, 0.72};

  const auto scaling = disort::delta_m_plus(settings.legendre_coefficients[0], retained_moments);
  settings.fractional_scattering[0] = scaling.fraction;
  settings.delta_m_peak_moments[0]  = scaling.moments;

  const auto solver = settings.init();
  if (solver.delta_m_peak_moments().ncols() != retained_moments)
    throw std::runtime_error("Delta-M-plus retained moments were not assigned to NLeg");
  if (solver.all_legendre_coeffs().ncols() != full_moments)
    throw std::runtime_error("Delta-M-plus full moments were not assigned to NLeg_all");
}
