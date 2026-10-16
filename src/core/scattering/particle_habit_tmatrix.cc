#include <tmatrix.h>

#include <cmath>
#include <stdexcept>

#include "particle_habit.h"

namespace scattering {
ParticleHabit ParticleHabit::tmatrix(const Vector&        t_grid,
                                     const Vector&        f_grid,
                                     const Vector&        diameters,
                                     const ComplexMatrix& refractive_index,
                                     Numeric              density,
                                     Numeric              aspect_ratio,
                                     int                  shape,
                                     int                  angles,
                                     Numeric              accuracy) {
  if (!::tmatrix::available()) throw std::runtime_error("T-matrix requires ENABLE_TMATRIX=ON");
  auto check_grid = [](const Vector& grid, const char* name) {
    Numeric previous = 0;
    for (Numeric x : grid) {
      if (!std::isfinite(x) || x <= previous)
        throw std::invalid_argument(std::string(name) + " must be positive and strictly increasing");
      previous = x;
    }
    if (grid.empty()) throw std::invalid_argument(std::string(name) + " must not be empty");
  };
  check_grid(t_grid, "t_grid");
  check_grid(f_grid, "f_grid");
  check_grid(diameters, "diameters");
  if (refractive_index.nrows() != t_grid.ncols() || refractive_index.ncols() != f_grid.ncols())
    throw std::invalid_argument("refractive_index must have shape (temperature, frequency)");
  for (auto row : refractive_index)
    for (Complex m : row)
      if (!std::isfinite(m.real()) || !std::isfinite(m.imag()) || m.real() <= 0 || m.imag() < 0)
        throw std::invalid_argument("refractive_index must be finite with real > 0 and imaginary >= 0");
  if (!std::isfinite(density) || density <= 0 || !std::isfinite(aspect_ratio) || aspect_ratio <= 0)
    throw std::invalid_argument("density and aspect_ratio must be finite and positive");
  if (shape != -1 && shape != -2) throw std::invalid_argument("shape must be -1 (spheroid) or -2 (cylinder)");
  if (angles < 2 || !std::isfinite(accuracy) || accuracy <= 0)
    throw std::invalid_argument("angles must be >= 2 and accuracy must be finite and positive");

  auto   temperatures = std::make_shared<Vector>(t_grid);
  auto   frequencies  = std::make_shared<Vector>(f_grid);
  Vector za(angles);
  for (int i = 0; i < angles; ++i) za[i] = 180.0 * i / (angles - 1);
  auto angular_grid = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(za));
  using SSD         = SingleScatteringData<Format::TRO, Representation::Gridded>;
  std::vector<SSD> particles;
  particles.reserve(diameters.size());
  for (Numeric diameter : diameters) {
    PhaseMatrixData<Format::TRO, Representation::Gridded>      phase(temperatures, frequencies, angular_grid);
    ExtinctionMatrixData<Format::TRO, Representation::Gridded> extinction(temperatures, frequencies);
    AbsorptionVectorData<Format::TRO, Representation::Gridded> absorption(temperatures, frequencies);
    for (Size it = 0; it < t_grid.size(); ++it) {
      for (Size jf = 0; jf < f_grid.size(); ++jf) {
        const Complex m = refractive_index[it, jf];
        // NKMAX=-1 gives one quadrature point at the centre of the symmetric
        // narrow interval: one size, with no second PSD integration here.
        const auto optical    = ::tmatrix::random(diameter / 2,
                                                  Constant::speed_of_light / f_grid[jf],
                                                  aspect_ratio,
                                                  m.real(),
                                                  m.imag(),
                                                  angles,
                                                  accuracy,
                                                  1.0,
                                                  shape);
        extinction[it, jf, 0] = optical.extinction;
        absorption[it, jf, 0] = optical.extinction - optical.scattering;
        const Numeric factor  = optical.scattering / (4 * Constant::pi);
        for (int ia = 0; ia < angles; ++ia) {
          const Muelmat& f  = optical.phase[ia];
          phase[it, jf, ia] = factor * CompactPlanarMuelmat{f[0, 0], f[0, 1], f[1, 1], f[2, 2], f[2, 3], f[3, 3]};
        }
      }
    }
    Numeric dmax;
    if (shape == -1) {
      dmax = diameter * std::max(std::cbrt(aspect_ratio), std::pow(aspect_ratio, -2.0 / 3.0));
    } else {
      const Numeric length = diameter * std::cbrt(2.0 / (3.0 * aspect_ratio * aspect_ratio));
      dmax                 = std::hypot(length, aspect_ratio * length);
    }
    ParticleProperties properties{.name             = shape == -1 ? "T-matrix spheroid" : "T-matrix cylinder",
                                  .source           = "ARTS Mishchenko T-matrix solver",
                                  .refractive_index = "User-supplied temperature/frequency grid",
                                  .mass             = density * Constant::pi / 6 * std::pow(diameter, 3),
                                  .d_veq            = diameter,
                                  .d_max            = dmax};
    auto               back    = phase.extract_backscatter_matrix();
    auto               forward = phase.extract_forwardscatter_matrix();
    particles.emplace_back(properties, phase, extinction, absorption, back, forward);
  }
  return ParticleHabit(std::move(particles));
}
}  // namespace scattering
