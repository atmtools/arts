#include <arts_constexpr_math.h>
#include <tmatrix.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <vector>

#include "particle_habit.h"

namespace scattering {
namespace {
void check_grid(const Vector& grid, const char* name) {
  if (grid.empty()) throw std::invalid_argument(std::format("{} must not be empty", name));
  Numeric previous = 0;
  for (Numeric x : grid) {
    if (!std::isfinite(x) || x <= previous)
      throw std::invalid_argument(std::format("{} must be positive and strictly increasing", name));
    previous = x;
  }
}

//! An angle grid of gridded ARO data, whose interpolation needs it strictly increasing
void check_angles(const Vector& grid, const char* name) {
  if (grid.empty()) throw std::invalid_argument(std::format("{} must not be empty", name));
  for (Size i = 0; i < grid.size(); ++i)
    if (!std::isfinite(grid[i]) || (i > 0 && grid[i] <= grid[i - 1]))
      throw std::invalid_argument(std::format("{} must be finite and strictly increasing", name));
}

//! The shape of every diameter: the aspect ratio and coefficients of one entry for all, or of one entry each
class Shapes {
  const Vector&        aspect_ratio;
  const ArrayOfVector& coefficients;
  static const Vector  none;

 public:
  ParticleShape shape;
  Index         degree;

  Shapes(const Vector&        aspect_ratio_,
         const ArrayOfVector& coefficients_,
         Size                 diameters,
         ParticleShape        shape_,
         Index                degree_)
      : aspect_ratio(aspect_ratio_), coefficients(coefficients_), shape(shape_), degree(degree_) {
    if (aspect_ratio.size() != 1 and aspect_ratio.size() != diameters)
      throw std::invalid_argument(
          std::format("aspect_ratio has {} elements for {} diameters; it must have 1 or one per diameter",
                      aspect_ratio.size(),
                      diameters));
    if (coefficients.size() > 1 and coefficients.size() != diameters)
      throw std::invalid_argument(
          std::format("coefficients has {} elements for {} diameters; it must have 0, 1 or one per diameter",
                      coefficients.size(),
                      diameters));
    if (shape != ParticleShape::generalized_chebyshev and not coefficients.empty())
      throw std::invalid_argument("coefficients are only for generalized Chebyshev particles");
    if (shape != ParticleShape::chebyshev and degree != 0)
      throw std::invalid_argument(
          std::format("degree is {}; it is only for Chebyshev particles and must be 0 for the others", degree));
  }

  [[nodiscard]] Numeric aspect(Size i) const { return aspect_ratio[aspect_ratio.size() == 1 ? 0 : i]; }

  [[nodiscard]] const Vector& coefficient(Size i) const {
    if (coefficients.empty()) return none;
    return coefficients[coefficients.size() == 1 ? 0 : i];
  }

  //! The shape of diameter i: aspect_ratio is the spheroid's aspect ratio, the cylinder's diameter-to-length ratio
  //! or the Chebyshev deformation, and the coefficients are the generalized Chebyshev particle's
  [[nodiscard]] ::tmatrix::Shape at(Size i) const {
    switch (shape) {
      case ParticleShape::spheroid:              return ::tmatrix::Shape::spheroid(aspect(i));
      case ParticleShape::cylinder:              return ::tmatrix::Shape::cylinder(aspect(i));
      case ParticleShape::chebyshev:             return ::tmatrix::Shape::chebyshev(degree, aspect(i));
      case ParticleShape::generalized_chebyshev: return ::tmatrix::Shape::generalized_chebyshev(coefficient(i));
    }
    throw std::invalid_argument("shape must be a ParticleShape");
  }

  //! The maximum over the volume-equivalent diameter of diameter i; checks the shape
  [[nodiscard]] Numeric dmax_ratio(Size i) const { return ::tmatrix::maximum_dimension(1.0, at(i)) / 2; }

  [[nodiscard]] const char* name() const {
    switch (shape) {
      case ParticleShape::spheroid:              return "T-matrix spheroid";
      case ParticleShape::cylinder:              return "T-matrix cylinder";
      case ParticleShape::chebyshev:             return "T-matrix Chebyshev particle";
      case ParticleShape::generalized_chebyshev: return "T-matrix generalized Chebyshev particle";
    }
    std::unreachable();
  }
};

const Vector Shapes::none{};

//! The checks of both factories, and the maximum dimension ratio of every diameter
std::vector<Numeric> check(const Vector&        t_grid,
                           const Vector&        f_grid,
                           const Vector&        diameters,
                           const ComplexMatrix& refractive_index,
                           Numeric              density,
                           const Shapes&        shapes,
                           Numeric              accuracy) {
  if (!::tmatrix::available()) throw std::runtime_error("T-matrix requires ENABLE_TMATRIX=ON");
  check_grid(t_grid, "t_grid");
  check_grid(f_grid, "f_grid");
  check_grid(diameters, "diameters");
  if (refractive_index.nrows() != t_grid.ncols() || refractive_index.ncols() != f_grid.ncols())
    throw std::invalid_argument(
        std::format("refractive_index has shape ({}, {}); it must be (temperature, frequency), ({}, {})",
                    refractive_index.nrows(),
                    refractive_index.ncols(),
                    t_grid.size(),
                    f_grid.size()));
  for (auto row : refractive_index)
    for (Complex m : row)
      if (!std::isfinite(m.real()) || !std::isfinite(m.imag()) || m.real() <= 0 || m.imag() < 0)
        throw std::invalid_argument("refractive_index must be finite with real > 0 and imaginary >= 0");
  if (!std::isfinite(density) || density <= 0)
    throw std::invalid_argument(std::format("density is {}; it must be finite and positive", density));
  if (!std::isfinite(accuracy) || accuracy <= 0)
    throw std::invalid_argument(std::format("accuracy is {}; it must be finite and positive", accuracy));
  // The shapes are checked as by the solver, before any solver call
  std::vector<Numeric> dmax_ratio(diameters.size());
  for (Size i = 0; i < diameters.size(); ++i) dmax_ratio[i] = shapes.dmax_ratio(i);
  return dmax_ratio;
}

ParticleProperties properties(const Shapes& shapes, Numeric density, Numeric diameter, Numeric dmax_ratio) {
  return {.name             = shapes.name(),
          .source           = "ARTS Mishchenko T-matrix solver",
          .refractive_index = "User-supplied temperature/frequency grid",
          .mass             = density * Constant::pi / 6 * Math::pow3(diameter),
          .d_veq            = diameter,
          .d_max            = dmax_ratio * diameter};
}
}  // namespace

ParticleHabit ParticleHabit::tmatrix(const Vector&        t_grid,
                                     const Vector&        f_grid,
                                     const Vector&        diameters,
                                     const ComplexMatrix& refractive_index,
                                     Numeric              density,
                                     const Vector&        aspect_ratio,
                                     ParticleShape        shape,
                                     Index                degree,
                                     const ArrayOfVector& coefficients,
                                     Numeric              accuracy,
                                     Index                surface_quadrature) {
  const Shapes shapes(aspect_ratio, coefficients, diameters.size(), shape, degree);
  const auto   dmax_ratio = check(t_grid, f_grid, diameters, refractive_index, density, shapes, accuracy);

  auto temperatures = std::make_shared<Vector>(t_grid);
  auto frequencies  = std::make_shared<Vector>(f_grid);
  using SSD         = SingleScatteringData<Format::TRO, Representation::Spectral>;
  std::vector<SSD> particles;
  particles.reserve(diameters.size());
  for (Size i = 0; i < diameters.size(); ++i) {
    const Numeric                                               diameter = diameters[i];
    ExtinctionMatrixData<Format::TRO, Representation::Spectral> extinction(temperatures, frequencies);
    AbsorptionVectorData<Format::TRO, Representation::Spectral> absorption(temperatures, frequencies);
    // The exact series of every temperature and frequency, as differential cross sections; the particle holds
    // them to the degree of the longest, the others zero beyond theirs
    std::vector<CompactPlanarMuelmatVector> series(t_grid.size() * f_grid.size());
    Index                                   series_degree = 0;
    for (Size it = 0; it < t_grid.size(); ++it) {
      for (Size jf = 0; jf < f_grid.size(); ++jf) {
        const Complex m = refractive_index[it, jf];
        // NKMAX=-1 gives one quadrature point at the centre of the symmetric
        // narrow interval: one size, with no second PSD integration here.
        // No gridded phase is needed: the series is the data.
        const auto optical      = ::tmatrix::random(diameter / 2,
                                                    Constant::speed_of_light / f_grid[jf],
                                                    shapes.at(i),
                                                    m.real(),
                                                    m.imag(),
                                                    0,
                                                    accuracy,
                                                    EquivalentRadius::volume,
                                                    TMatrixSizeDistribution::gamma,
                                                    0.1,
                                                    1.0,
                                                    -1,
                                                    surface_quadrature);
        extinction[it, jf].A()  = optical.extinction;
        absorption[it, jf].I()  = optical.extinction - optical.scattering;
        auto& s                 = series[it * f_grid.size() + jf];
        s                       = optical.legendre;
        s                      *= optical.scattering / (4 * Constant::pi);
        series_degree           = std::max<Index>(series_degree, s.size() - 1);
      }
    }
    PhaseMatrixData<Format::TRO, Representation::Spectral> phase(temperatures, frequencies, series_degree);
    for (Size it = 0; it < t_grid.size(); ++it)
      for (Size jf = 0; jf < f_grid.size(); ++jf) {
        const auto& s = series[it * f_grid.size() + jf];
        for (Size l = 0; l < s.size(); ++l) phase[it, jf, l] = s[l];
      }
    auto back    = phase.extract_backscatter_matrix();
    auto forward = phase.extract_forwardscatter_matrix();
    particles.emplace_back(
        properties(shapes, density, diameter, dmax_ratio[i]), phase, extinction, absorption, back, forward);
  }
  return ParticleHabit(std::move(particles));
}

ParticleHabit ParticleHabit::tmatrix_aro(const Vector&          t_grid,
                                         const Vector&          f_grid,
                                         const Vector&          diameters,
                                         const ComplexMatrix&   refractive_index,
                                         Numeric                density,
                                         const Vector&          za_inc_grid,
                                         const Vector&          delta_aa_grid,
                                         const ZenithAngleGrid& za_scat_grid,
                                         const Vector&          tilt_angles,
                                         const Vector&          tilt_weights,
                                         const Vector&          aspect_ratio,
                                         ParticleShape          shape,
                                         Index                  degree,
                                         const ArrayOfVector&   coefficients,
                                         Numeric                accuracy,
                                         Index                  surface_quadrature) {
  const Shapes shapes(aspect_ratio, coefficients, diameters.size(), shape, degree);
  const auto   dmax_ratio = check(t_grid, f_grid, diameters, refractive_index, density, shapes, accuracy);

  auto         temperatures = std::make_shared<Vector>(t_grid);
  auto         frequencies  = std::make_shared<Vector>(f_grid);
  auto         za_inc       = std::make_shared<Vector>(za_inc_grid);
  auto         delta_aa     = std::make_shared<Vector>(delta_aa_grid);
  auto         za_scat      = std::make_shared<ZenithAngleGrid>(za_scat_grid);
  const Vector za_scat_angles{grid_vector(za_scat_grid)};
  check_angles(za_inc_grid, "za_inc_grid");
  check_angles(delta_aa_grid, "delta_aa_grid");
  check_angles(za_scat_angles, "za_scat_grid");
  using SSD = SingleScatteringData<Format::ARO, Representation::Gridded>;
  std::vector<SSD> particles;
  particles.reserve(diameters.size());
  for (Size i = 0; i < diameters.size(); ++i) {
    const Numeric                                         diameter = diameters[i];
    PhaseMatrixData<Format::ARO, Representation::Gridded> phase(temperatures, frequencies, za_inc, delta_aa, za_scat);
    ExtinctionMatrixData<Format::ARO, Representation::Gridded> extinction(temperatures, frequencies, za_inc);
    AbsorptionVectorData<Format::ARO, Representation::Gridded> absorption(temperatures, frequencies, za_inc);
    for (Size it = 0; it < t_grid.size(); ++it) {
      for (Size jf = 0; jf < f_grid.size(); ++jf) {
        const Complex m = refractive_index[it, jf];
        // The tilts and the ranges of the angles are checked by the first call
        const auto optical = ::tmatrix::azimuthally_random(diameter / 2,
                                                           Constant::speed_of_light / f_grid[jf],
                                                           shapes.at(i),
                                                           m.real(),
                                                           m.imag(),
                                                           za_inc_grid,
                                                           delta_aa_grid,
                                                           za_scat_angles,
                                                           tilt_angles,
                                                           tilt_weights,
                                                           accuracy,
                                                           EquivalentRadius::volume,
                                                           surface_quadrature);
        for (Size k = 0; k < za_inc_grid.size(); ++k) {
          extinction[it, jf, k] = optical.extinction[k];
          absorption[it, jf, k] = optical.absorption[k];
          for (Size a = 0; a < delta_aa_grid.size(); ++a)
            for (Size s = 0; s < za_scat_angles.size(); ++s) phase[it, jf, k, a, s] = optical.phase[k, a, s];
        }
      }
    }
    auto back    = phase.extract_backscatter_matrix();
    auto forward = phase.extract_forwardscatter_matrix();
    particles.emplace_back(
        properties(shapes, density, diameter, dmax_ratio[i]), phase, extinction, absorption, back, forward);
  }
  return ParticleHabit(std::move(particles));
}
}  // namespace scattering
