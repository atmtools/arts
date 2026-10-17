#include "tmatrix.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>

#ifdef ARTS_HAS_TMATRIX
#include "tmatrix_solver.h"
#endif

namespace tmatrix {
namespace {
void require(bool valid, const char* message) {
  if (!valid) throw std::invalid_argument(message);
}
void positive(double x, const char* name) { require(std::isfinite(x) && x > 0, name); }
void check_shape(double aspect, ParticleShape shape, Index degree, const Vector& coefficients) {
  if (!available()) throw std::runtime_error("T-matrix requires ENABLE_TMATRIX=ON");
  require(good_enum(shape), "shape must be a ParticleShape");
  require(shape == ParticleShape::generalized_chebyshev or coefficients.empty(),
          "coefficients are only for generalized Chebyshev particles");
  if (shape == ParticleShape::generalized_chebyshev) {
    require(degree == 0, "degree must be 0 for generalized Chebyshev particles, whose coefficients give the shape");
    require(not coefficients.empty(), "a generalized Chebyshev particle needs its coefficients c_0, ..., c_N");
    // 1 + c_0 > sum_{n >= 1} |c_n| keeps the radius positive
    Numeric margin = 1.0 + coefficients[0];
    for (Size n = 1; n < coefficients.size(); ++n) margin -= std::abs(coefficients[n]);
    require(std::isfinite(margin) && margin > 0,
            "the generalized Chebyshev coefficients must be finite with 1 + c_0 > sum_{n >= 1} |c_n|");
  } else if (shape == ParticleShape::chebyshev) {
    require(degree >= 1, "a Chebyshev particle needs a degree of at least 1");
    require(std::isfinite(aspect) && std::abs(aspect) < 1,
            "aspect_ratio, the Chebyshev deformation, must be finite with magnitude below 1");
  } else {
    require(degree == 0, "degree must be 0 for spheroids and cylinders");
    positive(aspect, "aspect_ratio must be finite and positive");
  }
}
void check(double           radius,
           double           wavelength,
           double           aspect,
           double           real,
           double           imag,
           double           accuracy,
           EquivalentRadius equivalent_radius,
           ParticleShape    shape,
           Index            degree,
           const Vector&    coefficients) {
  check_shape(aspect, shape, degree, coefficients);
  positive(radius, "radius must be finite and positive");
  positive(wavelength, "wavelength must be finite and positive");
  positive(real, "refractive_real must be finite and positive");
  require(std::isfinite(imag) && imag >= 0, "refractive_imag must be finite and nonnegative");
  positive(accuracy, "accuracy must be finite and positive");
  require(good_enum(equivalent_radius), "equivalent_radius must be an EquivalentRadius");
}
#ifdef ARTS_HAS_TMATRIX

//! The volume-equivalent over the given equivalent radius
Numeric radius_factor(const solver::Particle& particle, EquivalentRadius equivalent_radius, Numeric accuracy) {
  return equivalent_radius == EquivalentRadius::surface_area ? solver::volume_to_surface_radius(particle, accuracy)
                                                             : 1.0;
}

//! The particle of a given volume-equivalent radius
solver::Particle particle(Numeric       radius,
                          Numeric       wavelength,
                          Numeric       aspect_ratio,
                          Numeric       refractive_real,
                          Numeric       refractive_imag,
                          ParticleShape shape,
                          Index         degree,
                          const Vector& coefficients) {
  return {.shape            = shape,
          .aspect_ratio     = aspect_ratio,
          .radius           = radius,
          .wavelength       = wavelength,
          .refractive_index = Complex{refractive_real, refractive_imag},
          .degree           = degree,
          .coefficients     = shape == ParticleShape::generalized_chebyshev
                                  ? solver::volume_normalized_coefficients(coefficients)
                                  : Vector{}};
}

#endif
}  // namespace

bool available() {
#ifdef ARTS_HAS_TMATRIX
  return true;
#else
  return false;
#endif
}

Numeric maximum_dimension(
    Numeric radius, Numeric aspect_ratio, ParticleShape shape, Index degree, const Vector& coefficients) {
  check_shape(aspect_ratio, shape, degree, coefficients);
  positive(radius, "radius must be finite and positive");
  if (shape == ParticleShape::spheroid) {
    const Numeric c = std::cbrt(aspect_ratio);
    return 2 * radius * std::max(c, 1.0 / (c * c));
  }
  if (shape == ParticleShape::cylinder) {
    const Numeric length = 2 * radius * std::cbrt(2.0 / (3.0 * aspect_ratio * aspect_ratio));
    return std::hypot(length, aspect_ratio * length);
  }
#ifdef ARTS_HAS_TMATRIX
  return solver::maximum_dimension(particle(radius, 1.0, aspect_ratio, 1.0, 0.0, shape, degree, coefficients));
#else
  return 0;
#endif
}

FixedResult fixed(Numeric          radius,
                  Numeric          wavelength,
                  Numeric          aspect_ratio,
                  Numeric          refractive_real,
                  Numeric          refractive_imag,
                  Numeric          theta_incident,
                  Numeric          theta_scattered,
                  Numeric          phi_incident,
                  Numeric          phi_scattered,
                  Numeric          alpha,
                  Numeric          beta,
                  Numeric          accuracy,
                  EquivalentRadius equivalent_radius,
                  ParticleShape    shape,
                  Index            degree,
                  const Vector&    coefficients,
                  Index            surface_quadrature) {
  const Geometry geometry{.theta_incident  = theta_incident,
                          .theta_scattered = theta_scattered,
                          .phi_incident    = phi_incident,
                          .phi_scattered   = phi_scattered,
                          .alpha           = alpha,
                          .beta            = beta};
  FixedResult    result;
  fixed_batch(std::span{&result, 1},
              radius,
              wavelength,
              aspect_ratio,
              refractive_real,
              refractive_imag,
              std::span{&geometry, 1},
              accuracy,
              equivalent_radius,
              shape,
              degree,
              coefficients,
              surface_quadrature);
  return result;
}

void fixed_batch(std::span<FixedResult>    results,
                 Numeric                   radius,
                 Numeric                   wavelength,
                 Numeric                   aspect_ratio,
                 Numeric                   refractive_real,
                 Numeric                   refractive_imag,
                 std::span<const Geometry> geometries,
                 Numeric                   accuracy,
                 EquivalentRadius          equivalent_radius,
                 ParticleShape             shape,
                 Index                     degree,
                 const Vector&             coefficients,
                 Index                     surface_quadrature) {
  check(radius,
        wavelength,
        aspect_ratio,
        refractive_real,
        refractive_imag,
        accuracy,
        equivalent_radius,
        shape,
        degree,
        coefficients);
  require(surface_quadrature >= 1, "surface_quadrature must be at least 1");
  require(results.size() == geometries.size(), "results must have one element per geometry");
  for (const Geometry& g : geometries) {
    for (Numeric angle : {g.theta_incident, g.theta_scattered, g.beta})
      require(std::isfinite(angle) && angle >= 0 && angle <= 180, "zenith angles and beta must be in [0, 180] degrees");
    for (Numeric angle : {g.phi_incident, g.phi_scattered, g.alpha})
      require(std::isfinite(angle) && angle >= 0 && angle <= 360,
              "azimuth angles and alpha must be in [0, 360] degrees");
  }
  if (geometries.empty()) return;
#ifdef ARTS_HAS_TMATRIX
  auto shaped =
      particle(radius, wavelength, aspect_ratio, refractive_real, refractive_imag, shape, degree, coefficients);
  shaped.radius                       *= radius_factor(shaped, equivalent_radius, accuracy);
  const auto T                         = solver::solve(shaped, accuracy, surface_quadrature, false);
  const auto [extinction, scattering]  = solver::cross_sections(T, wavelength, false);

  FixedResult out;
  out.order      = T.nmax;
  out.scattering = scattering;
  out.extinction = extinction;
  // Mishchenko's S is the Jones matrix in ARTS's bases (Amplitude matrix in a fixed orientation)
  const auto amplitudes = solver::amplitude(T, wavelength, geometries);
  for (Size i = 0; i < geometries.size(); ++i) {
    out.amplitude = amplitudes[i].S;
    out.phase     = rtepack::to_muelmat(out.amplitude);
    results[i]    = out;
  }
#endif
}

std::vector<FixedResult> fixed_batch(Numeric                   radius,
                                     Numeric                   wavelength,
                                     Numeric                   aspect_ratio,
                                     Numeric                   refractive_real,
                                     Numeric                   refractive_imag,
                                     std::span<const Geometry> geometries,
                                     Numeric                   accuracy,
                                     EquivalentRadius          equivalent_radius,
                                     ParticleShape             shape,
                                     Index                     degree,
                                     const Vector&             coefficients,
                                     Index                     surface_quadrature) {
  std::vector<FixedResult> results(geometries.size());
  fixed_batch(results,
              radius,
              wavelength,
              aspect_ratio,
              refractive_real,
              refractive_imag,
              geometries,
              accuracy,
              equivalent_radius,
              shape,
              degree,
              coefficients,
              surface_quadrature);
  return results;
}

RandomResult random(Numeric                 radius,
                    Numeric                 wavelength,
                    Numeric                 aspect_ratio,
                    Numeric                 refractive_real,
                    Numeric                 refractive_imag,
                    Index                   angles,
                    Numeric                 accuracy,
                    EquivalentRadius        equivalent_radius,
                    ParticleShape           shape,
                    Index                   degree,
                    const Vector&           coefficients,
                    TMatrixSizeDistribution distribution,
                    Numeric                 b,
                    Numeric                 gamma,
                    Index                   size_quadrature,
                    Index                   surface_quadrature,
                    Numeric                 lower_radius_ratio,
                    Numeric                 upper_radius_ratio) {
  check(radius,
        wavelength,
        aspect_ratio,
        refractive_real,
        refractive_imag,
        accuracy,
        equivalent_radius,
        shape,
        degree,
        coefficients);
  require(angles == 0 || angles >= 2, "angles must be 0 or at least 2");
  using Distribution = TMatrixSizeDistribution;
  require(good_enum(distribution), "distribution must be a TMatrixSizeDistribution");
  if (distribution == Distribution::modified_power_law)
    require(std::isfinite(b), "b must be finite");
  else
    positive(b, "b must be finite and positive");
  positive(gamma, "gamma must be finite and positive");
  require(size_quadrature >= -1, "size_quadrature must be at least -1");
  if (distribution == Distribution::modified_power_law)
    require(size_quadrature >= 1, "modified power law requires size_quadrature >= 1");
  require(surface_quadrature >= 1, "surface_quadrature must be at least 1");
  positive(lower_radius_ratio, "lower_radius_ratio must be finite and positive");
  require(std::isfinite(upper_radius_ratio) && upper_radius_ratio > lower_radius_ratio,
          "upper_radius_ratio must exceed lower_radius_ratio");
  // The gamma weights are powers 1 / b - 3 of the radius. Normalizing the
  // radius by its length a keeps them in range.
  const Numeric length_scale  = distribution == Distribution::gamma ? radius : 1;
  radius                     /= length_scale;
  wavelength                 /= length_scale;
  RandomResult out;
#ifdef ARTS_HAS_TMATRIX
  // N_k = size_quadrature + 2 Gauss-Legendre nodes, one for size_quadrature = -1 (Size distributions)
  const auto mean = solver::average(
      particle(radius, wavelength, aspect_ratio, refractive_real, refractive_imag, shape, degree, coefficients),
      radius_factor(
          particle(1.0, 1.0, aspect_ratio, 1.0, 0.0, shape, degree, coefficients), equivalent_radius, accuracy),
      {.kind        = distribution,
       .a           = radius,
       .b           = b,
       .gamma       = gamma,
       .lower_ratio = lower_radius_ratio,
       .upper_ratio = upper_radius_ratio,
       .nodes       = std::max<Index>(size_quadrature + 2, 1)},
      accuracy,
      surface_quadrature);
  const auto& sum  = mean.expansion;
  const auto& sums = mean.cross_sections;

  out.effective_radius   = mean.effective_radius * length_scale;
  out.effective_variance = mean.effective_variance;
  out.extinction         = sums.extinction * length_scale * length_scale;
  out.scattering         = sums.scattering * length_scale * length_scale;
  out.albedo             = sums.scattering / sums.extinction;
  out.asymmetry          = sum.alpha1.size() > 1 ? sum.alpha1[1] / 3.0 : 0.0;

  out.phase    = solver::scattering_matrix(sum, angles);
  out.legendre = solver::legendre_series(sum);
#endif
  return out;
}
AzimuthallyRandomResult azimuthally_random(Numeric          radius,
                                           Numeric          wavelength,
                                           Numeric          aspect_ratio,
                                           Numeric          refractive_real,
                                           Numeric          refractive_imag,
                                           const Vector&    za_inc,
                                           const Vector&    delta_aa,
                                           const Vector&    za_scat,
                                           const Vector&    tilt_angles,
                                           const Vector&    tilt_weights,
                                           Numeric          accuracy,
                                           EquivalentRadius equivalent_radius,
                                           ParticleShape    shape,
                                           Index            degree,
                                           const Vector&    coefficients,
                                           Index            surface_quadrature) {
  check(radius,
        wavelength,
        aspect_ratio,
        refractive_real,
        refractive_imag,
        accuracy,
        equivalent_radius,
        shape,
        degree,
        coefficients);
  require(surface_quadrature >= 1, "surface_quadrature must be at least 1");
  require(not za_inc.empty() and not delta_aa.empty() and not za_scat.empty(), "the angle grids must not be empty");
  for (const Vector* grid : {&za_inc, &za_scat})
    for (Numeric za : *grid)
      require(std::isfinite(za) && za >= 0 && za <= 180, "zenith angles must be in [0, 180] degrees");
  for (Numeric aa : delta_aa) require(std::isfinite(aa), "delta_aa must be finite");
  if (tilt_angles.empty() or tilt_angles.size() != tilt_weights.size())
    throw std::invalid_argument(
        std::format("tilt_angles has {} elements and tilt_weights {}; they must have the same, nonzero, number",
                    tilt_angles.size(),
                    tilt_weights.size()));
  Numeric total = 0.0;
  for (Size i = 0; i < tilt_angles.size(); ++i) {
    require(std::isfinite(tilt_angles[i]) && tilt_angles[i] >= 0 && tilt_angles[i] <= 180,
            "tilt angles must be in [0, 180] degrees");
    require(std::isfinite(tilt_weights[i]) && tilt_weights[i] >= 0, "tilt weights must be finite and nonnegative");
    total += tilt_weights[i];
  }
  require(total > 0, "the tilt weights must not all be zero");

  AzimuthallyRandomResult out;
#ifdef ARTS_HAS_TMATRIX
  auto shaped =
      particle(radius, wavelength, aspect_ratio, refractive_real, refractive_imag, shape, degree, coefficients);
  shaped.radius *= radius_factor(shaped, equivalent_radius, accuracy);
  const auto T   = solver::solve(shaped, accuracy, surface_quadrature, false);
  Vector     weights{tilt_weights};
  weights        /= total;
  auto aro        = solver::azimuthally_random(T, wavelength, tilt_angles, weights, za_inc, delta_aa, za_scat);
  out.order       = T.nmax;
  out.phase       = std::move(aro.phase);
  out.extinction  = std::move(aro.extinction);
  out.absorption  = std::move(aro.absorption);
#endif
  return out;
}
}  // namespace tmatrix
