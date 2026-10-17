/* The T-matrix of a sphere against the Mie solution of src/core/scattering/mie.h

   For a sphere the T-matrix is the diagonal of Mie coefficients, so the two
   independent codes must agree up to the T-matrix's own truncation.  The
   spheres are spheroids of aspect ratio 1 and Chebyshev particles of
   deformation 0, of even degree (integrated over a hemisphere) and of odd
   degree (integrated over the sphere):

   cross sections, the random-orientation scattering matrix, whose expansion
   series is complete, and the fixed-orientation phase matrix, to the requested
   accuracy, which bounds the truncation of the degrees.
*/
#include <mie.h>
#include <rtepack.h>
#include <tmatrix.h>
#include <wigner_functions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <print>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
constexpr Numeric lambda   = 6.283185307179586;
constexpr Numeric accuracy = 1e-9;
constexpr Numeric pi       = std::numbers::pi;

int failures = 0;

void require(bool ok, const std::string& what) {
  if (not ok) {
    std::println("FAIL {}", what);
    ++failures;
  }
}

Numeric relative(Numeric a, Numeric b) { return std::abs(a - b) / std::abs(b); }

void sphere(Numeric r, Complex m, ParticleShape shape, Index degree) {
  const Numeric     aspect = shape == ParticleShape::chebyshev ? 0.0 : 1.0;
  const std::string name   = std::format("sphere r={} m=({}, {}) {} {}", r, m.real(), m.imag(), shape, degree);

  // Random orientation: the six elements of F C_sca / (4 pi) at the scattering angles
  constexpr Index angles = 37;
  const auto      t =
      tmatrix::random(r, lambda, aspect, m.real(), m.imag(), angles, accuracy, EquivalentRadius::volume, shape, degree);
  Vector theta(angles);
  for (Index i = 0; i < angles; ++i) theta[i] = 180.0 * static_cast<Numeric>(i) / (angles - 1);
  scattering::MieSphere<Numeric> mie(lambda, r, m, theta);
  const auto                     F = mie.get_scattering_matrix_compact();

  require(relative(t.extinction, mie.get_extinction_coeff()) <= accuracy, name + ": extinction");
  require(relative(t.scattering, mie.get_scattering_coeff()) <= accuracy, name + ": scattering");

  Numeric f11max = 0.0;
  for (Index i = 0; i < angles; ++i) f11max = std::max(f11max, F[i, 0]);
  const Numeric tol_random = accuracy * f11max;

  // F11, F12, F22, F33, F34, F44 of the T-matrix and of Mie
  Numeric       worst = 0.0;
  const Numeric scale = t.scattering / (4.0 * pi);
  for (Index i = 0; i < angles; ++i)
    for (Index e = 0; e < 6; ++e) worst = std::max(worst, std::abs(t.phase[i][e] * scale - F[i, e]) / tol_random);
  require(worst <= 1.0,
          std::format("{}: random-orientation scattering matrix at {:.2f} of its tolerance", name, worst));

  // Fixed orientation: the Mueller matrix of the T-matrix's Jones matrix is the laboratory-frame
  // phase matrix that ARTS builds from Mie's F (rtepack::lab_frame_phase_matrix), element by
  // element, for any orientation and geometry: two codes and one Stokes convention.  The
  // geometries include the poles, exact backscattering and beta = 90 degrees.  At a pole, the
  // T-matrix's basis is the limit along the given azimuth, h = -phi = (sin phi, -cos phi, 0).
  using tmatrix::Geometry;
  const std::array g{
      Geometry{.theta_incident  = 56,
               .theta_scattered = 65,
               .phi_incident    = 114,
               .phi_scattered   = 128,
               .alpha           = 145,
               .beta            = 52},
      Geometry{
          .theta_incident = 0, .theta_scattered = 180, .phi_incident = 0, .phi_scattered = 0, .alpha = 0, .beta = 0},
      Geometry{.theta_incident  = 30,
               .theta_scattered = 100,
               .phi_incident    = 10,
               .phi_scattered   = 250,
               .alpha           = 20,
               .beta            = 90},
      Geometry{.theta_incident  = 90,
               .theta_scattered = 90,
               .phi_incident    = 0,
               .phi_scattered   = 180,
               .alpha           = 300,
               .beta            = 30}};
  const auto direction = [](Numeric t, Numeric p) {
    return Vector3{std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t)};
  };
  const auto                               pole = [](Numeric p) { return Vector3{std::sin(p), -std::cos(p), 0.0}; };
  Vector                                   scattering_angle(g.size());
  std::vector<rtepack::polarization_basis> basis_in(g.size()), basis_out(g.size());
  for (Size i = 0; i < g.size(); ++i) {
    const Numeric d      = pi / 180.0;
    const Numeric theta0 = g[i].theta_incident * d, theta1 = g[i].theta_scattered * d;
    const Numeric phi0 = g[i].phi_incident * d, phi1 = g[i].phi_scattered * d;
    const Vector3 k_in = direction(theta0, phi0), k_out = direction(theta1, phi1);
    basis_in[i]         = rtepack::meridional_basis(k_in, pole(phi0));
    basis_out[i]        = rtepack::meridional_basis(k_out, pole(phi1));
    const Numeric cos   = k_in[0] * k_out[0] + k_in[1] * k_out[1] + k_in[2] * k_out[2];
    scattering_angle[i] = std::acos(std::clamp(cos, -1.0, 1.0)) / d;
  }
  scattering::MieSphere<Numeric> mie_fixed(lambda, r, m, scattering_angle);
  const auto                     Ff = mie_fixed.get_scattering_matrix_compact();
  const auto                     fixed =
      tmatrix::fixed_batch(r, lambda, aspect, m.real(), m.imag(), g, accuracy, EquivalentRadius::volume, shape, degree);
  require(relative(fixed[0].extinction, mie.get_extinction_coeff()) <= accuracy, name + ": fixed extinction");
  require(relative(fixed[0].scattering, mie.get_scattering_coeff()) <= accuracy, name + ": fixed scattering");
  Numeric worst_fixed = 0.0;
  for (Index i = 0; i < 4; ++i) {
    const rtepack::compact_planar_muelmat F_mie{Ff[i, 0], Ff[i, 1], Ff[i, 2], Ff[i, 3], Ff[i, 4], Ff[i, 5]};
    const Muelmat                         Z = rtepack::lab_frame_phase_matrix(F_mie, basis_in[i], basis_out[i]);
    for (Index a = 0; a < 4; ++a)
      for (Index b = 0; b < 4; ++b)
        worst_fixed = std::max(worst_fixed, std::abs(fixed[i].phase[a, b] - Z[a, b]) / (accuracy * f11max));
  }
  require(worst_fixed <= 1.0,
          std::format("{}: fixed-orientation phase matrix at {:.2f} of its tolerance", name, worst_fixed));

  std::println("{}: random {:.3f} and fixed {:.3f} of their tolerances", name, worst, worst_fixed);
}
}  // namespace

int main() try {
  make_wigner_ready(250, 20000000, 6);
  for (Numeric r : {0.3, 2.0, 6.0})
    for (Complex m : {Complex{1.33, 1e-3}, Complex{1.78, 0.01}, Complex{1.5, 0.5}})
      for (auto [shape, degree] : {std::pair{ParticleShape::spheroid, Index{0}},
                                   std::pair{ParticleShape::chebyshev, Index{2}},
                                   std::pair{ParticleShape::chebyshev, Index{3}}})
        sphere(r, m, shape, degree);
  if (failures) throw std::runtime_error(std::format("{} comparisons failed", failures));
} catch (const std::exception& e) {
  std::println(stderr, "{}", e.what());
  return 1;
}
