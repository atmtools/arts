/* The maximum dimension of the T-matrix shapes

   The numerical maximum of the distance between two surface points
   (tmatrix_solver.h) against the closed forms of the interface (tmatrix.h):
   the longer axis of a spheroid and the diagonal of a cylinder.  The
   spheroid's maximum is smooth, the cylinder's is at the kink of its edge,
   which the one-dimensional maximizations locate only to their own precision.

   The volume- over the surface-equivalent radius of a Chebyshev particle
   against Gauss-Legendre rules of fixed order: their 200- and 400-point
   areas agree to rounding, so either is the area, and the volume is a
   polynomial of degree 3n in cos(theta), which they integrate exactly.  An
   accuracy whose error the area's quadrature cannot estimate throws.

   The closed-form radius ratio of spheroids against the same rules for the
   area of their surface r(theta), also for the sphere, where the closed forms
   are 0 / 0, and next to it, where they would cancel.
*/
#include <legendre.h>
#include <tmatrix.h>
#include <tmatrix_solver.h>

#include <cmath>
#include <limits>
#include <print>
#include <stdexcept>

namespace {
int failures = 0;

void compare(ParticleShape shape, Index degree, Numeric eps, Numeric tolerance) {
  const tmatrix::solver::Particle p{
      .shape = shape, .aspect_ratio = eps, .radius = 1.0, .wavelength = 1.0, .refractive_index = 1.0, .degree = degree};
  const Numeric numeric = tmatrix::solver::maximum_dimension(p);
  // A Chebyshev particle of deformation 0 is a sphere
  const Numeric closed = shape == ParticleShape::chebyshev ? 2.0
                         : shape == ParticleShape::spheroid
                             ? tmatrix::maximum_dimension(1.0, tmatrix::Shape::spheroid(eps))
                             : tmatrix::maximum_dimension(1.0, tmatrix::Shape::cylinder(eps));
  const Numeric error  = std::abs(numeric / closed - 1.0);
  std::println("{} {} eps={}: maximum dimension {:.15f}, relative error {:.1e}", shape, degree, eps, numeric, error);
  if (error > tolerance) {
    std::println("FAIL {} {} eps={}: {} against {}", shape, degree, eps, numeric, closed);
    ++failures;
  }
}

//! The area over 2 pi and the volume over 4 pi / 3 of 1 + eps cos(n theta), by an N-point Gauss-Legendre rule
std::pair<Numeric, Numeric> chebyshev_integrals(Index n, Numeric eps, Index N) {
  Vector x(N), w(N);
  Legendre::GaussLegendre(x, w);
  Numeric area = 0.0, volume = 0.0;
  for (Index i = 0; i < N; ++i) {
    const Numeric nt   = static_cast<Numeric>(n) * std::acos(x[i]);
    const Numeric a    = 1.0 + eps * std::cos(nt);
    const Numeric ens  = eps * static_cast<Numeric>(n) * std::sin(nt);
    area              += w[i] * a * std::sqrt(a * a + ens * ens);
    volume            += 0.5 * w[i] * a * a * a;
  }
  return {area, volume};
}

void surface_radius(Index n, Numeric eps) {
  const tmatrix::solver::Particle p{.shape            = ParticleShape::chebyshev,
                                    .aspect_ratio     = eps,
                                    .radius           = 1.0,
                                    .wavelength       = 1.0,
                                    .refractive_index = 1.0,
                                    .degree           = n};
  const auto [area, volume] = chebyshev_integrals(n, eps, 400);
  const Numeric change      = std::abs(chebyshev_integrals(n, eps, 200).first / area - 1.0);
  // q = r_v / r_s with r_v^3 the volume and r_s^2 half the area over 2 pi; the accuracy 1e-12 asks for the area to
  // 0.1 times it, which bounds the error of q
  const Numeric q      = std::cbrt(volume) / std::sqrt(0.5 * area);
  const Numeric result = tmatrix::solver::volume_to_surface_radius(p, 1e-12);
  const Numeric error  = std::abs(result / q - 1.0);
  std::println("chebyshev {} eps={}: radius ratio {:.15f}, relative error {:.1e}, quadratures agree to {:.1e}",
               n,
               eps,
               result,
               error,
               change);
  if (not(change <= 16 * std::numeric_limits<Numeric>::epsilon() and error <= 1e-13)) {
    std::println("FAIL chebyshev {} eps={}: radius ratio {} against {}", n, eps, result, q);
    ++failures;
  }
  try {
    tmatrix::solver::volume_to_surface_radius(p, 1e-15);
    std::println("FAIL chebyshev {} eps={}: an area to 1e-16 did not throw", n, eps);
    ++failures;
  } catch (const std::exception&) {}
}

/* A generalized Chebyshev particle of the single coefficient c_n = eps is the Chebyshev particle of degree n: the
   surface agrees to rounding, the Gauss-Legendre volume normalization against the closed form, and the maximum
   dimension and radius ratio to the precision of their searches and quadratures.  Its normalized shape has unit
   volume by a far larger Gauss-Legendre rule, and only odd coefficients break its mirror symmetry. */
void generalized_chebyshev(Index n, Numeric eps) {
  namespace solver = tmatrix::solver;
  Vector c(n + 1, 0.0);
  c[n] = eps;
  const solver::Particle general{.shape            = ParticleShape::generalized_chebyshev,
                                 .aspect_ratio     = 0.0,
                                 .radius           = 1.0,
                                 .wavelength       = 1.0,
                                 .refractive_index = 1.0,
                                 .coefficients     = solver::volume_normalized_coefficients(c)};
  const solver::Particle chebyshev{.shape            = ParticleShape::chebyshev,
                                   .aspect_ratio     = eps,
                                   .radius           = 1.0,
                                   .wavelength       = 1.0,
                                   .refractive_index = 1.0,
                                   .degree           = n};
  constexpr Numeric      ulp     = std::numeric_limits<Numeric>::epsilon();
  Numeric                surface = 0.0;
  for (Numeric x : {-0.97, -0.5, 0.0, 0.3, 0.81}) {
    const auto [r2g, rhog] = solver::surface(general, x);
    const auto [r2c, rhoc] = solver::surface(chebyshev, x);
    surface                = std::max({surface, std::abs(r2g / r2c - 1.0), std::abs(rhog - rhoc)});
  }
  const Numeric dmax  = std::abs(solver::maximum_dimension(general) / solver::maximum_dimension(chebyshev) - 1.0);
  const Numeric ratio = std::abs(
      solver::volume_to_surface_radius(general, 1e-12) / solver::volume_to_surface_radius(chebyshev, 1e-12) - 1.0);
  Vector x(200), w(200);
  Legendre::GaussLegendre(x, w);
  Numeric volume = 0.0;
  for (Index i = 0; i < 200; ++i) {
    Numeric a = 0.0;
    for (Index k = 0; k <= n; ++k) a += general.coefficients[k] * std::cos(static_cast<Numeric>(k) * std::acos(x[i]));
    volume += 0.5 * w[i] * a * a * a;
  }
  const bool mirror = solver::mirror_symmetric(general) == (n % 2 == 0);
  std::println("generalized chebyshev c_{}={}: surface {:.1e}, d_max {:.1e}, radius ratio {:.1e}, volume {:.1e}",
               n,
               eps,
               surface,
               dmax,
               ratio,
               std::abs(volume - 1.0));
  if (not(surface <= 16 * ulp and dmax <= 1e-13 and ratio <= 1e-13 and std::abs(volume - 1.0) <= 16 * ulp and mirror)) {
    std::println("FAIL generalized chebyshev c_{}={}", n, eps);
    ++failures;
  }
}

//! The closed-form radius ratio of a spheroid against Gauss-Legendre rules for its area over 2 pi, the integral of
//! r^2 sqrt(1 + rho^2) over cos(theta), for the volume-equivalent radius 1
void spheroid_radius(Numeric eps) {
  const tmatrix::solver::Particle p{
      .shape = ParticleShape::spheroid, .aspect_ratio = eps, .radius = 1.0, .wavelength = 1.0, .refractive_index = 1.0};
  const auto area = [&p](Index N) {
    Vector x(N), w(N);
    Legendre::GaussLegendre(x, w);
    Numeric sum = 0.0;
    for (Index i = 0; i < N; ++i) {
      const auto [r2, rho]  = tmatrix::solver::surface(p, x[i]);
      sum                  += w[i] * r2 * std::sqrt(1.0 + rho * rho);
    }
    return sum;
  };
  const Numeric fine = area(400), change = std::abs(area(200) / fine - 1.0);
  const Numeric q = 1.0 / std::sqrt(0.5 * fine), result = tmatrix::solver::volume_to_surface_radius(p, 1e-12);
  const Numeric error = std::abs(result / q - 1.0);
  std::println("spheroid eps={:.17g}: radius ratio {:.15f}, relative error {:.1e}, quadratures agree to {:.1e}",
               eps,
               result,
               error,
               change);
  if (not(change <= 16 * std::numeric_limits<Numeric>::epsilon() and
          error <= 16 * std::numeric_limits<Numeric>::epsilon())) {
    std::println("FAIL spheroid eps={}: radius ratio {} against {}", eps, result, q);
    ++failures;
  }
}
}  // namespace

int main() try {
  for (Numeric eps : {0.3, 0.7, 1.0, 1.5, 3.0}) compare(ParticleShape::spheroid, 0, eps, 1e-13);
  for (Numeric eps : {0.3, 1.0, 3.0}) compare(ParticleShape::cylinder, 0, eps, 1e-7);
  for (Index degree : {2, 3}) compare(ParticleShape::chebyshev, degree, 0.0, 1e-13);
  for (Index degree : {2, 3, 8})
    for (Numeric eps : {-0.15, 0.1}) surface_radius(degree, eps);
  for (Index degree : {2, 3, 6})
    for (Numeric eps : {-0.15, 0.1}) generalized_chebyshev(degree, eps);
  constexpr Numeric ulp = std::numeric_limits<Numeric>::epsilon();
  for (Numeric eps : {0.5, 1.0 - ulp / 2, 1.0, 1.0 + ulp, 1.0 + 1e-8, 1.7}) spheroid_radius(eps);
  if (failures) throw std::runtime_error(std::format("{} comparisons failed", failures));
} catch (const std::exception& e) {
  std::println(stderr, "{}", e.what());
  return 1;
}
