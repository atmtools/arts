/* The amplitude matrix towards its singular directions

   The amplitude matrix is evaluated at the given angles, also along the
   particle axis, at the laboratory poles and at beta = 90 degrees
   (concept.tmatrix).  Along a smooth path towards such a direction,
   S(delta) = S(0) + delta S' + O(delta^2).  With S' from delta_r and the
   second derivative C from delta_r and 2 delta_r, the residual
   S(delta) - S(0) - delta S' of every smaller delta is at most C delta delta_r,
   plus the rounding of the four evaluations it combines: 16 epsilon times the
   magnitude of the terms of the amplitude series, and |S'| times the rounding
   of an angle of at most 360 degrees in degrees and in radians.  A limit taken
   away from the singular direction, or an angle that loses accuracy next to
   it, leaves a residual of order delta S'.
*/
#include <tmatrix_solver.h>

#include <array>
#include <cmath>
#include <limits>
#include <print>
#include <stdexcept>
#include <vector>

namespace {
using namespace tmatrix;

constexpr Numeric epsilon = std::numeric_limits<Numeric>::epsilon();
constexpr Numeric lambda  = 6.283185307179586;

int failures = 0;

//! A singular geometry, and the angle that moves towards it
struct Approach {
  const char* name;
  Geometry    singular;
  Numeric Geometry::* angle;
  Numeric             sign;
};

constexpr std::array approaches{
    Approach{
        "incident along the axis",
        Geometry{
            .theta_incident = 0, .theta_scattered = 65, .phi_incident = 0, .phi_scattered = 128, .alpha = 0, .beta = 0},
        &Geometry::theta_incident,
        1.0},
    Approach{"scattered against the axis",
             Geometry{.theta_incident  = 56,
                      .theta_scattered = 180,
                      .phi_incident    = 114,
                      .phi_scattered   = 0,
                      .alpha           = 0,
                      .beta            = 0},
             &Geometry::theta_scattered,
             -1.0},
    Approach{"incident on the axis at beta = 90 from below",
             Geometry{.theta_incident  = 90,
                      .theta_scattered = 65,
                      .phi_incident    = 20,
                      .phi_scattered   = 128,
                      .alpha           = 20,
                      .beta            = 90},
             &Geometry::beta,
             -1.0},
    Approach{"incident on the axis at beta = 90 from above",
             Geometry{.theta_incident  = 90,
                      .theta_scattered = 65,
                      .phi_incident    = 20,
                      .phi_scattered   = 128,
                      .alpha           = 20,
                      .beta            = 90},
             &Geometry::beta,
             1.0},
    Approach{"incident at the laboratory pole",
             Geometry{.theta_incident  = 0,
                      .theta_scattered = 65,
                      .phi_incident    = 30,
                      .phi_scattered   = 128,
                      .alpha           = 145,
                      .beta            = 52},
             &Geometry::theta_incident,
             1.0},
    Approach{"scattered at the laboratory pole",
             Geometry{.theta_incident  = 56,
                      .theta_scattered = 180,
                      .phi_incident    = 114,
                      .phi_scattered   = 200,
                      .alpha           = 145,
                      .beta            = 52},
             &Geometry::theta_scattered,
             -1.0},
    Approach{"generic at beta = 90",
             Geometry{.theta_incident  = 56,
                      .theta_scattered = 65,
                      .phi_incident    = 114,
                      .phi_scattered   = 128,
                      .alpha           = 145,
                      .beta            = 90},
             &Geometry::beta,
             1.0},
};

void approach(const solver::Particle& p, const char* particle) {
  const auto           T         = solver::solve(p, 1e-3, 2, false);
  constexpr Numeric    reference = 1e-2;
  std::vector<Numeric> deltas{0.0, reference, 2.0 * reference};
  for (Index k = 3; k <= 12; ++k) deltas.push_back(std::pow(10.0, -static_cast<Numeric>(k)));
  const Index n = static_cast<Index>(deltas.size());

  // The moves are the differences of the rounded angles
  std::vector<Geometry> g;
  std::vector<Numeric>  moved;
  for (const Approach& A : approaches)
    for (const Numeric delta : deltas) {
      Geometry geometry  = A.singular;
      geometry.*A.angle += A.sign * delta;
      moved.push_back(std::abs(geometry.*A.angle - A.singular.*A.angle));
      g.push_back(geometry);
    }
  const auto S = solver::amplitude(T, lambda, g);

  for (Index a = 0; a < static_cast<Index>(approaches.size()); ++a) {
    const Index   base = a * n;
    const Numeric dr = moved[base + 1], d2 = moved[base + 2];
    Numeric       magnitude = 0.0;
    for (Index i = 0; i < n; ++i) magnitude = std::max(magnitude, S[base + i].magnitude);

    std::array<Complex, 4> slope;
    std::array<Numeric, 4> curvature;
    Numeric                largest_slope = 0.0;
    for (Index e = 0; e < 4; ++e) {
      const Index   r = e / 2, c = e % 2;
      const Complex s0 = S[base].S[r, c], sr = S[base + 1].S[r, c], s2 = S[base + 2].S[r, c];
      slope[e]      = (sr - s0) / dr;
      curvature[e]  = std::abs(2.0 * ((s2 - s0) / d2 - slope[e]) / (d2 - dr));
      largest_slope = std::max(largest_slope, std::abs(slope[e]));
    }
    const Numeric rounding = 16.0 * epsilon * magnitude + 2.0 * 360.0 * epsilon * largest_slope;

    Numeric worst = 0.0;
    for (Index i = 3; i < n; ++i)
      for (Index e = 0; e < 4; ++e) {
        const Index   r = e / 2, c = e % 2;
        const Numeric delta    = moved[base + i];
        const Numeric residual = std::abs(S[base + i].S[r, c] - S[base].S[r, c] - delta * slope[e]);
        worst                  = std::max(worst, residual / (curvature[e] * delta * dr + 4.0 * rounding));
      }
    std::println("{}, {}: residual at {:.3f} of its bound", particle, approaches[a].name, worst);
    if (not(worst <= 1.0)) {
      std::println("FAIL {}, {}", particle, approaches[a].name);
      ++failures;
    }
  }
}

solver::Particle particle(ParticleShape shape, Numeric eps, Numeric radius, Index degree) {
  return {.shape            = shape,
          .aspect_ratio     = eps,
          .radius           = radius,
          .wavelength       = lambda,
          .refractive_index = Complex{1.5, 0.02},
          .degree           = degree};
}
}  // namespace

int main() try {
  approach(particle(ParticleShape::spheroid, 0.5, 2.0, 0), "oblate spheroid");
  approach(particle(ParticleShape::cylinder, 2.0, 1.5, 0), "cylinder");
  approach(particle(ParticleShape::chebyshev, 0.1, 2.0, 3), "Chebyshev particle of degree 3");
  if (failures) throw std::runtime_error(std::format("{} approaches failed", failures));
} catch (const std::exception& e) {
  std::println(stderr, "{}", e.what());
  return 1;
}
