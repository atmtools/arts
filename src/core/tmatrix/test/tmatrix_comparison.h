#pragma once

/* The C++ T-matrix against outputs of the original Fortran it ports

   The two compute the same equations (doc/arts/concept.tmatrix.rst) in double
   precision, except that the Fortran stores the T-matrix in fixed orientation,
   and the T-matrix and the intermediate arrays of the expansion in random
   orientation, in single precision; so do its extended-precision variants.
   The tolerances are therefore

   - cross sections: rounding, 16 epsilon, times n_max for the replaced
     Gauss-Legendre nodes, times the componentwise (Skeel) condition number
     of Q, as both perturb Q elementwise and its elements span many orders of
     magnitude;
   - amplitudes: the single-precision unit roundoff 2^-24 times the sum of the
     magnitudes of the terms of the amplitude series, plus the cross-section term;
   - expansion coefficients: 2^-24 relative, which the stored intermediates
     give each, plus the terms the Fortran drops (below);
   - scattering matrices: 2^-24 times the sum of the magnitudes of the
     expansion coefficients, which bounds the Legendre sums, plus the terms the
     Fortran drops: it ends the series of each size at the first l with
     |c_l g_1| = |alpha_1 + alpha_4| / 2 < 1e-6, where ARTS keeps all of them,
     and every |P^l_{mn}| <= 1.

   The Fortran moves the angles of an amplitude matrix 1e-7 rad off the
   singular directions, and the port does not, so the port is given the moved
   angles.

   The outputs are recorded from the Fortran as the data of
   tmatrix_fortran_cases.h, or printed by the original programs
   (tmatrix-reference-test).
*/
#include <arts_constants.h>
#include <arts_conversions.h>
#include <lin_alg.h>
#include <tmatrix.h>
#include <tmatrix_solver.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace tmatrix::comparison {
constexpr Numeric epsilon  = std::numeric_limits<Numeric>::epsilon();
constexpr Numeric single   = 0x1p-24;
constexpr Numeric rounding = 16 * epsilon;

inline int failures = 0;

inline void require(bool ok, const std::string& what) {
  if (not ok) {
    std::println("FAIL {}", what);
    ++failures;
  }
}

inline Numeric relative(Numeric a, Numeric b) { return std::abs(a - b) / std::abs(b); }

//! The geometry the Fortran evaluates: polar angles moved 1e-7 rad towards pi / 2, azimuths towards pi, and beta
//! away from pi / 2 when within 1e-7 of it
inline Geometry fortran_geometry(Geometry g) {
  constexpr Numeric pi = Constant::pi, eps = 1e-7;
  const auto        polar   = [](Numeric t) { return t < pi / 2 ? t + eps : (t > pi / 2 ? t - eps : t); };
  const auto        azimuth = [](Numeric p) { return p < pi ? p + eps : (p > pi ? p - eps : p); };
  const auto        beta    = [](Numeric b) {
    if (b <= pi / 2 and pi / 2 - b <= eps) return b - eps;
    if (b > pi / 2 and b - pi / 2 <= eps) return b + eps;
    return b;
  };
  const auto move = [](Numeric& angle, auto&& f) { angle = Conversion::rad2deg(f(Conversion::deg2rad(angle))); };
  move(g.theta_incident, polar);
  move(g.theta_scattered, polar);
  move(g.phi_incident, azimuth);
  move(g.phi_scattered, azimuth);
  move(g.beta, beta);
  return g;
}

//! The componentwise condition number || |Q^{-1}| |Q| ||_inf of the order m at the converged truncation
inline Numeric condition(const solver::Particle& p, const solver::TMatrix& T, Index m) {
  const solver::Surface surface(p, T.nmax, T.nodes);
  const auto [Q, RgQ] = solver::order_q(surface, m);
  ComplexMatrix Qinv(Q.shape());
  inv(Qinv, Q);
  const Index n   = Q.nrows();
  Numeric     out = 0.0;
  for (Index i = 0; i < n; ++i) {
    Numeric row = 0.0;
    for (Index j = 0; j < n; ++j)
      for (Index k = 0; k < n; ++k) row += std::abs(Qinv[i, k]) * std::abs(Q[k, j]);
    out = std::max(out, row);
  }
  return out;
}

//! A shape of the interface
struct Shape {
  ParticleShape shape;
  Index         degree{0};
  Vector        coefficients{};
};

//! The particle of a case: of volume-equivalent radius 1 when the radius is the surface-equivalent one
struct Particle {
  Numeric          radius, wavelength, aspect_ratio;
  Complex          refractive_index;
  Shape            shape;
  EquivalentRadius equivalent;
  Numeric          accuracy;
  Index            surface_quadrature;

  [[nodiscard]] solver::Particle solver_particle() const {
    return {.shape            = shape.shape,
            .aspect_ratio     = aspect_ratio,
            .radius           = 0.0,
            .wavelength       = wavelength,
            .refractive_index = refractive_index,
            .degree           = shape.degree,
            .coefficients =
                shape.coefficients.empty() ? Vector{} : solver::volume_normalized_coefficients(shape.coefficients)};
  }

  //! The volume-equivalent over the given equivalent radius, as tmatrix.h selects it
  [[nodiscard]] Numeric factor() const {
    return equivalent == EquivalentRadius::surface_area ? solver::volume_to_surface_radius(solver_particle(), accuracy)
                                                        : 1.0;
  }

  [[nodiscard]] std::string describe() const {
    return std::format("r={} eps={} m=({}, {}) {} {} {}",
                       radius,
                       aspect_ratio,
                       refractive_index.real(),
                       refractive_index.imag(),
                       shape.shape,
                       shape.degree,
                       equivalent);
  }
};

//! The size distribution of a random-orientation case
struct Distribution {
  TMatrixSizeDistribution kind;
  Numeric                 b, gamma, lower_ratio, upper_ratio;
  Index                   size_quadrature;
};

/** The port at the geometries a reference program evaluates, and the bounds of its difference from that program

  The amplitude bound is per geometry, for every element of S; the phase bound per geometry, for every element of Z.
  Cross is the relative bound of the cross sections. */
struct FixedPort {
  std::vector<FixedResult> results;
  Numeric                  kappa, cross;
  std::vector<Numeric>     amplitude, phase;
};

inline FixedPort fixed_port(const Particle& c, std::span<const Geometry> geometries) {
  FixedPort out;
  out.results.resize(geometries.size());
  fixed_batch(out.results,
              c.radius,
              c.wavelength,
              c.aspect_ratio,
              c.refractive_index.real(),
              c.refractive_index.imag(),
              geometries,
              c.accuracy,
              c.equivalent,
              c.shape.shape,
              c.shape.degree,
              c.shape.coefficients,
              c.surface_quadrature);

  solver::Particle p = c.solver_particle();
  p.radius           = c.factor() * c.radius;
  const auto T       = solver::solve(p, c.accuracy, c.surface_quadrature, false);
  out.kappa          = condition(p, T, 0);
  out.cross          = rounding * static_cast<Numeric>(T.nmax) * out.kappa;
  const auto S       = solver::amplitude(T, c.wavelength, geometries);
  for (Size i = 0; i < geometries.size(); ++i) {
    Numeric smax = 0.0;
    for (Index a = 0; a < 2; ++a)
      for (Index b = 0; b < 2; ++b) smax = std::max(smax, std::abs(out.results[i].amplitude[a, b]));
    // The laboratory rotation mixes at most four particle-frame elements of the shared magnitude
    const Numeric s = 4.0 * single * S[i].magnitude + out.cross * smax;
    out.amplitude.push_back(s);
    // Each Mueller element sums at most four products of two amplitudes
    out.phase.push_back(4.0 * 2.0 * (smax + s) * s);
  }
  return out;
}

//! Compare the port with a reference program's results at the same geometries
inline void compare(const std::string& name, const FixedPort& port, std::span<const FixedResult> ref) {
  require(port.results[0].order == ref[0].order, name + ": n_max");
  require(relative(port.results[0].extinction, ref[0].extinction) <= port.cross, name + ": extinction");
  require(relative(port.results[0].scattering, ref[0].scattering) <= port.cross, name + ": scattering");
  Numeric worst_s = 0.0, worst_z = 0.0;
  for (Size i = 0; i < ref.size(); ++i) {
    for (Index a = 0; a < 2; ++a)
      for (Index b = 0; b < 2; ++b)
        worst_s =
            std::max(worst_s, std::abs(port.results[i].amplitude[a, b] - ref[i].amplitude[a, b]) / port.amplitude[i]);
    for (Index a = 0; a < 4; ++a)
      for (Index b = 0; b < 4; ++b)
        worst_z = std::max(worst_z, std::abs(port.results[i].phase[a, b] - ref[i].phase[a, b]) / port.phase[i]);
  }
  require(worst_s <= 1.0, std::format("{}: amplitude at {:.2f} of its tolerance", name, worst_s));
  require(worst_z <= 1.0, std::format("{}: phase matrix at {:.2f} of its tolerance", name, worst_z));
  std::println("{}: kappa {:.1e}, cross sections {:.3f}, amplitude {:.3f} and phase {:.3f} of their tolerances",
               name,
               port.kappa,
               std::max(relative(port.results[0].extinction, ref[0].extinction),
                        relative(port.results[0].scattering, ref[0].scattering)) /
                   port.cross,
               worst_s,
               worst_z);
}

/** The port in random orientation, and the bounds of its difference from a reference program

  Cross is the relative bound of the cross sections, matrix the absolute bound of the asymmetry and the scattering
  matrix, and tail the absolute bound of the terms the Fortran drops, which bounds them in every coefficient. */
struct RandomPort {
  RandomResult           result;
  solver::Average        mean;
  solver::SizeQuadrature nodes;
  Numeric                cross, matrix, tail, spread;
};

inline RandomPort random_port(const Particle& c, const Distribution& d, Index angles) {
  RandomPort out;
  out.result = tmatrix::random(c.radius,
                               c.wavelength,
                               c.aspect_ratio,
                               c.refractive_index.real(),
                               c.refractive_index.imag(),
                               angles,
                               c.accuracy,
                               c.equivalent,
                               c.shape.shape,
                               c.shape.degree,
                               c.shape.coefficients,
                               d.kind,
                               d.b,
                               d.gamma,
                               d.size_quadrature,
                               c.surface_quadrature,
                               d.lower_ratio,
                               d.upper_ratio);

  // The port's own averages; the gamma distribution is normalized by the radius, as tmatrix.h does
  const Numeric    scale = d.kind == TMatrixSizeDistribution::gamma ? c.radius : 1.0;
  solver::Particle p     = c.solver_particle();
  p.wavelength           = c.wavelength / scale;
  const solver::SizeDistribution sd{
      d.kind, c.radius / scale, d.b, d.gamma, d.lower_ratio, d.upper_ratio, std::max<Index>(d.size_quadrature + 2, 1)};
  const Numeric factor = c.factor();
  out.mean             = solver::average(p, factor, sd, c.accuracy, c.surface_quadrature);
  out.nodes            = solver::size_quadrature(sd);
  out.spread           = max(out.nodes.radius) / out.mean.effective_radius;

  // The conditioning of the largest particle bounds that of the cross sections
  solver::Particle largest = p;
  largest.radius           = factor * max(out.nodes.radius);
  const auto T             = solver::solve(largest, c.accuracy, c.surface_quadrature, true);
  out.cross                = rounding * static_cast<Numeric>(T.nmax) * condition(largest, T, 0);
  Numeric sum              = 0.0;
  for (const auto* x : {&out.mean.expansion.alpha1,
                        &out.mean.expansion.alpha2,
                        &out.mean.expansion.alpha3,
                        &out.mean.expansion.alpha4,
                        &out.mean.expansion.beta1,
                        &out.mean.expansion.beta2})
    for (auto v : *x) sum += std::abs(v);

  // The Fortran's dropped terms, size by size and weighted as the average weights them
  Numeric tail = 0.0, weights = 0.0;
  for (Size i = 0; i < out.nodes.radius.size(); ++i) {
    solver::Particle node = p;
    node.radius           = factor * out.nodes.radius[i];
    const auto  Ti        = solver::solve(node, c.accuracy, c.surface_quadrature, true);
    const auto  ci        = solver::cross_sections(Ti, node.wavelength, true);
    const auto  e         = solver::expansion(Ti, ci.scattering, node.wavelength);
    const Index degrees   = e.alpha1.size();
    Index       first     = 0;
    while (first < degrees and std::abs(0.5 * (e.alpha1[first] + e.alpha4[first])) >= 1e-6) ++first;
    Numeric dropped = 0.0;
    for (Index l = first + 1; l < degrees; ++l)
      dropped += std::abs(e.alpha1[l]) + std::abs(e.alpha2[l]) + std::abs(e.alpha3[l]) + std::abs(e.alpha4[l]) +
                 std::abs(e.beta1[l]) + std::abs(e.beta2[l]);
    tail    += out.nodes.weight[i] * ci.scattering * dropped;
    weights += out.nodes.weight[i] * ci.scattering;
  }
  out.tail   = tail / weights;
  out.matrix = single * sum + out.tail;
  return out;
}

//! The bound of the effective variance: the size quadratures differ by the Gauss-Legendre nodes, which perturb
//! r - r_eff by about 2 epsilon r
inline Numeric variance_bound(const RandomPort& port, Numeric variance) {
  return 4.0 * rounding * (std::sqrt(variance) + rounding) * port.spread * port.spread;
}

//! Compare the port with a reference program's results
inline void compare(const std::string& name, const RandomPort& port, const RandomResult& ref) {
  const auto& r = port.result;
  require(relative(r.extinction, ref.extinction) <= port.cross, name + ": extinction");
  require(relative(r.scattering, ref.scattering) <= port.cross, name + ": scattering");
  require(relative(r.albedo, ref.albedo) <= 2 * port.cross, name + ": albedo");
  require(relative(r.effective_radius, ref.effective_radius) <= rounding, name + ": effective radius");
  require(std::abs(r.effective_variance - ref.effective_variance) <= variance_bound(port, ref.effective_variance),
          name + ": effective variance");
  require(std::abs(r.asymmetry - ref.asymmetry) <= port.matrix, name + ": asymmetry");
  require(r.phase.size() == ref.phase.size(), name + ": scattering angles");
  Numeric worst = 0.0;
  for (Size i = 0; i < std::min(r.phase.size(), ref.phase.size()); ++i)
    for (Index e = 0; e < 6; ++e) worst = std::max(worst, std::abs(r.phase[i][e] - ref.phase[i][e]) / port.matrix);
  require(worst <= 1.0, std::format("{}: scattering matrix at {:.2f} of its tolerance", name, worst));
  std::println("{}: cross sections {:.3f} and scattering matrix {:.3f} of their tolerances",
               name,
               std::max(relative(r.extinction, ref.extinction), relative(r.scattering, ref.scattering)) / port.cross,
               worst);
}

//! The Fortran's results in fixed orientation at the geometries it evaluates
struct FixedRecord {
  struct Point {
    Geometry geometry;
    Jones    amplitude;
  };

  Particle           particle;
  Index              order;
  Numeric            extinction, scattering;
  std::vector<Point> points;

  [[nodiscard]] std::vector<Geometry> geometries() const {
    std::vector<Geometry> out;
    for (const auto& p : points) out.push_back(p.geometry);
    return out;
  }

  //! As the reference interface returns them: the phase matrix is the Mueller matrix of the Fortran's amplitude
  [[nodiscard]] std::vector<FixedResult> results() const {
    std::vector<FixedResult> out;
    for (const auto& p : points)
      out.push_back({.order      = order,
                     .scattering = scattering,
                     .extinction = extinction,
                     .amplitude  = p.amplitude,
                     .phase      = rtepack::to_muelmat(p.amplitude)});
    return out;
  }
};

//! The Fortran's results in random orientation, the scattering matrix in ARTS's Stokes convention
struct RandomRecord {
  Particle                          particle;
  Distribution                      distribution;
  Numeric                           extinction, scattering, albedo, asymmetry, effective_radius, effective_variance;
  std::vector<CompactPlanarMuelmat> phase;

  [[nodiscard]] RandomResult result() const {
    RandomResult out;
    out.extinction         = extinction;
    out.scattering         = scattering;
    out.albedo             = albedo;
    out.asymmetry          = asymmetry;
    out.effective_radius   = effective_radius;
    out.effective_variance = effective_variance;
    out.phase.resize(phase.size());
    for (Size i = 0; i < phase.size(); ++i) out.phase[i] = phase[i];
    return out;
  }
};
}  // namespace tmatrix::comparison
