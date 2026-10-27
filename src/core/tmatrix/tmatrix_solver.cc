#include "tmatrix_solver.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <arts_omp.h>
#include <debug.h>
#include <legendre.h>
#include <lin_alg.h>
#include <wigner_functions.h>

#include <algorithm>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>

namespace tmatrix::solver {
namespace {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#else
#define WIGNER3 wig3jj
#endif

constexpr Numeric pi = Constant::pi;

Numeric parity(Index n) { return n % 2 == 0 ? 1.0 : -1.0; }

//! i^n
Complex i_power(Index n) {
  switch (((n % 4) + 4) % 4) {
    case 0:  return {1.0, 0.0};
    case 1:  return {0.0, 1.0};
    case 2:  return {-1.0, 0.0};
    default: return {0.0, -1.0};
  }
}

//! C(j1 m1 j2 m2 | J M), Condon-Shortley phase, from wigxjpf; Wigner scratch must be initialized by the caller
Numeric clebsch_gordan(Index j1, Index m1, Index j2, Index m2, Index J, Index M) {
  if (m1 + m2 != M) return 0.0;
  const Numeric w3j = WIGNER3(static_cast<int>(2 * j1),
                              static_cast<int>(2 * j2),
                              static_cast<int>(2 * J),
                              static_cast<int>(2 * m1),
                              static_cast<int>(2 * m2),
                              static_cast<int>(-2 * M));
  return parity(j1 - j2 + M) * std::sqrt(static_cast<Numeric>(2 * J + 1)) * w3j;
}

/** Spherical Bessel functions of the first kind and their tilde forms, n = 1, ..., nmax (Special functions)
 *
 * Downward recurrence of the ratio j_n / j_{n-1} from degree L = ratio.size() - 1,
 * normalized by j_{-1}(z) = cos(z) / z.  ratio is the caller's scratch.
 */
template <typename T> void spherical_bessel_j(T z, Index nmax, auto&& ratio, auto&& j, auto&& dj) {
  const Index L = ratio.size() - 1;
  ratio[L]      = z / static_cast<Numeric>(2 * L + 1);
  for (Index n = L - 1; n >= 1; --n) ratio[n] = T{1.0} / (static_cast<Numeric>(2 * n + 1) / z - ratio[n + 1]);

  const T j0 = std::cos(z) / z / (T{1.0} / z - ratio[1]);
  T       previous{j0};
  for (Index n = 1; n <= nmax; ++n) {
    const T jn = previous * ratio[n];
    j[n]       = jn;
    dj[n]      = previous - static_cast<Numeric>(n) * jn / z;
    previous   = jn;
  }
}

//! Spherical Bessel functions of the second kind and their tilde forms, n = 1, ..., nmax, by upward recurrence
void spherical_bessel_y(Numeric x, Index nmax, auto&& y, auto&& dy) {
  const Numeric c = std::cos(x), s = std::sin(x);
  const Numeric y0 = -c / x;
  y[1]             = -c / (x * x) - s / x;
  y[2]             = (-3.0 / (x * x * x) + 1.0 / x) * c - 3.0 * s / (x * x);
  for (Index n = 2; n < nmax; ++n) y[n + 1] = static_cast<Numeric>(2 * n + 1) / x * y[n] - y[n - 1];

  dy[1] = y0 - y[1] / x;
  for (Index n = 2; n <= nmax; ++n) dy[n] = y[n - 1] - static_cast<Numeric>(n) * y[n] / x;
}

//! d^n_{0m}(theta) and tau^m_n(theta) into [nmax + 1, nmax + 1] from the Schmidt functions (Special functions)
void wigner_d(Numeric theta, MatrixView d, MatrixView dd) {
  Legendre::schmidt(d, dd, theta);
  for (Index n = 1; n < d.nrows(); ++n) {
    for (Index m = 1; m <= n; ++m) {
      d[n, m]  *= Constant::inv_sqrt_2;
      dd[n, m] *= Constant::inv_sqrt_2;
    }
  }
}

//! The block [[T^{11}, T^{12}], [T^{21}, T^{22}]] at the degrees n, n' of the T-matrix of one order, from degree n0
Jones order_block(const ComplexMatrix& Tm, Index n0, Index n, Index nn) {
  const Index N = Tm.nrows() / 2, i = n - n0, j = nn - n0;
  return {Tm[i, j], Tm[i, N + j], Tm[N + i, j], Tm[N + i, N + j]};
}

//! The terms of degree n of the m = 0 sums q_ext and q_sca of the convergence test, from the diagonal of the block
std::pair<Numeric, Numeric> order0_terms(const ComplexMatrix& T0, Index n) {
  const Jones   t  = order_block(T0, 1, n, n);
  const Numeric dn = static_cast<Numeric>(2 * n + 1);
  return {dn * trace(t).real(), dn * (std::norm(t[0, 0]) + std::norm(t[1, 1]))};
}

//! The m = 0 sums q_ext and q_sca of the convergence test (Truncation and convergence)
std::pair<Numeric, Numeric> order0_sums(const ComplexMatrix& T0, Index nmax) {
  Numeric qext = 0.0, qsca = 0.0;
  for (Index n = 1; n <= nmax; ++n) {
    const auto [ext, sca]  = order0_terms(T0, n);
    qext                  += ext;
    qsca                  += sca;
  }
  return {qext, qsca};
}

//! n_1: the first degree from nmax / 2 + 1 whose terms are below delta relative to the sums, else nmax
Index truncation_degree(const ComplexMatrix& T0, Index nmax, Numeric qext, Numeric qsca, Numeric delta) {
  for (Index n = nmax / 2 + 1; n <= nmax; ++n) {
    const auto [ext, sca] = order0_terms(T0, n);
    if (std::abs(sca / qsca) <= delta and std::abs(ext / qext) <= delta) return n;
  }
  return nmax;
}
//! (r_v / r_0)^3 of a Chebyshev particle r_0 (1 + eps cos(n theta)), exactly (Particle shapes)
Numeric chebyshev_volume(const Particle& particle) {
  const Numeric eps = particle.aspect_ratio;
  const Numeric n   = static_cast<Numeric>(particle.degree);
  const Numeric n2 = n * n, e2 = eps * eps;
  Numeric       A = 1.0 + 1.5 * e2 * (4.0 * n2 - 2.0) / (4.0 * n2 - 1.0);
  if (particle.degree % 2 == 0) A = A - 3.0 * eps * (1.0 + 0.25 * e2) / (n2 - 1.0) - 0.25 * e2 * eps / (9.0 * n2 - 1.0);
  return A;
}

//! a = sum_n b_n cos(n theta) and da / dtheta of a generalized Chebyshev shape (Particle shapes)
std::pair<Numeric, Numeric> generalized_shape(const Vector& b, Numeric theta) {
  Numeric a = 0.0, da = 0.0;
  for (Size n = 0; n < b.size(); ++n) {
    const Numeric nt  = static_cast<Numeric>(n) * theta;
    a                += b[n] * std::cos(nt);
    da               -= static_cast<Numeric>(n) * b[n] * std::sin(nt);
  }
  return {a, da};
}

/** The area over 2 pi of the surface r = a(theta), the integral of a sqrt(a^2 + (da/dtheta)^2) over cos(theta), and
 * its error estimate, to delta.  Boost accepts a subinterval at its share of the tolerance, so the estimates of the
 * subintervals add to at most twice it. */
std::pair<Numeric, Numeric> surface_area(auto&& shape, Numeric delta) {
  const auto f = [&shape](Numeric x) {
    const auto [a, da] = shape(std::acos(x));
    return a * std::sqrt(a * a + da * da);
  };
  Numeric       error = 0.0;
  const Numeric area  = boost::math::quadrature::gauss_kronrod<Numeric, 61>::integrate(
      f, -1.0, 1.0, 15 /* Boost's default depth */, 0.5 * delta, &error);
  return {area, error};
}
}  // namespace

Vector volume_normalized_coefficients(const Vector& c) {
  // a(x) = 1 + c_0 + sum_n c_n T_n(x) has degree N in x = cos(theta), so its volume over that of the unit sphere, half
  // the integral of a^3, is exact with ceil((3 N + 1) / 2) Gauss-Legendre nodes
  const Size N = c.size() - 1, nodes = (3 * N + 2) / 2;
  Vector     x(nodes), w(nodes);
  Legendre::GaussLegendre(x, w);
  Numeric volume = 0.0;
  for (Size i = 0; i < nodes; ++i) {
    Numeric t_previous = 1.0, t = x[i], a = 1.0 + c[0];
    for (Size n = 1; n <= N; ++n) {
      a                  += c[n] * t;
      const Numeric next  = 2.0 * x[i] * t - t_previous;
      t_previous          = std::exchange(t, next);
    }
    volume += 0.5 * w[i] * a * a * a;
  }
  Vector        b    = c;
  const Numeric r_v  = std::cbrt(volume);
  b[0]              += 1.0;
  b                 /= r_v;
  return b;
}

bool mirror_symmetric(const Particle& particle) {
  switch (particle.shape) {
    case ParticleShape::spheroid:
    case ParticleShape::cylinder:  return true;
    case ParticleShape::chebyshev: return particle.degree % 2 == 0;
    case ParticleShape::generalized_chebyshev:
      // cos(n (pi - theta)) = (-1)^n cos(n theta): mirror symmetric without odd terms
      for (Size n = 1; n < particle.coefficients.size(); n += 2)
        if (particle.coefficients[n] != 0.0) return false;
      return true;
  }
  std::unreachable();
}

std::pair<Numeric, Numeric> surface(const Particle& particle, Numeric x) {
  const Numeric eps = particle.aspect_ratio;
  // 1 - x and 1 + x are exact where the other cancels, so s is sin(theta) to rounding also next to the poles
  const Numeric s = std::sqrt((1.0 - x) * (1.0 + x));
  switch (particle.shape) {
    case ParticleShape::spheroid: {
      const Numeric a       = particle.radius * std::cbrt(eps);
      const Numeric inverse = 1.0 / (s * s + eps * eps * x * x);
      return {a * a * inverse, (eps * eps - 1.0) * s * x * inverse};
    }
    case ParticleShape::cylinder: {
      // r(pi - theta) = r(theta) and rho(pi - theta) = -rho(theta)
      const Numeric h    = particle.radius * std::cbrt(2.0 / (3.0 * eps * eps));
      const Numeric ax   = std::abs(x);
      const bool    side = s / ax > eps;
      const Numeric r    = side ? eps * h / s : h / ax;
      const Numeric rho  = side ? -ax / s : s / ax;
      return {r * r, x < 0.0 ? -rho : rho};
    }
    case ParticleShape::chebyshev: {
      // r_0 from the volume of r_0 (1 + eps cos(n theta)) and r_v
      const Numeric n  = static_cast<Numeric>(particle.degree);
      const Numeric r0 = particle.radius / std::cbrt(chebyshev_volume(particle));
      const Numeric nt = std::acos(x) * n;
      const Numeric r  = r0 * (1.0 + eps * std::cos(nt));
      return {r * r, -r0 * eps * n * std::sin(nt) / r};
    }
    case ParticleShape::generalized_chebyshev: {
      // The coefficients are those of the volume-equivalent radius 1
      const auto [a, da] = generalized_shape(particle.coefficients, std::acos(x));
      const Numeric r    = particle.radius * a;
      return {r * r, da / a};
    }
  }
  std::unreachable();
}

Numeric volume_to_surface_radius(const Particle& particle, Numeric accuracy) {
  const Numeric eps = particle.aspect_ratio;
  switch (particle.shape) {
    case ParticleShape::cylinder: return std::cbrt(1.5 / eps) / std::sqrt((eps + 2.0) / (2.0 * eps));
    case ParticleShape::spheroid: {
      // A sphere has e = 0, where both expressions are 0 / 0 with the limit 1.  ln((1 + e) / (1 - e)) is 2 atanh(e),
      // which does not cancel as e -> 0.
      if (eps == 1.0) return 1.0;
      const Numeric c = std::cbrt(eps);
      if (eps < 1.0) {
        const Numeric e = std::sqrt((1.0 - eps) * (1.0 + eps));
        return 1.0 / std::sqrt(0.5 * (c * c + std::asin(e) / (c * e)));
      }

      const Numeric e = std::sqrt(1.0 - 1.0 / (eps * eps));
      return 1.0 / std::sqrt(0.25 * (2.0 * c * c + 2.0 * std::atanh(e) / (eps * c * e)));
    }
    case ParticleShape::chebyshev: {
      // The volume of 1 + eps cos(n theta) exactly, and its area to 0.1 accuracy as the convergence criterion of the
      // T-matrix
      const Numeric n          = static_cast<Numeric>(particle.degree);
      const Numeric delta      = 0.1 * accuracy;
      const auto [area, error] = surface_area(
          [eps, n](Numeric theta) {
            return std::pair{1.0 + eps * std::cos(n * theta), -eps * n * std::sin(n * theta)};
          },
          delta);
      if (not(error <= delta * area))
        arts::user_error(
            "The surface area of the Chebyshev particle of degree {} and deformation {} does not converge "
            "to the accuracy {}: its relative error estimate is {}",
            particle.degree,
            eps,
            accuracy,
            error / area);
      return std::cbrt(chebyshev_volume(particle)) / std::sqrt(0.5 * area);
    }
    case ParticleShape::generalized_chebyshev: {
      // The coefficients are those of the volume-equivalent radius 1; the area as for Chebyshev particles
      const Numeric delta = 0.1 * accuracy;
      const auto [area, error] =
          surface_area([&particle](Numeric theta) { return generalized_shape(particle.coefficients, theta); }, delta);
      if (not(error <= delta * area))
        arts::user_error(
            "The surface area of the generalized Chebyshev particle of {} coefficients does not converge "
            "to the accuracy {}: its relative error estimate is {}",
            particle.coefficients.size(),
            accuracy,
            error / area);
      return 1.0 / std::sqrt(0.5 * area);
    }
  }
  std::unreachable();
}

Numeric maximum_dimension(const Particle& particle) {
  // Two points at opposite azimuths and polar angles t1, t2 are sqrt(D2) apart
  const auto point = [&particle](Numeric theta) {
    const Numeric r = std::sqrt(surface(particle, std::cos(theta)).first);
    return std::pair{r * std::sin(theta), r * std::cos(theta)};
  };
  const auto D2 = [&point](Numeric t1, Numeric t2) {
    const auto [rho1, z1] = point(t1);
    const auto [rho2, z2] = point(t2);
    return (rho1 + rho2) * (rho1 + rho2) + (z1 - z2) * (z1 - z2);
  };

  // The surface oscillates degree times over [0, pi], a spheroid or cylinder not at all, and D2 once more in
  // each angle from the sines and cosines.  Four grid steps per half-oscillation put grid points on both sides
  // of every maximum of D2, and every local maximum of the grid is refined.
  const Index   highest = particle.shape == ParticleShape::chebyshev ? particle.degree
                          : particle.shape == ParticleShape::generalized_chebyshev
                              ? static_cast<Index>(particle.coefficients.size()) - 1
                              : 0;
  const Index   half    = highest + 1;
  const Index   K       = 4 * half + 1;
  const Numeric step    = pi / static_cast<Numeric>(K - 1);
  Vector        rho(K), z(K);
  for (Index k = 0; k < K; ++k) std::tie(rho[k], z[k]) = point(step * static_cast<Numeric>(k));
  Matrix grid(K, K);
  for (Index a = 0; a < K; ++a)
    for (Index b = 0; b < K; ++b) grid[a, b] = (rho[a] + rho[b]) * (rho[a] + rho[b]) + (z[a] - z[b]) * (z[a] - z[b]);

  // Maximized in each angle in turn, within a grid step, until D2 no longer grows
  constexpr auto bits   = std::numeric_limits<Numeric>::digits / 2;
  const auto     refine = [&](Numeric t1, Numeric t2) {
    Numeric best = D2(t1, t2);
    for (;;) {
      const Numeric previous = best;
      t1                     = boost::math::tools::brent_find_minima(
                                   [&](Numeric t) { return -D2(t, t2); }, std::max(0.0, t1 - step), std::min(pi, t1 + step), bits)
                                   .first;
      t2                     = boost::math::tools::brent_find_minima(
                                   [&](Numeric t) { return -D2(t1, t); }, std::max(0.0, t2 - step), std::min(pi, t2 + step), bits)
                                   .first;
      best                   = std::max(best, D2(t1, t2));
      if (not(best > previous)) return best;
    }
  };

  Numeric best = 0.0;
  for (Index a = 0; a < K; ++a) {
    for (Index b = a; b < K; ++b) {
      bool local = true;
      for (Index da = -1; da <= 1 and local; ++da)
        for (Index db = -1; db <= 1 and local; ++db) {
          const Index aa = a + da, bb = b + db;
          if (aa >= 0 and aa < K and bb >= 0 and bb < K) local = grid[a, b] >= grid[aa, bb];
        }
      if (local) best = std::max(best, refine(step * static_cast<Numeric>(a), step * static_cast<Numeric>(b)));
    }
  }
  return std::sqrt(best);
}

Surface::Surface(const Particle& particle, Index nmax_, Index hemisphere)
    : refractive_index(particle.refractive_index),
      nmax(nmax_),
      mirror(mirror_symmetric(particle)),
      d(mirror ? hemisphere : 2 * hemisphere, nmax + 1, nmax + 1),
      dd(d.npages(), nmax + 1, nmax + 1),
      x(d.npages()),
      s(x.size()),
      W(x.size()),
      rho(x.size()),
      xi(x.size()),
      j(x.size(), nmax + 1, 0.0),
      dj(x.size(), nmax + 1, 0.0),
      y(x.size(), nmax + 1, 0.0),
      dy(x.size(), nmax + 1, 0.0),
      J(x.size(), nmax + 1, 0.0),
      dJ(x.size(), nmax + 1, 0.0) {
  const Index   nodes = x.size();
  const Numeric k     = 2.0 * pi / particle.wavelength;

  // The quadrature in x = cos(theta) on the hemisphere 0 < x < 1, or on -1 < x < 1 without mirror symmetry
  Vector w(nodes);
  if (not mirror) {
    Legendre::GaussLegendre(x, w);
  } else if (particle.shape != ParticleShape::cylinder) {
    Legendre::PositiveGaussLegendre(x, w);
  } else {
    const Numeric eps = particle.aspect_ratio;
    // Gauss-Legendre rules on the end cap [cos(theta_e), 1] and the side [0, cos(theta_e)]
    const Numeric edge = std::cos(std::atan(eps));
    const Index   cap  = nodes / 2;
    Vector        t(cap), wt(cap);
    Legendre::GaussLegendre(t, wt);
    for (Index i = 0; i < cap; ++i) {
      x[i] = edge + 0.5 * (1.0 - edge) * (t[i] + 1.0);
      w[i] = 0.5 * (1.0 - edge) * wt[i];
    }
    t.resize(nodes - cap);
    wt.resize(nodes - cap);
    Legendre::GaussLegendre(t, wt);
    for (Index i = cap; i < nodes; ++i) {
      x[i] = 0.5 * edge * (t[i - cap] + 1.0);
      w[i] = 0.5 * edge * wt[i - cap];
    }
  }

  // The surface r(theta) and rho = r'/r at the nodes (Particle shapes); s is the sin(theta) of the Schmidt
  // recurrences behind d, which carry it as a factor, so that d / s is exact to rounding
  Vector r2(nodes);
  for (Index i = 0; i < nodes; ++i) {
    s[i]                    = std::sin(std::acos(x[i]));
    std::tie(r2[i], rho[i]) = surface(particle, x[i]);
  }

  // The Bessel functions at k r and eta k r, started Delta n degrees above nmax (Special functions)
  Numeric xi_max = 0.0;
  for (Index i = 0; i < nodes; ++i) {
    W[i]   = w[i] * r2[i];
    xi[i]  = k * std::sqrt(r2[i]);
    xi_max = std::max(xi_max, xi[i]);
  }
  const Numeric nm      = static_cast<Numeric>(nmax);
  const Numeric t       = std::max(xi_max * std::abs(refractive_index), nm);
  const auto    extra_x = static_cast<Index>(1.2 * std::sqrt(std::max(xi_max, nm)) + 3.0);
  const auto    extra_z = static_cast<Index>(t + 4.0 * std::cbrt(t) + 1.2 * std::sqrt(t)) - nmax + 5;
  Vector        ratio_x(nmax + extra_x + 1);
  ComplexVector ratio_z(nmax + extra_z + 1);
  for (Index i = 0; i < nodes; ++i) {
    spherical_bessel_j(xi[i], nmax, ratio_x, j[i], dj[i]);
    spherical_bessel_y(xi[i], nmax, y[i], dy[i]);
    spherical_bessel_j(refractive_index * xi[i], nmax, ratio_z, J[i], dJ[i]);
    wigner_d(std::acos(x[i]), d[i], dd[i]);
  }
}

namespace {
/** The surface integrals and T-matrices of the orders of one surface (The T-matrix of one azimuthal order)
 *
 * The storage is allocated once for the surface.  Order m uses the degrees
 * n0 = max(m, 1), ..., nmax: the rows n - 1 of the per-degree arrays and the
 * leading (2N)^2 elements of the Q storage, N = nmax - n0 + 1.
 */
class Orders {
  const Surface& surface;
  Index          nodes;
  // Per node: W rho, 1 / xi; then, per order, m^2 / s^2, W m / s and W m rho / s
  Vector Wr, ixi, q, Wo, Wor;
  // Per degree n and node, contiguous in the node: the Bessel functions, J / (eta xi), and per order d and tau
  Matrix         j, dj, d, dp;
  ComplexMatrix  h, dh, J, dJ, Jz;
  ComplexVector  storage;
  solve_workdata pivots;

 public:
  explicit Orders(const Surface& surface_)
      : surface(surface_),
        nodes(surface.x.size()),
        Wr(nodes),
        ixi(nodes),
        q(nodes),
        Wo(nodes),
        Wor(nodes),
        j(surface.nmax, nodes),
        dj(surface.nmax, nodes),
        d(surface.nmax, nodes),
        dp(surface.nmax, nodes),
        h(surface.nmax, nodes),
        dh(surface.nmax, nodes),
        J(surface.nmax, nodes),
        dJ(surface.nmax, nodes),
        Jz(surface.nmax, nodes),
        storage(4 * surface.nmax * surface.nmax),
        pivots(2 * surface.nmax) {
    for (Index i = 0; i < nodes; ++i) {
      Wr[i]            = surface.W[i] * surface.rho[i];
      ixi[i]           = 1.0 / surface.xi[i];
      const Complex iz = 1.0 / (surface.refractive_index * surface.xi[i]);
      for (Index a = 0; a < surface.nmax; ++a) {
        const Index n = a + 1;
        j[a, i]       = surface.j[i, n];
        dj[a, i]      = surface.dj[i, n];
        h[a, i]       = Complex{surface.j[i, n], surface.y[i, n]};
        dh[a, i]      = Complex{surface.dj[i, n], surface.dy[i, n]};
        J[a, i]       = surface.J[i, n];
        dJ[a, i]      = surface.dJ[i, n];
        Jz[a, i]      = surface.J[i, n] * iz;
      }
    }
  }

  //! The size 2N of the matrices of order m
  [[nodiscard]] Index size(Index m) const { return 2 * (surface.nmax - std::max<Index>(m, 1) + 1); }

  //! Q and Rg Q of order m into [2N, 2N]
  void q_matrices(ComplexMatrixView Q, ComplexMatrixView RgQ, Index m) {
    const Index   nmax   = surface.nmax;
    const Index   n0     = std::max<Index>(m, 1);
    const Index   N      = nmax - n0 + 1;
    const Complex eta    = surface.refractive_index;
    const Numeric order  = static_cast<Numeric>(m);
    const bool    mirror = surface.mirror;

    for (Index i = 0; i < nodes; ++i) {
      const Numeric si = surface.s[i];
      q[i]             = order * order / (si * si);
      Wo[i]            = surface.W[i] * order / si;
      Wor[i]           = Wo[i] * surface.rho[i];
      for (Index n = n0; n <= nmax; ++n) {
        d[n - 1, i]  = surface.d[i, n, m];
        dp[n - 1, i] = surface.dd[i, n, m];
      }
    }

    // The surface integrals K^{pq} and Rg K^{pq}, and the blocks of Q and Rg Q from them
    constexpr Complex imaginary{0.0, 1.0};
    for (Index n = n0; n <= nmax; ++n) {
      const Index   a  = n - 1;
      const Numeric an = static_cast<Numeric>(n * (n + 1));
      for (Index nn = n0; nn <= nmax; ++nn) {
        const Index   b   = nn - 1;
        const Numeric ann = static_cast<Numeric>(nn * (nn + 1));
        // 2 a_{nn'} on a hemisphere, a_{nn'} on the sphere
        const Numeric a2 = (mirror ? 1.0 : 0.5) * std::sqrt(static_cast<Numeric>(2 * n + 1) / an) *
                           std::sqrt(static_cast<Numeric>(2 * nn + 1) / ann);
        // Mirror symmetry leaves K^{12} and K^{21} for n + n' even and K^{11} and K^{22} for n + n' odd
        const bool even  = (n + nn) % 2 == 0;
        const bool first = not mirror or even, second = m > 0 and (not mirror or not even);

        Complex k11{}, k12{}, k21{}, k22{}, r11{}, r12{}, r21{}, r22{};
        if (first) {
          for (Index i = 0; i < nodes; ++i) {
            const Numeric f1  = surface.W[i] * (q[i] * d[a, i] * d[b, i] + dp[a, i] * dp[b, i]);
            const Numeric f2  = Wr[i] * an * d[a, i] * dp[b, i] * ixi[i];
            const Numeric f3  = Wr[i] * ann * dp[a, i] * d[b, i];
            k12              += J[b, i] * (f1 * dh[a, i] + f2 * h[a, i]);
            r12              += J[b, i] * (f1 * dj[a, i] + f2 * j[a, i]);
            const Complex g   = f1 * dJ[b, i] + f3 * Jz[b, i];
            k21              += g * h[a, i];
            r21              += g * j[a, i];
          }
        }
        if (second) {
          for (Index i = 0; i < nodes; ++i) {
            const Numeric e1  = Wo[i] * (d[a, i] * dp[b, i] + dp[a, i] * d[b, i]);
            const Numeric e2  = Wor[i] * d[a, i] * d[b, i];
            const Complex x   = an * ixi[i] * dJ[b, i];
            const Complex y   = ann * Jz[b, i];
            k11              += e1 * J[b, i] * h[a, i];
            r11              += e1 * J[b, i] * j[a, i];
            k22              += e1 * dJ[b, i] * dh[a, i] + e2 * (x * h[a, i] + y * dh[a, i]);
            r22              += e1 * dJ[b, i] * dj[a, i] + e2 * (x * j[a, i] + y * dj[a, i]);
          }
        }

        const Index   u = n - n0, v = nn - n0;
        const Complex K11 = a2 * k11, K12 = a2 * k12, K21 = a2 * k21, K22 = a2 * k22;
        const Complex R11 = a2 * r11, R12 = a2 * r12, R21 = a2 * r21, R22 = a2 * r22;
        Q[u, v]           = imaginary * (eta * K21 - K12);
        Q[u, v + N]       = -(eta * K11 + K22);
        Q[u + N, v]       = -(eta * K22 + K11);
        Q[u + N, v + N]   = imaginary * (K21 - eta * K12);
        RgQ[u, v]         = imaginary * (eta * R21 - R12);
        RgQ[u, v + N]     = -(eta * R11 + R22);
        RgQ[u + N, v]     = -(eta * R22 + R11);
        RgQ[u + N, v + N] = imaginary * (R21 - eta * R12);
      }
    }
  }

  //! T = -Rg Q Q^{-1} of order m, by the LU solve of -Rg Q = T Q in the storage of T
  ComplexMatrix tmatrix(Index m) {
    const Index       size = this->size(m);
    ComplexMatrix     T(size, size);
    ComplexMatrixView Q{matpack::mdview_t<Complex, 2>{storage.data_handle(), std::array{size, size}}};
    q_matrices(Q, T, m);
    // A Q or Rg Q that is not finite, such as of Bessel functions that overflow inside a strongly absorbing particle,
    // has no T-matrix, and what LAPACK makes of it differs between implementations
    const auto finite = [size](const ComplexMatrixView& M) {
      for (Index a = 0; a < size; ++a)
        for (Index b = 0; b < size; ++b)
          if (not(std::isfinite(M[a, b].real()) and std::isfinite(M[a, b].imag()))) return false;
      return true;
    };
    if (not(finite(Q) and finite(T)))
      arts::user_error(
          "The T-matrix does not converge: its Q matrices of order {} are not finite at {} degrees and {} "
          "surface nodes",
          m,
          surface.nmax,
          nodes);
    T *= -1.0;
    right_solve_inplace(T, Q, pivots);
    return T;
  }
};
}  // namespace

std::pair<ComplexMatrix, ComplexMatrix> order_q(const Surface& surface, Index m) {
  Orders        orders(surface);
  const Index   size = orders.size(m);
  ComplexMatrix Q(size, size), RgQ(size, size);
  orders.q_matrices(Q, RgQ, m);
  return {std::move(Q), std::move(RgQ)};
}

Jones TMatrix::block(Index m, Index n, Index nn) const { return order_block(orders[m], std::max<Index>(m, 1), n, nn); }

TMatrix solve(const Particle& particle, Numeric accuracy, Index nodes_per_degree, bool random) {
  // The relative changes of the sums are rounding below the machine epsilon (Truncation and convergence)
  const Numeric delta = 0.1 * accuracy;
  if (not(delta >= std::numeric_limits<Numeric>::epsilon()))
    arts::user_error(
        "The T-matrix cannot converge to the accuracy {}: its criterion, 0.1 accuracy, is below the "
        "double-precision epsilon {}",
        accuracy,
        std::numeric_limits<Numeric>::epsilon());
  const Numeric xv    = 2.0 * pi * particle.radius / particle.wavelength;
  Index         nmax  = std::max<Index>(4, static_cast<Index>(xv + 4.05 * std::cbrt(xv)));
  Index         nodes = nodes_per_degree * nmax, n1 = 0;

  // The m = 0 sums at nmax and nodes; the convergence fails once they are not finite
  std::optional<Surface> surface;
  ComplexMatrix          T0;
  const auto             order0 = [&] {
    surface.emplace(particle, nmax, nodes);
    T0              = Orders(*surface).tmatrix(0);
    const auto sums = order0_sums(T0, nmax);
    if (not std::isfinite(sums.first) or not std::isfinite(sums.second))
      arts::user_error(
          "The T-matrix does not converge to the accuracy {}: it is not finite at {} degrees and {} "
          "quadrature nodes",
          accuracy,
          nmax,
          nodes);
    return sums;
  };
  // The larger relative change of the two sums
  const auto change = [](std::pair<Numeric, Numeric> now, std::pair<Numeric, Numeric> before) {
    return std::max(std::abs((before.first - now.first) / now.first),
                    std::abs((before.second - now.second) / now.second));
  };

  // Converge n_max on the m = 0 sums (Truncation and convergence).  n_1 is at least n_max / 2 + 1, so random
  // orientation fails as soon as that exceeds the Wigner tables of its expansion.
  std::pair<Numeric, Numeric> sums{0.0, 0.0};
  for (;; ++nmax) {
    if (random and not is_wigner_ready(static_cast<int>(4 * (nmax / 2 + 1) + 1)))
      arts::user_error(
          "The random-orientation T-matrix needs n_1 >= {} at {} degrees, which the Wigner tables do not "
          "cover: initialize larger ones (WignerInit) for 3j symbols with j1 + j2 + j3 up to {}",
          nmax / 2 + 1,
          nmax,
          4 * (nmax / 2 + 1));
    nodes               = nodes_per_degree * nmax;
    const auto previous = std::exchange(sums, order0());
    if (random) n1 = truncation_degree(T0, nmax, sums.first, sums.second, delta);
    if (change(sums, previous) <= delta) break;
  }

  // Then N_G at that n_max, until converged; the convergence fails once an added node no longer reduces the change
  Numeric previous_change = std::numeric_limits<Numeric>::infinity();
  for (;;) {
    ++nodes;
    const auto    previous = std::exchange(sums, order0());
    const Numeric current  = change(sums, previous);
    if (current <= delta) break;
    if (current >= previous_change)
      arts::user_error(
          "The T-matrix does not converge to the accuracy {}: at {} degrees, the relative change of its "
          "sums is {} at {} quadrature nodes and {} at {}",
          accuracy,
          nmax,
          previous_change,
          nodes - 1,
          current,
          nodes);
    previous_change = current;
  }

  if (not random) n1 = nmax;

  TMatrix T{.nmax = nmax, .nodes = nodes, .n1 = n1, .orders = {}};
  T.orders.reserve(n1 + 1);
  T.orders.push_back(std::move(T0));
  Orders orders(*surface);
  for (Index m = 1; m <= n1; ++m) T.orders.push_back(orders.tmatrix(m));
  return T;
}

CrossSections cross_sections(const TMatrix& T, Numeric wavelength, bool random) {
  const Index top  = random ? T.n1 : T.nmax;
  Numeric     qext = 0.0, qsca = 0.0;
  for (Index m = 0; m < static_cast<Index>(T.orders.size()); ++m) {
    const Numeric weight = m == 0 ? 1.0 : 2.0;
    const Index   n0     = std::max<Index>(m, 1);
    for (Index n = n0; n <= T.nmax; ++n) qext += weight * trace(T.block(m, n, n)).real();
    for (Index n = n0; n <= top; ++n)
      for (Index nn = n0; nn <= top; ++nn)
        for (const Complex& t : T.block(m, n, nn).data) qsca += weight * std::norm(t);
  }
  const Numeric scale = wavelength * wavelength / (2.0 * pi);
  return {.extinction = -scale * qext, .scattering = scale * qsca};
}

namespace {
/** pi^m_n and tau^m_n at a particle-frame polar angle, 0 <= m <= n <= nmax (Special functions)
 *
 * The Schmidt recurrences carry sin^m(theta) as a factor and never divide by sin(theta), so m d / sin(theta), with
 * their sin(theta), is accurate to rounding up to the poles.  Only where sin(theta) = 0 is it 0 / 0, and there the
 * functions take their limits.  The storage is reused from angle to angle.
 */
class AngularFunctions {
  Matrix d, dd;

 public:
  Matrix pi, tau;

  explicit AngularFunctions(Index nmax)
      : d(nmax + 1, nmax + 1), dd(nmax + 1, nmax + 1), pi(nmax + 1, nmax + 1), tau(nmax + 1, nmax + 1) {}

  void at(Numeric theta) {
    const Index nmax = pi.nrows() - 1;
    pi               = 0.0;
    tau              = 0.0;
    const Numeric s  = std::sin(theta);
    if (s == 0.0) {
      const bool south = std::cos(theta) < 0.0;
      for (Index n = 1; n <= nmax; ++n) {
        const Numeric c = 0.5 * std::sqrt(static_cast<Numeric>(n * (n + 1)));
        pi[n, 1]        = south ? parity(n + 1) * c : c;
        tau[n, 1]       = south ? parity(n) * c : c;
      }
      return;
    }
    wigner_d(theta, d, dd);
    for (Index n = 1; n <= nmax; ++n)
      for (Index m = 0; m <= n; ++m) {
        pi[n, m]  = static_cast<Numeric>(m) * d[n, m] / s;
        tau[n, m] = dd[n, m];
      }
  }

  //! P^m_n = [[pi, tau], [tau, pi]]
  [[nodiscard]] Jones P(Index m, Index n) const { return {pi[n, m], tau[n, m], tau[n, m], pi[n, m]}; }

  //! |pi| + |tau|, the sum of the moduli of P^m_n's distinct elements
  [[nodiscard]] Numeric modulus(Index m, Index n) const { return std::abs(pi[n, m]) + std::abs(tau[n, m]); }
};

//! The particle frame of the Euler angles alpha and beta [rad]: B, by its rows, takes laboratory components to it
struct ParticleFrame {
  std::array<Vector3, 3> B;

  ParticleFrame(Numeric alpha, Numeric beta) {
    const Numeric ca = std::cos(alpha), sa = std::sin(alpha), cb = std::cos(beta), sb = std::sin(beta);
    B = {Vector3{ca * cb, sa * cb, -sb}, Vector3{-sa, ca, 0.0}, Vector3{ca * sb, sa * sb, cb}};
  }

  [[nodiscard]] Vector3 operator()(const Vector3& v) const { return {dot(B[0], v), dot(B[1], v), dot(B[2], v)}; }
};

//! A laboratory direction in the particle frame: its angles, and R = A_p B A_l, the components of the laboratory
//! basis (theta, phi) in the particle-frame basis (theta_p, phi_p)
struct FrameDirection {
  Numeric theta, phi;
  Jones   R;
};

FrameDirection frame_direction(const ParticleFrame& frame, Numeric theta, Numeric phi) {
  const Numeric ct = std::cos(theta), st = std::sin(theta), cp = std::cos(phi), sp = std::sin(phi);
  // The polar angle from the particle-frame components is accurate to rounding also next to the axis, where the
  // arccos of the z component would turn the rounding of that component into sqrt(eps)
  const Vector3  k = frame({st * cp, st * sp, ct});
  FrameDirection dir;
  dir.theta = std::atan2(std::hypot(k[0], k[1]), k[2]);
  dir.phi   = std::atan2(k[1], k[0]);

  const Numeric ctp = std::cos(dir.theta), stp = std::sin(dir.theta);
  const Numeric cpp = std::cos(dir.phi), spp = std::sin(dir.phi);
  const Vector3 theta_l = frame({ct * cp, ct * sp, -st}), phi_l = frame({-sp, cp, 0.0});
  const Vector3 theta_p{ctp * cpp, ctp * spp, -stp}, phi_p{-spp, cpp, 0.0};
  dir.R = Jones{dot(theta_p, theta_l), dot(theta_p, phi_l), dot(phi_p, theta_l), dot(phi_p, phi_l)};
  return dir;
}

//! c_{nn'} = i^{n'-n-1} sqrt((2n+1)(2n'+1) / (n n' (n+1) (n'+1))) of the amplitude series at [n, n'], 1 <= n, n'
ComplexMatrix amplitude_coefficients(Index nmax) {
  ComplexMatrix c(nmax + 1, nmax + 1, Complex{});
  for (Index n = 1; n <= nmax; ++n)
    for (Index nn = 1; nn <= nmax; ++nn)
      c[n, nn] = i_power(nn - n - 1) * std::sqrt(static_cast<Numeric>((2 * n + 1) * (2 * nn + 1)) /
                                                 static_cast<Numeric>(n * nn * (n + 1) * (nn + 1)));
  return c;
}

/** The T-matrix contracted with the angular functions of one direction (Amplitude matrix in a fixed orientation)
 *
 * Contracted with the incident direction: W_{m n} = sum_{n'} c_{nn'} T_{nn'}(m) P^0_{n'}; with the scattered one:
 * V_{m n'} = sum_n c_{nn'} P_n T_{nn'}(m).  Either leaves the sum over the other direction's degrees, so each further
 * direction costs O(nmax^2) instead of O(nmax^3).  It also keeps the sums of the moduli of the terms, |c| (|pi| +
 * |tau|) times those of T.
 */
class Contraction {
  Index              nmax;
  std::vector<Jones> blocks;
  Vector             moduli;

 public:
  explicit Contraction(Index nmax_)
      : nmax(nmax_), blocks((nmax_ + 1) * nmax_, Jones{Complex{}}), moduli((nmax_ + 1) * nmax_) {}

  //! At order m and the remaining degree n
  [[nodiscard]] Jones&       at(Index m, Index n) { return blocks[m * nmax + n - 1]; }
  [[nodiscard]] const Jones& at(Index m, Index n) const { return blocks[m * nmax + n - 1]; }
  [[nodiscard]] Numeric&     modulus(Index m, Index n) { return moduli[m * nmax + n - 1]; }
  [[nodiscard]] Numeric      modulus(Index m, Index n) const { return moduli[m * nmax + n - 1]; }

  void incident(const TMatrix& T, const ComplexMatrix& c, const AngularFunctions& incident) {
    for (Index m = 0; m <= nmax; ++m)
      for (Index n = std::max<Index>(m, 1); n <= nmax; ++n) {
        Jones   W{Complex{}};
        Numeric sum = 0.0;
        for (Index nn = std::max<Index>(m, 1); nn <= nmax; ++nn) {
          const Jones Tb  = T.block(m, n, nn);
          W              += c[n, nn] * (Tb * incident.P(m, nn));
          Numeric terms   = 0.0;
          for (const Complex& z : Tb.data) terms += std::abs(z);
          sum += std::abs(c[n, nn]) * incident.modulus(m, nn) * terms;
        }
        at(m, n)      = W;
        modulus(m, n) = sum;
      }
  }

  void scattered(const TMatrix& T, const ComplexMatrix& c, const AngularFunctions& scattered) {
    for (Index m = 0; m <= nmax; ++m)
      for (Index nn = std::max<Index>(m, 1); nn <= nmax; ++nn) {
        Jones   V{Complex{}};
        Numeric sum = 0.0;
        for (Index n = std::max<Index>(m, 1); n <= nmax; ++n) {
          const Jones Tb  = T.block(m, n, nn);
          V              += c[n, nn] * (scattered.P(m, n) * Tb);
          Numeric terms   = 0.0;
          for (const Complex& z : Tb.data) terms += std::abs(z);
          sum += std::abs(c[n, nn]) * scattered.modulus(m, n) * terms;
        }
        at(m, nn)      = V;
        modulus(m, nn) = sum;
      }
  }
};

/** The particle-frame amplitude matrix times k, and the sum of the moduli of its terms, from a contraction and the
 * other direction's angular functions: incident_side tells which direction the contraction has.  Order m enters with
 * (2 - delta_m0) cos(m dphi) on its co-polar and (2 - delta_m0) sin(m dphi), negated below the diagonal, on its
 * cross-polar elements, dphi the scattered minus the incident particle-frame azimuth. */
std::pair<Jones, Numeric> particle_amplitude(const Contraction&      contraction,
                                             const AngularFunctions& other,
                                             Numeric                 dphi,
                                             bool                    incident_side) {
  const Index nmax = other.pi.nrows() - 1;
  Jones       Sp{Complex{}};
  Numeric     magnitude = 0.0;
  for (Index m = 0; m <= nmax; ++m) {
    const Numeric w = m == 0 ? 1.0 : 2.0;
    Jones         Sm{Complex{}};
    for (Index n = std::max<Index>(m, 1); n <= nmax; ++n) {
      Sm        += incident_side ? other.P(m, n) * contraction.at(m, n) : contraction.at(m, n) * other.P(m, n);
      magnitude += w * other.modulus(m, n) * contraction.modulus(m, n);
    }
    const Numeric md = static_cast<Numeric>(m) * dphi, fc = w * std::cos(md), fs = w * std::sin(md);
    Sp += Jones{fc * Sm[0, 0], fs * Sm[0, 1], -fs * Sm[1, 0], fc * Sm[1, 1]};
  }
  return {Sp, magnitude};
}
}  // namespace

std::vector<Amplitude> amplitude(const TMatrix& T, Numeric wavelength, std::span<const Geometry> geometries) {
  const Index         nmax = T.nmax;
  const ComplexMatrix c    = amplitude_coefficients(nmax);
  AngularFunctions    incident_functions(nmax), scattered_functions(nmax);
  Contraction         contraction(nmax);
  const Numeric       k = 2.0 * pi / wavelength;

  std::vector<Amplitude> out(geometries.size());
  for (Size row = 0; row < geometries.size(); ++row) {
    // The given angles, also at the singular directions (Amplitude matrix in a fixed orientation)
    const Geometry&      g = geometries[row];
    const ParticleFrame  frame(Conversion::deg2rad(g.alpha), Conversion::deg2rad(g.beta));
    const FrameDirection incident =
        frame_direction(frame, Conversion::deg2rad(g.theta_incident), Conversion::deg2rad(g.phi_incident));
    const FrameDirection scattered =
        frame_direction(frame, Conversion::deg2rad(g.theta_scattered), Conversion::deg2rad(g.phi_scattered));
    incident_functions.at(incident.theta);
    scattered_functions.at(scattered.theta);
    contraction.incident(T, c, incident_functions);
    const auto [Sp, magnitude] =
        particle_amplitude(contraction, scattered_functions, scattered.phi - incident.phi, true);

    // S = R_s^{-1} S^p R_0
    out[row] = {.S = (1.0 / k) * (inv(scattered.R) * Sp * incident.R), .magnitude = magnitude / k};
  }
  return out;
}

AzimuthallyRandom azimuthally_random(const TMatrix& T,
                                     Numeric        wavelength,
                                     const Vector&  tilts,
                                     const Vector&  weights,
                                     const Vector&  za_inc,
                                     const Vector&  delta_aa,
                                     const Vector&  za_scat) {
  const Index         nmax = T.nmax;
  const ComplexMatrix c    = amplitude_coefficients(nmax);
  const Numeric       k    = 2.0 * pi / wavelength;

  // Exact for the truncated T-matrix: its amplitude matrix has azimuthal orders |m| <= nmax in each direction, so the
  // phase matrix has frequencies up to 4 nmax in alpha (both directions turn with the particle) and the absorption
  // integrand, summed over the incident polarization, is a spherical polynomial of degree 2 nmax in the incident
  // direction: 4 nmax + 1 nodes in alpha, nmax + 1 Gauss-Legendre nodes in cos(theta') and 2 nmax + 1 in phi'.  A
  // vertical axis is invariant under alpha.
  const Index n_theta = nmax + 1, n_phi = 2 * nmax + 1;
  Vector      x(n_theta), wx(n_theta);
  Legendre::GaussLegendre(x, wx);
  Vector cos_theta(n_theta), sin_theta(n_theta);
  for (Index q = 0; q < n_theta; ++q) {
    cos_theta[q] = x[q];
    sin_theta[q] = std::sqrt((1.0 - x[q]) * (1.0 + x[q]));
  }

  AzimuthallyRandom out{.phase      = MuelmatTensor3(za_inc.size(), delta_aa.size(), za_scat.size(), Muelmat{0.0}),
                        .extinction = CompactPlanarPropmatVector(za_inc.size(), rtepack::compact_planar_propmat{}),
                        .absorption = CompactPlanarStokvecVector(za_inc.size(), rtepack::compact_planar_stokvec{})};

  // ARTS's azimuth aa, clockwise from north, is phi = 90 deg - aa of the T-matrix's frame with x east and y north,
  // whose bases are then ARTS's (concept.tmatrix): the incident aa = 0 and the scattered aa = delta_aa
  const Numeric phi_incident = pi / 2;
#pragma omp parallel if (not arts_omp_in_parallel())
  {
    AngularFunctions incident_functions(nmax), scattered_functions(nmax), node_functions(nmax);
    Contraction      from_incident(nmax), from_outgoing(nmax);
#pragma omp for schedule(dynamic)
    for (Size i = 0; i < za_inc.size(); ++i) {
      for (Size t = 0; t < tilts.size(); ++t) {
        const bool    vertical = tilts[t] == 0.0 or tilts[t] == 180.0;
        const Index   n_alpha  = vertical ? 1 : 4 * nmax + 1;
        const Numeric weight   = weights[t] / static_cast<Numeric>(n_alpha);
        for (Index j = 0; j < n_alpha; ++j) {
          const ParticleFrame  frame(2.0 * pi * static_cast<Numeric>(j) / static_cast<Numeric>(n_alpha),
                                     Conversion::deg2rad(tilts[t]));
          const FrameDirection incident = frame_direction(frame, Conversion::deg2rad(za_inc[i]), phi_incident);
          incident_functions.at(incident.theta);
          from_incident.incident(T, c, incident_functions);
          from_outgoing.scattered(T, c, incident_functions);
          const auto lab = [k](const FrameDirection& scattered, const Jones& Sp, const FrameDirection& from) {
            return (1.0 / k) * (inv(scattered.R) * Sp * from.R);
          };

          // Extinction from the forward amplitude: Mishchenko's K11, K12 and K34, the last with ARTS's sign of U
          const Jones forward =
              lab(incident, particle_amplitude(from_incident, incident_functions, 0.0, true).first, incident);
          const rtepack::compact_planar_propmat K{wavelength * (forward[0, 0] + forward[1, 1]).imag(),
                                                  wavelength * (forward[0, 0] - forward[1, 1]).imag(),
                                                  wavelength * (forward[0, 0] - forward[1, 1]).real()};
          out.extinction[i] += weight * K;

          // The phase matrix
          for (Size a = 0; a < delta_aa.size(); ++a)
            for (Size s = 0; s < za_scat.size(); ++s) {
              const FrameDirection scattered = frame_direction(
                  frame, Conversion::deg2rad(za_scat[s]), phi_incident - Conversion::deg2rad(delta_aa[a]));
              scattered_functions.at(scattered.theta);
              const Jones S =
                  lab(scattered,
                      particle_amplitude(from_incident, scattered_functions, scattered.phi - incident.phi, true).first,
                      incident);
              out.phase[i, a, s] += weight * rtepack::to_muelmat(S);
            }

          // The absorption vector of the outgoing direction za_inc: K e_1 minus the integral of Z e_1 over the
          // incident directions
          Numeric scattered_I = 0.0, scattered_Q = 0.0;
          for (Index q = 0; q < n_theta; ++q)
            for (Index r = 0; r < n_phi; ++r) {
              const Numeric        phi  = 2.0 * pi * static_cast<Numeric>(r) / static_cast<Numeric>(n_phi);
              const FrameDirection node = frame_direction(frame, std::atan2(sin_theta[q], cos_theta[q]), phi);
              node_functions.at(node.theta);
              const Jones S =
                  lab(incident,
                      particle_amplitude(from_outgoing, node_functions, incident.phi - node.phi, false).first,
                      node);
              const Muelmat Z  = rtepack::to_muelmat(S);
              const Numeric w  = wx[q] * 2.0 * pi / static_cast<Numeric>(n_phi);
              scattered_I     += w * Z[0, 0];
              scattered_Q     += w * Z[1, 0];
            }
          out.absorption[i] += weight * rtepack::compact_planar_stokvec{K.A() - scattered_I, K.B() - scattered_Q};
        }
      }
    }
  }
  return out;
}

Expansion expansion(const TMatrix& T, Numeric scattering, Numeric wavelength) {
  const Index N = T.n1;

  // T^{(1)} and T^{(2)} of the orders +m and -m: the first column of U t U, with U = [[1, 1], [1, -1]] and t the block
  // of order m, or for -m that block with its cross terms negated, Z t Z with Z = diag(1, -1) (Expansion coefficients
  // in random orientation)
  const auto tk = [&T](Index m, Index n, Index nn) {
    constexpr Jones U{1.0, 1.0, 1.0, -1.0}, Z{1.0, 0.0, 0.0, -1.0};
    const Jones     t = T.block(std::abs(m), n, nn);
    const Jones     H = U * (m < 0 ? Z * t * Z : t) * U;
    return std::pair{H[0, 0], H[1, 0]};
  };

  // B^{(k)}_{n_s m n} at [n_s, m + N, n - 1], n_s = 0..2N, m = -N..N, n = 1..N
  const Index    count = 2 * N + 1;
  ComplexTensor3 B1(count, count, N, Complex{}), B2(count, count, N, Complex{});

  // D^{(k)}_{m nu n} at [m + N, nu - 1, n - 1], m = -N..N, nu, n = 1..N
  Tensor3        D1(count, N, N, 0.0), D2(count, N, N, 0.0), D3(count, N, N, 0.0), D4(count, N, N, 0.0);
  ComplexTensor3 D5(count, N, N, Complex{});

  // The symbols have j1 + j2 + j3 <= 4 n_1, so wigxjpf needs factorials up to 4 n_1 + 1
  if (not is_wigner_ready(static_cast<int>(4 * N + 1)))
    arts::user_error(
        "Initialize the Wigner tables (WignerInit) before random-orientation T-matrix calculations; "
        "n_1 = {} needs 3j symbols with j1 + j2 + j3 up to {}",
        N,
        4 * N);

  // Each thread computes whole terms in the order of their sums, so the results do not depend on the threads.
  // The wigxjpf scratch and the per-thread sums are allocated once per thread, around its loops.
  struct scratch_guard {
    scratch_guard(Index n1) {
      errno = 0;
      arts_wigner_thread_init(static_cast<int>(6 * n1 + 1));
    }
    ~scratch_guard() { arts_wigner_thread_free(); }
  };
  bool failed = false;

#pragma omp parallel if (not arts_omp_in_parallel()) reduction(|| : failed)
  {
    const scratch_guard guard(N);
    ComplexVector       A1(N + 1), A2(N + 1);
#pragma omp for schedule(dynamic)
    for (Index n = 1; n <= N; ++n) {
      for (Index ns = 0; ns <= N + n; ++ns) {
        const Index lo = std::max<Index>(1, std::abs(n - ns)), hi = std::min<Index>(N, n + ns);
        for (Index nu = lo; nu <= hi; ++nu) {
          Complex sum1{}, sum2{};
          for (Index m = 0; m <= std::min(n, nu); ++m) {
            const Numeric c      = clebsch_gordan(n, m, ns, 0, nu, m);
            const auto [t1, t2]  = tk(m, n, nu);
            sum1                += c * t1;
            sum2                += c * t2;
            if (m > 0) {
              const Numeric sigma    = parity(n + ns + nu);
              const auto [tm1, tm2]  = tk(-m, n, nu);
              sum1                  += c * sigma * tm1;
              sum2                  += c * sigma * tm2;
            }
          }
          const Complex f = i_power(nu - n) / std::sqrt(static_cast<Numeric>(2 * nu + 1));
          A1[nu]          = f * sum1;
          A2[nu]          = f * sum2;
        }

        for (Index m = std::max<Index>(1 - ns, -n); m <= std::min<Index>(ns + 1, n); ++m) {
          Complex sum1{}, sum2{};
          for (Index nu = lo; nu <= hi; ++nu) {
            const Numeric c  = clebsch_gordan(n, m, ns, 1 - m, nu, 1);
            sum1            += c * A1[nu];
            sum2            += c * A2[nu];
          }
          B1[ns, m + N, n - 1] = sum1;
          B2[ns, m + N, n - 1] = sum2;
        }
      }
    }
    failed = failed or errno == EDOM;
  }

#pragma omp parallel for if (not arts_omp_in_parallel()) schedule(dynamic)
  for (Index n = 1; n <= N; ++n) {
    for (Index nu = 1; nu <= N; ++nu) {
      const Index top = N + std::min(n, nu);
      for (Index m = -std::min(n, nu); m <= std::min(n, nu); ++m) {
        Numeric d1 = 0.0, d2 = 0.0;
        for (Index ns = std::abs(m - 1); ns <= top; ++ns) {
          const Numeric dn  = static_cast<Numeric>(2 * ns + 1);
          d1               += dn * (B1[ns, m + N, n - 1] * std::conj(B1[ns, m + N, nu - 1])).real();
          d2               += dn * (B2[ns, m + N, n - 1] * std::conj(B2[ns, m + N, nu - 1])).real();
        }
        D1[m + N, nu - 1, n - 1] = d1;
        D2[m + N, nu - 1, n - 1] = d2;
      }
      for (Index m = std::max<Index>(-n, 2 - nu); m <= std::min<Index>(n, nu + 2); ++m) {
        Numeric d3 = 0.0, d4 = 0.0;
        Complex d5{};
        for (Index ns = std::abs(m - 1); ns <= top; ++ns) {
          const Numeric dn  = static_cast<Numeric>(2 * ns + 1);
          d3               += dn * (B1[ns, m + N, n - 1] * std::conj(B1[ns, 2 - m + N, nu - 1])).real();
          d4               += dn * (B2[ns, m + N, n - 1] * std::conj(B2[ns, 2 - m + N, nu - 1])).real();
          d5               += dn * B2[ns, m + N, n - 1] * std::conj(B1[ns, 2 - m + N, nu - 1]);
        }
        D3[m + N, nu - 1, n - 1] = d3;
        D4[m + N, nu - 1, n - 1] = d4;
        D5[m + N, nu - 1, n - 1] = d5;
      }
    }
  }

  // The g sums of degree l into g[l] = [g1, g2, g3, g4] and g5[l];
  // C(n -m; l 0 | nu -m) = (-1)^(n + l + nu) C(n m; l 0 | nu m), as in the original
  const Index   lmax = 2 * N;
  Matrix        g(lmax + 1, 4, 0.0);
  ComplexVector g5(lmax + 1, Complex{});
  const auto    sums = [&](Index l, VectorView cg0, VectorView cg2) {
    Numeric g1 = 0.0, g2 = 0.0, g3 = 0.0, g4 = 0.0;
    Complex g5l{};
    for (Index n = 1; n <= N; ++n) {
      for (Index nu = std::max<Index>(1, std::abs(n - l)); nu <= std::min<Index>(N, n + l); ++nu) {
        const Numeric sigma = parity(n + l + nu);
        const Numeric f     = std::sqrt(static_cast<Numeric>(2 * n + 1) / static_cast<Numeric>(2 * nu + 1));
        const Index   mmax  = std::min(n, nu);
        for (Index m = 0; m <= mmax; ++m) cg0[m] = clebsch_gordan(n, m, l, 0, nu, m);
        const Numeric c1 = f * cg0[1];

        Numeric m1 = 0.0, m2 = 0.0;
        for (Index m = -mmax; m <= mmax; ++m) {
          const Numeric c  = m < 0 ? sigma * cg0[-m] : cg0[m];
          m1              += c * D1[m + N, nu - 1, n - 1];
          m2              += c * D2[m + N, nu - 1, n - 1];
        }
        g1 += c1 * m1;
        g2 += c1 * sigma * m2;
        if (l < 2) continue;

        // C(n -m; l 2 | nu 2 - m) at m + n
        const Index lo = std::max<Index>(-n, 2 - nu), hi = std::min<Index>(n, nu + 2);
        for (Index m = lo; m <= hi; ++m) cg2[m + n] = clebsch_gordan(n, -m, l, 2, nu, 2 - m);
        Numeric m3 = 0.0, m4 = 0.0;
        Complex m5{};
        for (Index m = lo; m <= hi; ++m) {
          const Numeric c  = cg2[m + n];
          m3              += c * D3[m + N, nu - 1, n - 1];
          m4              += c * D4[m + N, nu - 1, n - 1];
          m5              += c * D5[m + N, nu - 1, n - 1];
        }
        const Numeric c2  = f * cg2[1 + n];
        g3               += c2 * m3;
        g4               += c2 * sigma * m4;
        g5l              -= c1 * m5;
      }
    }
    g[l, 0] = g1;
    g[l, 1] = g2;
    g[l, 2] = g3;
    g[l, 3] = g4;
    g5[l]   = g5l;
  };

  // Every degree l = 0, ..., 2 n_1 of the series, the degrees in parallel
#pragma omp parallel if (not arts_omp_in_parallel()) reduction(|| : failed)
  {
    const scratch_guard guard(N);
    Vector              cg0(N + 1), cg2(2 * N + 3);
#pragma omp for schedule(dynamic)
    for (Index l = 0; l <= lmax; ++l) sums(l, cg0, cg2);
    failed = failed or errno == EDOM;
  }
  if (failed) arts::user_error("Cannot compute the Wigner symbols of the random-orientation T-matrix");

  const Numeric dk = wavelength * wavelength / (4.0 * pi * scattering);
  Expansion     out{.alpha1 = Vector(lmax + 1),
                    .alpha2 = Vector(lmax + 1),
                    .alpha3 = Vector(lmax + 1),
                    .alpha4 = Vector(lmax + 1),
                    .beta1  = Vector(lmax + 1),
                    .beta2  = Vector(lmax + 1)};
  for (Index l = 0; l <= lmax; ++l) {
    const Numeric cl = static_cast<Numeric>(2 * l + 1) * dk;
    out.alpha1[l]    = cl * (g[l, 0] + g[l, 1]);
    out.alpha2[l]    = cl * (g[l, 2] + g[l, 3]);
    out.alpha3[l]    = cl * (g[l, 2] - g[l, 3]);
    out.alpha4[l]    = cl * (g[l, 0] - g[l, 1]);
    out.beta1[l]     = 2.0 * cl * g5[l].real();
    out.beta2[l]     = 2.0 * cl * g5[l].imag();
  }
  return out;
}

namespace {
//! The scattering matrix at u = cos(Theta) in ARTS's Stokes convention; P is scratch of l_max + 1 (Scattering matrix from the expansion)
CompactPlanarMuelmat scattering_matrix_at(const Expansion& expansion, Numeric u, VectorView P) {
  const Index lmax = expansion.alpha1.size() - 1;
  Legendre::legendre_polynomials(P, u);

  Numeric f11 = 0.0, f44 = 0.0, f2 = 0.0, f3 = 0.0, f12 = 0.0, f34 = 0.0;
  Numeric p22 = 0.25 * (1.0 + u) * (1.0 + u), p2m2 = 0.25 * (1.0 - u) * (1.0 - u);
  Numeric p02          = 0.25 * std::sqrt(6.0) * (u * u - 1.0);
  Numeric p22_previous = 0.0, p2m2_previous = 0.0, p02_previous = 0.0;
  for (Index l = 0; l <= lmax; ++l) {
    f11 += expansion.alpha1[l] * P[l];
    f44 += expansion.alpha4[l] * P[l];
    if (l < 2) continue;

    f2  += (expansion.alpha2[l] + expansion.alpha3[l]) * p22;
    f3  += (expansion.alpha2[l] - expansion.alpha3[l]) * p2m2;
    f12 += expansion.beta1[l] * p02;
    f34 += expansion.beta2[l] * p02;

    // The generalized spherical functions of degree l + 1
    const Numeric dl      = static_cast<Numeric>(l);
    const Numeric pl1     = static_cast<Numeric>(2 * l + 1);
    const Numeric pl2     = dl * (dl + 1.0) * u;
    const Numeric pl3     = (dl + 1.0) * (dl * dl - 4.0);
    const Numeric pl4     = 1.0 / (dl * ((dl + 1.0) * (dl + 1.0) - 4.0));
    const Numeric next22  = (pl1 * (pl2 - 4.0) * p22 - pl3 * p22_previous) * pl4;
    const Numeric next2m2 = (pl1 * (pl2 + 4.0) * p2m2 - pl3 * p2m2_previous) * pl4;
    const Numeric next02 =
        (pl1 * u * p02 - std::sqrt(dl * dl - 4.0) * p02_previous) / std::sqrt((dl + 1.0) * (dl + 1.0) - 4.0);
    p22_previous  = std::exchange(p22, next22);
    p2m2_previous = std::exchange(p2m2, next2m2);
    p02_previous  = std::exchange(p02, next02);
  }

  // In ARTS's Stokes convention, whose U is the negative of Mishchenko's: F34 changes sign
  return {f11, f12, 0.5 * (f2 + f3), 0.5 * (f2 - f3), -f34, f44};
}
}  // namespace

CompactPlanarMuelmatVector scattering_matrix(const Expansion& expansion, Index angles) {
  CompactPlanarMuelmatVector F(angles);
  Vector                     P(expansion.alpha1.size());
  for (Index j = 0; j < angles; ++j)
    F[j] =
        scattering_matrix_at(expansion, std::cos(pi * static_cast<Numeric>(j) / static_cast<Numeric>(angles - 1)), P);
  return F;
}

CompactPlanarMuelmatVector legendre_series(const Expansion& expansion) {
  // a_l = 2 pi int F Y_l0 du, Y_l0 = sqrt((2 l + 1) / 4 pi) P_l(u), by the l_max + 1 node Gauss-Legendre rule
  const Index                n = expansion.alpha1.size();
  Vector                     u(n), w(n), F_P(n), P(n);
  CompactPlanarMuelmatVector a(n);
  Legendre::GaussLegendre(u, w);
  for (Index i = 0; i < n; ++i) {
    const CompactPlanarMuelmat F = scattering_matrix_at(expansion, u[i], F_P);
    Legendre::legendre_polynomials(P, u[i]);
    for (Index l = 0; l < n; ++l)
      a[l] += (2.0 * pi * w[i] * std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * pi)) * P[l]) * F;
  }
  return a;
}

SizeQuadrature size_quadrature(const SizeDistribution& sd) {
  Numeric r1 = sd.lower_ratio * sd.a, r2 = sd.upper_ratio * sd.a;
  using enum TMatrixSizeDistribution;
  if (sd.kind == power_law) {
    // The power law: (r2 - r1) / ln(r2 / r1) = a with r2 = 2 a (1 + b) - r1, for r1 in (0, a].  The left side is
    // -a in the limit r1 -> 0, where the logarithm diverges, and positive at r1 = a for b > 0.
    // ln(y / x) as a difference of logarithms where y / x overflows, which would otherwise end the logarithm's
    // divergence at the overflow and give a false sign change there
    const auto f = [&sd](Numeric x) {
      const Numeric y = 2.0 * sd.a * (1.0 + sd.b) - x, q = y / x;
      return (y - x) / (std::isfinite(q) ? std::log(q) : std::log(y) - std::log(x)) - sd.a;
    };
    const Numeric lo = 0.0, hi = sd.a, flo = -sd.a, fhi = f(hi);
    if (not(fhi > 0.0))
      arts::user_error(
          "The power law has no radii for the effective variance b = {}: (r_2 - r_1) / ln(r_2 / r_1) "
          "= a with r_2 = 2 a (1 + b) - r_1 has no root r_1 in (0, a] that double precision resolves",
          sd.b);
    // Boost's TOMS 748 at least halves the bracket in each pass of at most four evaluations, after two initial ones,
    // so the halvings from the largest double to the subnormal spacing bound the evaluations to any resolvable root
    using limits                        = std::numeric_limits<Numeric>;
    constexpr std::uintmax_t bound      = 2 + 4 * (limits::max_exponent - limits::min_exponent + limits::digits);
    std::uintmax_t           iterations = bound;
    const auto [r1lo, r1hi]             = boost::math::tools::toms748_solve(
        f, lo, hi, flo, fhi, boost::math::tools::eps_tolerance<Numeric>(), iterations);
    if (iterations >= bound)
      arts::user_error(
          "The power-law radius r_1 for the effective variance b = {} is not resolvable in double "
          "precision: it is in [{}, {}] after the {} evaluations that resolve any root above the subnormal "
          "range",
          sd.b,
          r1lo,
          r1hi,
          bound);
    r1 = 0.5 * (r1lo + r1hi);
    r2 = 2.0 * sd.a * (1.0 + sd.b) - r1;
  }

  // Gauss-Legendre on [r1, r2], or on [0, r1] and [r1, r2] for the modified power law
  Vector t(sd.nodes), wt(sd.nodes);
  Legendre::GaussLegendre(t, wt);
  const auto interval = [&](Numeric lo, Numeric hi, Index offset, SizeQuadrature& q) {
    for (Index i = 0; i < sd.nodes; ++i) {
      q.radius[offset + i] = 0.5 * (hi - lo) * t[i] + 0.5 * (hi + lo);
      q.weight[offset + i] = 0.5 * (hi - lo) * wt[i];
    }
  };
  SizeQuadrature q;
  q.lower = r1;
  q.upper = r2;
  if (sd.kind == modified_power_law) {
    q.radius.resize(2 * sd.nodes);
    q.weight.resize(2 * sd.nodes);
    interval(0.0, r1, 0, q);
    interval(r1, r2, sd.nodes, q);
  } else {
    q.radius.resize(sd.nodes);
    q.weight.resize(sd.nodes);
    interval(r1, r2, 0, q);
  }

  const Index count = q.radius.size();
  for (Index i = 0; i < count; ++i) {
    const Numeric r = q.radius[i];
    switch (sd.kind) {
      case modified_gamma:
        q.weight[i] *= std::pow(r, sd.a) * std::exp(-sd.a / sd.gamma * std::pow(r / sd.b, sd.gamma));
        break;
      case lognormal:          q.weight[i] *= std::exp(-Math::pow2(std::log(r / sd.a)) * 0.5 / sd.b) / r; break;
      case power_law:          q.weight[i] /= r * r * r; break;
      case gamma:              q.weight[i] *= std::pow(r, (1.0 - 3.0 * sd.b) / sd.b) * std::exp(-r / (sd.a * sd.b)); break;
      case modified_power_law: q.weight[i] *= r > r1 ? std::pow(r / r1, sd.b) : 1.0; break;
    }
  }
  Numeric total = 0.0;
  for (Index i = 0; i < count; ++i) total += q.weight[i];
  q.weight /= total;

  Numeric g = 0.0, reff = 0.0;
  for (Index i = 0; i < count; ++i) {
    g    += q.radius[i] * q.radius[i] * q.weight[i];
    reff += q.radius[i] * q.radius[i] * q.radius[i] * q.weight[i];
  }
  reff         /= g;
  Numeric veff  = 0.0;
  for (Index i = 0; i < count; ++i)
    veff += (q.radius[i] - reff) * (q.radius[i] - reff) * q.radius[i] * q.radius[i] * q.weight[i];
  q.effective_radius   = reff;
  q.effective_variance = veff / (g * reff * reff);
  return q;
}
namespace {
//! The series of an expansion
constexpr std::array expansion_series{&Expansion::alpha1,
                                      &Expansion::alpha2,
                                      &Expansion::alpha3,
                                      &Expansion::alpha4,
                                      &Expansion::beta1,
                                      &Expansion::beta2};
}  // namespace

void Expansion::add(const Expansion& other, Numeric weight) {
  const Size length = other.size();
  for (const auto series : expansion_series) {
    Vector& mine = this->*series;
    if (length > mine.size()) {
      Vector longer(length, 0.0);
      longer[Range(0, mine.size())] = mine;
      mine                          = std::move(longer);
    }
    const Vector& theirs = other.*series;
    for (Size l = 0; l < length; ++l) mine[l] += theirs[l] * weight;
  }
}

Expansion& Expansion::operator/=(Numeric x) {
  for (const auto series : expansion_series) this->*series /= x;
  return *this;
}

Average average(Particle                particle,
                Numeric                 radius_factor,
                const SizeDistribution& distribution,
                Numeric                 accuracy,
                Index                   nodes_per_degree) {
  const auto quadrature = size_quadrature(distribution);

  // The cross sections and the expansion coefficients weighted by w_i C_sca,i, largest radius first: its
  // series is the longest, unless a smaller particle needs more degrees
  Average out{.cross_sections     = {.extinction = 0.0, .scattering = 0.0},
              .expansion          = {},
              .effective_radius   = quadrature.effective_radius,
              .effective_variance = quadrature.effective_variance};
  for (Index i = static_cast<Index>(quadrature.radius.size()) - 1; i >= 0; --i) {
    particle.radius         = radius_factor * quadrature.radius[i];
    const auto    T         = solve(particle, accuracy, nodes_per_degree, true);
    const auto    cross     = cross_sections(T, particle.wavelength, true);
    const auto    expansion = solver::expansion(T, cross.scattering, particle.wavelength);
    const Numeric weight    = quadrature.weight[i] * cross.scattering;

    out.expansion.add(expansion, weight);
    out.cross_sections.scattering += weight;
    out.cross_sections.extinction += cross.extinction * quadrature.weight[i];
  }

  out.expansion /= out.cross_sections.scattering;
  return out;
}
}  // namespace tmatrix::solver
