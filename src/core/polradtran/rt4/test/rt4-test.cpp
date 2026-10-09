// Tests of the RT4 wrapper against closed forms derived here, independently
// of RT4's own formulation (it uses a doubling-adding scheme with a
// B-and-slope source; the references below integrate the transfer equation
// along each stream analytically).
#include <arts_constants.h>
#include <physics_funcs.h>
#include <rt4.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rt4 = polradtran::rt4;

namespace {
constexpr Numeric pi        = Constant::pi;
constexpr Numeric frequency = 50e9;

using rt4::down;
using rt4::up;

Index size(const Vector& v) { return static_cast<Index>(v.size()); }

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

//! P(n + 1, x) = 1 - exp(-x) sum_{m <= n} x^m / m!, by its positive series
Numeric lower_gamma_ratio(int n, Numeric x) {
  Numeric t = 1.0;
  for (int m = 1; m <= n; m++) t *= x / m;
  Numeric sum = 0.0;
  for (int m = n + 1;; m++) {
    t   *= x / m;
    sum += t;
    if (t < 1e-17 * sum) break;
  }
  return std::exp(-x) * sum;
}

struct iq {
  Numeric I{}, Q{};
};

/** Exact solution of dV/ds = -K V + a B(s) over a path of length L, with
 *  K = [[k, 0], [kappa, k]], a = [aI, aQ] and B linear from Bs (start) to
 *  Be (end).  Using exp(-K u) = exp(-k u) (1 - kappa u N), N = [[0, 0],
 *  [1, 0]], and E_n = int_0^L u^n exp(-k u) du = n! / k^(n+1) P(n+1, k L):
 *    I = I0 e^{-kL} + aI (Be E0 - g E1)
 *    Q = e^{-kL} (Q0 - kappa L I0) + aQ Be E0 - (aQ g + kappa aI Be) E1
 *        + kappa aI g E2,                       g = (Be - Bs) / L.
 */
iq propagate(iq v0, Numeric L, Numeric k, Numeric kappa, Numeric aI, Numeric aQ, Numeric Bs, Numeric Be) {
  const Numeric x  = k * L;
  const Numeric E0 = lower_gamma_ratio(0, x) / k;
  const Numeric E1 = lower_gamma_ratio(1, x) / (k * k);
  const Numeric E2 = 2.0 * lower_gamma_ratio(2, x) / (k * k * k);
  const Numeric g  = (Be - Bs) / L;
  const Numeric ex = std::exp(-x);
  return {.I = v0.I * ex + aI * (Be * E0 - g * E1),
          .Q = ex * (v0.Q - kappa * L * v0.I) + aQ * Be * E0 - (aQ * g + kappa * aI * Be) * E1 + kappa * aI * g * E2};
}

//! Per-stream, per-layer optics of a non-scattering-equivalent medium:
//! effective extinction [[k, 0], [kappa, k]] per length and emission a per length.
struct stream_optics {
  Numeric k{}, kappa{}, aI{}, aQ{};
};

//! Three-layer atmosphere, top-down, temperatures varying.
struct atmosphere {
  Vector  height{3000.0, 2000.0, 1000.0, 0.0};
  Vector  temperature{220.0, 240.0, 265.0, 285.0};
  Vector  gas{1e-4, 3e-4, 5e-4};  // tau = 0.1, 0.3, 0.5
  Numeric sky     = Constant::cosmic_microwave_background_temperature;
  Numeric surface = 290.0;

  Index   nlay() const { return size(height) - 1; }
  Numeric dz(Index l) const { return std::abs(height[l] - height[l + 1]); }
};

rt4::problem base_problem(const atmosphere& atm, Index nstokes) {
  rt4::problem p;
  p.nstokes             = nstokes;
  p.nmu                 = 8;
  p.quad                = polradtran::quadrature_type::double_gauss;
  p.extra_mu            = Vector{1.0};
  p.frequency           = frequency;
  p.height              = atm.height;
  p.temperature         = atm.temperature;
  p.gas_extinction      = atm.gas;
  p.layer_optics_index  = ArrayOfIndex(atm.nlay(), -1);
  p.sky_temperature     = atm.sky;
  p.surface_temperature = atm.surface;
  p.ground              = polradtran::lambertian_surface{.albedo = 0.0};
  return p;
}

//! Downwelling closed form at every level for each stream, from an
//! unpolarized sky.  optics(l, i) gives the stream optics of layer l.
std::vector<std::vector<iq>> closed_form_down(const atmosphere&                                 atm,
                                              const Vector&                                     mu,
                                              const std::function<stream_optics(Index, Index)>& optics) {
  std::vector<std::vector<iq>> down(atm.nlay() + 1, std::vector<iq>(mu.size()));
  for (Index i = 0; i < size(mu); i++) {
    down[0][i] = {.I = planck(frequency, atm.sky), .Q = 0.0};
    for (Index l = 0; l < atm.nlay(); l++) {
      const auto o   = optics(l, i);
      down[l + 1][i] = propagate(down[l][i],
                                 atm.dz(l) / mu[i],
                                 o.k,
                                 o.kappa,
                                 o.aI,
                                 o.aQ,
                                 planck(frequency, atm.temperature[l]),
                                 planck(frequency, atm.temperature[l + 1]));
    }
  }
  return down;
}

//! Upwelling closed form at every level for each stream, given the
//! upwelling radiance just above the surface.
std::vector<std::vector<iq>> closed_form_up(const atmosphere&                                 atm,
                                            const Vector&                                     mu,
                                            const std::vector<iq>&                            surface_up,
                                            const std::function<stream_optics(Index, Index)>& optics) {
  std::vector<std::vector<iq>> up(atm.nlay() + 1, std::vector<iq>(mu.size()));
  for (Index i = 0; i < size(mu); i++) {
    up[atm.nlay()][i] = surface_up[i];
    for (Index l = atm.nlay() - 1; l >= 0; l--) {
      const auto o = optics(l, i);
      up[l][i]     = propagate(up[l + 1][i],
                               atm.dz(l) / mu[i],
                               o.k,
                               o.kappa,
                               o.aI,
                               o.aQ,
                               planck(frequency, atm.temperature[l + 1]),
                               planck(frequency, atm.temperature[l]));
    }
  }
  return up;
}

//! max over levels, streams of |I - I_ref| / I_ref and |Q - Q_ref| / I_ref
struct deviation {
  Numeric I{}, Q{};
  Numeric max() const { return std::max(I, Q); }
};

deviation compare(const Tensor3& rt4, const std::vector<std::vector<iq>>& ref) {
  deviation d;
  for (Index l = 0; l < rt4.extent(0); l++) {
    for (Index i = 0; i < rt4.extent(1); i++) {
      const auto& r = ref[l][i];
      d.I           = std::max(d.I, std::abs(rt4[l, i, 0] - r.I) / r.I);
      if (rt4.extent(2) > 1) d.Q = std::max(d.Q, std::abs(rt4[l, i, 1] - r.Q) / r.I);
    }
  }
  return d;
}

deviation combine(deviation a, deviation b) { return {std::max(a.I, b.I), std::max(a.Q, b.Q)}; }

std::vector<iq> black_surface(const atmosphere& atm, Index nmu) {
  return std::vector<iq>(nmu, iq{.I = planck(frequency, atm.surface), .Q = 0.0});
}

void report(std::string_view name, deviation d, Numeric tol) {
  std::cout << std::format("{:<72} max rel dev I {:9.3e}, Q {:9.3e} (tolerance {:.1e})\n", name, d.I, d.Q, tol);
  require(d.max() <= tol, std::format("{}: deviation {:.3e} exceeds {:.1e}", name, d.max(), tol));
}

/** (g) RT4's quadratures against their defining exactness:
 *  D: nmu-point Gauss on [0, 1], exact for mu^k, k <= 2 nmu - 1;
 *  G: half of the 2 nmu-point Gauss rule on [-1, 1], exact for even
 *     mu^(2k), 2k <= 4 nmu - 1;
 *  L: half of the 2 nmu-point Lobatto rule, exact for even mu^(2k),
 *     2k <= 4 nmu - 3, and with a node at mu = 1. */
void test_quadrature() {
  Numeric worst = 0.0;
  for (Index n : {1, 2, 5, 8, 16}) {
    const auto D = polradtran::get_quadrature(n, polradtran::quadrature_type::double_gauss);
    const auto G = polradtran::get_quadrature(n, polradtran::quadrature_type::gauss);
    const auto L = polradtran::get_quadrature(n, polradtran::quadrature_type::lobatto);
    for (const auto* q : {&D, &G, &L}) {
      require(size(q->mu) == n and size(q->weights) == n, "quadrature size");
      for (Index i = 0; i < n; i++)
        require(q->mu[i] > 0 and q->mu[i] <= 1 and (i == 0 or q->mu[i] > q->mu[i - 1]),
                "quadrature nodes must be ascending in (0, 1]");
    }
    const auto moment = [](const polradtran::quadrature& q, Index k) {
      Numeric s = 0.0;
      for (Index i = 0; i < size(q.mu); i++) s += q.weights[i] * std::pow(q.mu[i], k);
      return s;
    };
    for (Index k = 0; k <= 2 * n - 1; k++)
      worst = std::max(worst, std::abs(moment(D, k) - 1.0 / static_cast<Numeric>(k + 1)));
    for (Index k = 0; 2 * k <= 4 * n - 1; k++)
      worst = std::max(worst, std::abs(moment(G, 2 * k) - 1.0 / static_cast<Numeric>(2 * k + 1)));
    for (Index k = 0; 2 * k <= 4 * n - 3; k++)
      worst = std::max(worst, std::abs(moment(L, 2 * k) - 1.0 / static_cast<Numeric>(2 * k + 1)));
    require(L.mu[n - 1] == 1.0, "Lobatto must include mu = 1");
  }
  std::cout << std::format("{:<72} max moment error {:9.3e}\n", "(g) quadrature exactness D/G/L", worst);
  require(worst < 1e-13, "quadrature moments are not exact");
}

/** (e) Isothermal blackbody: every radiance equals B(T).  Compares RT4's
 *  Planck function and its per-micrometre to per-Hz conversion with ARTS'
 *  planck(). */
void test_planck() {
  Numeric worst = 0.0;
  for (Numeric f : {10e9, 89e9, 664e9}) {
    for (Numeric T : {Constant::cosmic_microwave_background_temperature, 150.0, 300.0}) {
      atmosphere atm;
      atm.temperature = Vector{T, T, T, T};
      atm.sky         = T;
      atm.surface     = T;
      auto p          = base_problem(atm, 2);
      p.frequency     = f;
      const auto r    = rt4::solve(p);
      const auto B    = planck(f, T);
      for (const auto* t : {&r.up, &r.down}) {
        for (Index l = 0; l < t->extent(0); l++) {
          for (Index i = 0; i < t->extent(1); i++) {
            worst = std::max(worst, std::abs((*t)[l, i, 0] - B) / B);
            require((*t)[l, i, 1] == 0.0, "isothermal blackbody Q must be 0");
          }
        }
      }
    }
  }
  report("(e) isothermal blackbody == ARTS planck, 3 f x 3 T", {worst, 0.0}, 1e-12);
}

/** (a) Gas-only atmosphere over a black surface, nstokes 1 and 2. */
void test_gas_only() {
  const atmosphere atm;
  const auto       gas = [&](Index l, Index) {
    return stream_optics{.k = atm.gas[l], .kappa = 0.0, .aI = atm.gas[l], .aQ = 0.0};
  };
  for (Index ns : {1, 2}) {
    const auto r    = rt4::solve(base_problem(atm, ns));
    const auto down = closed_form_down(atm, r.mu, gas);
    const auto up   = closed_form_up(atm, r.mu, black_surface(atm, r.mu.size()), gas);
    require(r.mu.size() == 9 and r.mu[8] == 1.0 and r.weights[8] == 0.0, "extra angle mu = 1 with weight 0");
    report(std::format("(a) gas-only, black surface, nstokes {}", ns),
           combine(compare(r.up, up), compare(r.down, down)),
           1e-12);
  }
}

/** (b) The same atmosphere as "scattering" layers with zero phase matrix,
 *  K = k 1 and a = [k, 0]: exercises the doubling, whose initial layer is
 *  first order in the sublayer thickness. */
void test_doubling_convergence() {
  const atmosphere atm;
  const auto       gas = [&](Index l, Index) {
    return stream_optics{.k = atm.gas[l], .kappa = 0.0, .aI = atm.gas[l], .aQ = 0.0};
  };
  std::vector<Numeric> err;
  Numeric              mu_min = 0.0;
  for (Numeric mdt : {1e-5, 1e-6, 1e-7}) {
    auto        p   = base_problem(atm, 2);
    const Index nmu = p.nmu + static_cast<Index>(p.extra_mu.size());
    for (Index l = 0; l < atm.nlay(); l++) {
      rt4::layer_optics o{.extinction = Tensor4(2, nmu, 2, 2, 0.0),
                          .absorption = Tensor3(2, nmu, 2, 0.0),
                          .phase      = Tensor6(2, 2, nmu, nmu, 2, 2, 0.0)};
      for (Index h = 0; h < 2; h++) {
        for (Index i = 0; i < nmu; i++) {
          o.extinction[h, i, 0, 0] = o.extinction[h, i, 1, 1] = atm.gas[l];
          o.absorption[h, i, 0]                               = atm.gas[l];
        }
      }
      p.optics.push_back(o);
      p.layer_optics_index[l] = l;
    }
    p.gas_extinction = Vector(atm.nlay(), 0.0);
    p.max_delta_tau  = mdt;
    const auto r     = rt4::solve(p);
    mu_min           = r.mu[0];
    const auto down  = closed_form_down(atm, r.mu, gas);
    const auto up    = closed_form_up(atm, r.mu, black_surface(atm, r.mu.size()), gas);
    const auto d     = combine(compare(r.up, up), compare(r.down, down));
    err.push_back(d.max());
    // The first-order initial layer has a slant thickness of at most
    // max_delta_tau / mu_min, which bounds the relative error.
    report(std::format("(b) doubling, zero phase, max_delta_tau {:.0e}", mdt), d, mdt / mu_min);
  }
  // The number of doublings is int(log2(tau / max_delta_tau)) + 1, so each
  // decade changes the sublayer thickness by a factor 8 or 16.
  for (std::size_t i = 0; i + 1 < err.size(); i++) {
    const Numeric ratio = err[i] / err[i + 1];
    std::cout << std::format("    error ratio between successive max_delta_tau: {:.2f}\n", ratio);
    require(ratio > 4.0 and ratio < 25.0, "doubling error is not first order in max_delta_tau");
  }
}

/** (b2) Stream and Stokes layout.  Scattering layers with direction
 *  dependent, lower-triangular K(mu) = [[k(mu), 0], [kappa(mu), k(mu)]],
 *  absorption a(mu) = [aI(mu), aQ(mu)], and a phase matrix that only scatters
 *  forward into the same stream, P(h <- h)(i, j) = c delta_ij / (2 pi w_j),
 *  for I and Q, with aI = k - c so that energy is conserved.  On a
 *  quadrature stream this is a medium with extinction K - c 1; the
 *  zero-weight extra stream gets no in-scattering.  A
 *  permutation of streams, a Stokes transpose of K, or forward/backward
 *  quadrants swapped all change the answer. */
void test_layout() {
  const atmosphere atm;
  auto             p   = base_problem(atm, 2);
  const Index      nmu = p.nmu + static_cast<Index>(p.extra_mu.size());
  const auto       qw  = polradtran::get_quadrature(p.nmu, p.quad);
  Vector           mu(nmu);
  for (Index i = 0; i < p.nmu; i++) mu[i] = qw.mu[i];
  mu[p.nmu] = p.extra_mu[0];

  const auto k     = [&](Index l, Index i) { return atm.gas[l] * (1.0 + 0.4 * mu[i]); };
  const auto kappa = [&](Index l, Index i) { return 0.3 * atm.gas[l] * mu[i]; };
  const auto c     = [&](Index l, Index i) { return i < p.nmu ? 0.35 * atm.gas[l] : 0.0; };
  // Energy conservation (rt4::solve checks it): the absorption is the extinction minus what the stream scatters
  const auto aI = [&](Index l, Index i) { return k(l, i) - c(l, i); };
  const auto aQ = [&](Index l, Index i) { return -0.2 * atm.gas[l] * (1.0 - mu[i]); };

  for (Index l = 0; l < atm.nlay(); l++) {
    rt4::layer_optics o{.extinction = Tensor4(2, nmu, 2, 2, 0.0),
                        .absorption = Tensor3(2, nmu, 2, 0.0),
                        .phase      = Tensor6(2, 2, nmu, nmu, 2, 2, 0.0)};
    for (Index h = 0; h < 2; h++) {
      for (Index i = 0; i < nmu; i++) {
        o.extinction[h, i, 0, 0] = o.extinction[h, i, 1, 1] = k(l, i);
        o.extinction[h, i, 1, 0]                            = kappa(l, i);
        o.absorption[h, i, 0]                               = aI(l, i);
        o.absorption[h, i, 1]                               = aQ(l, i);
        if (i < p.nmu) o.phase[h, h, i, i, 0, 0] = o.phase[h, h, i, i, 1, 1] = c(l, i) / (2 * pi * qw.weights[i]);
      }
    }
    p.optics.push_back(o);
    p.layer_optics_index[l] = l;
  }
  p.gas_extinction = Vector(atm.nlay(), 0.0);
  p.max_delta_tau  = 1e-7;

  const auto optics = [&](Index l, Index i) {
    return stream_optics{.k = k(l, i) - c(l, i), .kappa = kappa(l, i), .aI = aI(l, i), .aQ = aQ(l, i)};
  };
  const auto r    = rt4::solve(p);
  const auto down = closed_form_down(atm, r.mu, optics);
  const auto up   = closed_form_up(atm, r.mu, black_surface(atm, r.mu.size()), optics);
  report("(b2) layout: K(mu) lower-triangular, forward-only phase",
         combine(compare(r.up, up), compare(r.down, down)),
         p.max_delta_tau / r.mu[0]);
}

//! Fresnel reflectivities from Snell's law, n1 = 1 above
std::pair<Numeric, Numeric> fresnel_vh(Complex n, Numeric mu) {
  const Complex sin2 = 1.0 - mu * mu;
  const Complex cost = std::sqrt(1.0 - sin2 / (n * n));
  const Complex rv   = (n * mu - cost) / (n * mu + cost);
  const Complex rh   = (mu - n * cost) / (mu + n * cost);
  return {std::norm(rv), std::norm(rh)};
}

/** (c) Gas-only atmosphere over specular and discrete surfaces.  The
 *  downwelling at the surface is unpolarized (closed form as in (a)); the
 *  surface returns e + R [I_d, 0] with e = (1 - R) [B_s, 0]; above it the
 *  gas only attenuates Q. */
void test_surfaces() {
  const atmosphere atm;
  const auto       gas = [&](Index l, Index) {
    return stream_optics{.k = atm.gas[l], .kappa = 0.0, .aI = atm.gas[l], .aQ = 0.0};
  };
  const Numeric Bs = planck(frequency, atm.surface);

  for (Complex n : {Complex{1.5, 0.0}, Complex{3.0, 0.2}}) {
    auto p             = base_problem(atm, 2);
    p.ground           = polradtran::fresnel_surface{.refractive_index = n};
    const auto      r  = rt4::solve(p);
    const auto      dn = closed_form_down(atm, r.mu, gas);
    std::vector<iq> sfc(r.mu.size());
    for (Index i = 0; i < size(r.mu); i++) {
      const auto [Rv, Rh] = fresnel_vh(n, r.mu[i]);
      const Numeric Id    = dn[atm.nlay()][i].I;
      // In I = Iv + Ih, Q = Iv - Ih: each polarization is (1 - R_p) B_s / 2 + R_p I_d / 2
      sfc[i] = {.I = 0.5 * ((1 - Rv) + (1 - Rh)) * Bs + 0.5 * (Rv + Rh) * Id,
                .Q = 0.5 * ((1 - Rv) - (1 - Rh)) * Bs + 0.5 * (Rv - Rh) * Id};
    }
    // Sign: the surface is warmer than the sky, and oblique emission is
    // vertically polarized, Iv > Ih.
    for (Index i = 0; i < size(r.mu); i++)
      if (r.mu[i] < 0.9)
        require(r.up[atm.nlay(), i, 1] > 0.0 and sfc[i].Q > 0.0, "Fresnel surface emission must have Q > 0");
    const auto upr = closed_form_up(atm, r.mu, sfc, gas);
    report(std::format("(c) Fresnel surface n = {}{:+}i", n.real(), n.imag()),
           combine(compare(r.up, upr), compare(r.down, dn)),
           1e-12);
  }

  {
    // Non-symmetric R(out, in) detects a transpose: R(Q, I) != R(I, Q).
    Matrix R(2, 2);
    R[0, 0]            = 0.3;
    R[0, 1]            = 0.05;
    R[1, 0]            = -0.1;
    R[1, 1]            = 0.25;
    auto p             = base_problem(atm, 2);
    p.ground           = rt4::specular_surface{.reflectivity = R};
    const auto      r  = rt4::solve(p);
    const auto      dn = closed_form_down(atm, r.mu, gas);
    std::vector<iq> sfc(r.mu.size());
    for (Index i = 0; i < size(r.mu); i++) {
      const Numeric Id = dn[atm.nlay()][i].I;
      sfc[i]           = {.I = (1 - R[0, 0]) * Bs + R[0, 0] * Id, .Q = -R[1, 0] * Bs + R[1, 0] * Id};
    }
    const auto upr = closed_form_up(atm, r.mu, sfc, gas);
    report("(c) specular surface, R = [[0.3, 0.05], [-0.1, 0.25]]",
           combine(compare(r.up, upr), compare(r.down, dn)),
           1e-12);
  }

  for (bool discrete : {false, true}) {
    // Lambertian albedo A as RT4 'L' and as a discrete_surface with the
    // operator 2 A mu_j w_j (I to I) and emission (1 - A) B_s in
    // W m-2 Hz-1 sr-1.  The reference applies that operator to the closed
    // form downwelling at the quadrature nodes.
    constexpr Numeric A   = 0.3;
    auto              p   = base_problem(atm, 2);
    const Index       nmu = p.nmu + static_cast<Index>(p.extra_mu.size());
    const auto        qw  = polradtran::get_quadrature(p.nmu, p.quad);
    if (discrete) {
      rt4::discrete_surface s{.reflection = Tensor4(nmu, nmu, 2, 2, 0.0), .emission = Matrix(nmu, 2, 0.0)};
      for (Index i = 0; i < nmu; i++) {
        s.emission[i, 0] = (1 - A) * Bs;
        for (Index j = 0; j < p.nmu; j++) s.reflection[i, j, 0, 0] = 2 * A * qw.mu[j] * qw.weights[j];
      }
      p.ground = s;
    } else {
      p.ground = polradtran::lambertian_surface{.albedo = A};
    }
    const auto r    = rt4::solve(p);
    const auto dn   = closed_form_down(atm, r.mu, gas);
    Numeric    flux = 0.0;
    for (Index j = 0; j < p.nmu; j++) flux += qw.weights[j] * qw.mu[j] * dn[atm.nlay()][j].I;
    std::vector<iq> sfc(r.mu.size(), iq{.I = (1 - A) * Bs + 2 * A * flux, .Q = 0.0});
    const auto      upr = closed_form_up(atm, r.mu, sfc, gas);
    report(std::format("(c) Lambertian A = 0.3 as {}", discrete ? "discrete_surface" : "lambertian_surface"),
           combine(compare(r.up, upr), compare(r.down, dn)),
           1e-12);
  }
}

/** (d) Isothermal Kirchhoff: sky, layers and surface all at T, a scattering
 *  layer with gas, a gas-only layer, a Fresnel surface.  The exact discrete
 *  solution is I = B, Q = 0 when the emission vector is the Kirchhoff one,
 *    a_s(mu_j) = K_sI - 2 pi sum_{i, h} w_i Z_sI(mu_j <- mu_i, h).
 *
 *  reciprocal: Rayleigh, Z = sigma P / (4 pi) with the m = 0 forms
 *    P_II = 3/8 (3 - mu^2 - mu'^2 + 3 mu^2 mu'^2), P_IQ = 3/8 (1 - 3 mu^2)(1 - mu'^2),
 *    P_QI = 3/8 (1 - mu^2)(1 - 3 mu'^2),          P_QQ = 9/8 (1 - mu^2)(1 - mu'^2)
 *  (mu outgoing, mu' incident) and a = [a1, 0] with a1 from discrete energy
 *  conservation, a1(mu_j) = K - 2 pi sum_i w_i [Z11(up <- h) + Z11(down <- h)](i, j),
 *  which equals the Kirchhoff a1 as P_II is symmetric; a_Q = 0 is Kirchhoff
 *  because the quadrature integrates 1 - 3 mu'^2 to 0.  Q = 0 thus detects
 *  a Stokes or a stream transpose of Z, but not both at once.
 *
 *  non-reciprocal: the Q row of P (P_QI, P_QQ) multiplied by
 *  (1 + 0.3 mu_out) and a from the Kirchhoff formula, so a_Q != 0.  This
 *  detects the full transpose too.  P_II stays reciprocal, so its row and
 *  column sums agree and the medium conserves energy (rt4::solve checks it):
 *  a medium whose P_II were non-reciprocal could not obey both Kirchhoff's
 *  law and energy conservation with one absorption. */
void test_kirchhoff(bool reciprocal) {
  constexpr Numeric T = 260.0;
  atmosphere        atm;
  atm.height      = Vector{2000.0, 1000.0, 0.0};
  atm.temperature = Vector{T, T, T};
  atm.gas         = Vector{5e-5, 2e-4};
  atm.sky         = T;
  atm.surface     = T;

  auto p               = base_problem(atm, 2);
  p.ground             = polradtran::fresnel_surface{.refractive_index = Complex{3.0, 0.2}};
  p.layer_optics_index = ArrayOfIndex{0, -1};
  p.max_delta_tau      = 1e-7;
  const Index nmu      = p.nmu + static_cast<Index>(p.extra_mu.size());
  const auto  qw       = polradtran::get_quadrature(p.nmu, p.quad);
  Vector      mu(nmu), w(nmu, 0.0);
  for (Index i = 0; i < p.nmu; i++) {
    mu[i] = qw.mu[i];
    w[i]  = qw.weights[i];
  }
  mu[p.nmu] = p.extra_mu[0];

  constexpr Numeric sigma = 6e-4, kabs = 4e-4;
  const auto        P = [reciprocal](Index s, Index t, Numeric m, Numeric mp) {
    const Numeric a = m * m, b = mp * mp, f = reciprocal or s == 0 ? 1.0 : 1.0 + 0.3 * m;
    if (s == 0 and t == 0) return f * 3.0 / 8.0 * (3 - a - b + 3 * a * b);
    if (s == 0 and t == 1) return f * 3.0 / 8.0 * (1 - 3 * a) * (1 - b);
    if (s == 1 and t == 0) return f * 3.0 / 8.0 * (1 - a) * (1 - 3 * b);
    return f * 9.0 / 8.0 * (1 - a) * (1 - b);
  };
  rt4::layer_optics o{.extinction = Tensor4(2, nmu, 2, 2, 0.0),
                      .absorption = Tensor3(2, nmu, 2, 0.0),
                      .phase      = Tensor6(2, 2, nmu, nmu, 2, 2, 0.0)};
  for (Index ho = 0; ho < 2; ho++)
    for (Index hi = 0; hi < 2; hi++)
      for (Index i = 0; i < nmu; i++)
        for (Index j = 0; j < nmu; j++)
          for (Index s = 0; s < 2; s++)
            for (Index t = 0; t < 2; t++) o.phase[ho, hi, i, j, s, t] = sigma * P(s, t, mu[i], mu[j]) / (4 * pi);
  for (Index h = 0; h < 2; h++) {
    for (Index j = 0; j < nmu; j++) {
      o.extinction[h, j, 0, 0] = o.extinction[h, j, 1, 1] = sigma + kabs;
      if (reciprocal) {
        Numeric scat = 0.0;
        for (Index i = 0; i < nmu; i++) scat += w[i] * (o.phase[up, h, i, j, 0, 0] + o.phase[down, h, i, j, 0, 0]);
        o.absorption[h, j, 0] = sigma + kabs - 2 * pi * scat;
      } else {
        for (Index s = 0; s < 2; s++) {
          Numeric scat = 0.0;
          for (Index i = 0; i < nmu; i++) scat += w[i] * (o.phase[h, up, j, i, s, 0] + o.phase[h, down, j, i, s, 0]);
          o.absorption[h, j, s] = (s == 0 ? sigma + kabs : 0.0) - 2 * pi * scat;
        }
      }
    }
  }
  p.optics.push_back(o);

  const auto    r = rt4::solve(p);
  const Numeric B = planck(frequency, T);
  deviation     d;
  for (const auto* t : {&r.up, &r.down}) {
    for (Index l = 0; l < t->extent(0); l++) {
      for (Index i = 0; i < t->extent(1); i++) {
        d.I = std::max(d.I, std::abs((*t)[l, i, 0] - B) / B);
        d.Q = std::max(d.Q, std::abs((*t)[l, i, 1]) / B);
      }
    }
  }
  // The doubling is exact for this fixed point, so only round-off remains;
  // it grows with the number of sublayers 2^n, n = int(log2(tau / max_delta_tau)) + 1.
  const Numeric tau     = (sigma + kabs + atm.gas[0]) * atm.dz(0);
  const Numeric n_doubl = std::floor(std::log2(tau / p.max_delta_tau)) + 1;
  report(std::format("(d) isothermal Kirchhoff, {} + gas + Fresnel 3+0.2i",
                     reciprocal ? "Rayleigh" : "non-reciprocal Rayleigh"),
         d,
         std::exp2(n_doubl) * std::numeric_limits<Numeric>::epsilon());
}

void expect_throw(std::string_view what, const std::function<void()>& f) {
  bool threw = false;
  try {
    f();
  } catch (const std::exception& e) {
    threw = true;
    std::cout << std::format("    {:<44} throws: {}\n", what, std::string_view{e.what()}.substr(0, 110));
  }
  require(threw, std::format("{} did not throw", what));
}

/** (f) Error paths; each must throw before the Fortran code (which would
 *  STOP the process) is called. */
void test_errors() {
  const atmosphere atm;
  const auto       good = [&] {
    auto              p   = base_problem(atm, 2);
    const Index       nmu = p.nmu + static_cast<Index>(p.extra_mu.size());
    const auto        q   = polradtran::get_quadrature(p.nmu, p.quad);
    rt4::layer_optics o{.extinction = Tensor4(2, nmu, 2, 2, 0.0),
                        .absorption = Tensor3(2, nmu, 2, 0.0),
                        .phase      = Tensor6(2, 2, nmu, nmu, 2, 2, 0.0)};
    // Energy conservation: a1 = K11 - 2 pi sum_i w_i (1e-6 + 2e-6) over the quadrature streams
    const Numeric scattered = 2 * pi * 3e-6 * std::accumulate(q.weights.begin(), q.weights.end(), 0.0);
    for (Index h = 0; h < 2; h++)
      for (Index i = 0; i < nmu; i++) {
        o.extinction[h, i, 0, 0] = o.extinction[h, i, 1, 1] = 1e-4;
        o.absorption[h, i, 0]                               = 1e-4 - scattered;
        for (Index j = 0; j < nmu; j++) {
          o.phase[h, h, i, j, 0, 0]     = 1e-6;
          o.phase[h, 1 - h, i, j, 0, 0] = 2e-6;
        }
      }
    p.optics.push_back(o);
    p.layer_optics_index[1] = 0;
    return p;
  };
  rt4::solve(good());

  expect_throw("optics that do not conserve energy", [&] {
    auto p                           = good();
    p.optics[0].absorption[0, 0, 0] *= 1.01;
    p.optics[0].absorption[1, 0, 0] *= 1.01;  // mirror symmetric, so only the energy balance fails
    rt4::solve(p);
  });
  {
    // An infinite tolerance accepts anything: the same optics, and an all-zero optics set (K11 = 0)
    auto p                           = good();
    p.optics[0].absorption[0, 0, 0] *= 1.01;
    p.optics[0].absorption[1, 0, 0] *= 1.01;
    auto zero                        = p.optics[0];
    zero.extinction                  = 0.0;
    zero.absorption                  = 0.0;
    zero.phase                       = 0.0;
    p.optics.push_back(zero);
    p.layer_optics_index[0]   = 1;
    p.normalisation_tolerance = std::numeric_limits<Numeric>::infinity();
    rt4::solve(p);
  }
  expect_throw("nstokes = 3", [&] {
    auto p    = good();
    p.nstokes = 3;
    rt4::solve(p);
  });
  expect_throw("nstokes * nmu_total = 2 * 40 > 64", [&] {
    auto p = base_problem(atm, 2);
    p.nmu  = 39;
    rt4::solve(p);
  });
  expect_throw("nstokes * nmu_total = 2 * (32 + 1) > 64", [&] {
    auto p = base_problem(atm, 2);
    p.nmu  = 32;
    rt4::solve(p);
  });
  expect_throw("401 layers", [&] {
    auto p               = base_problem(atm, 1);
    p.nmu                = 1;
    p.extra_mu           = Vector{};
    p.height             = Vector(402, 0.0);
    p.temperature        = Vector(402, 250.0);
    p.gas_extinction     = Vector(401, 0.0);
    p.layer_optics_index = ArrayOfIndex(401, -1);
    for (Index i = 0; i < 402; i++) p.height[i] = static_cast<Numeric>(402 - i);
    rt4::solve(p);
  });
  expect_throw("(nlay + 1) * 64^2 > 301 * 4096", [&] {
    auto p               = base_problem(atm, 2);
    p.nmu                = 31;
    p.height             = Vector(302, 0.0);
    p.temperature        = Vector(302, 250.0);
    p.gas_extinction     = Vector(301, 0.0);
    p.layer_optics_index = ArrayOfIndex(301, -1);
    for (Index i = 0; i < 302; i++) p.height[i] = static_cast<Numeric>(302 - i);
    rt4::solve(p);
  });
  expect_throw("asymmetric extinction", [&] {
    auto p                               = good();
    p.optics[0].extinction[up, 3, 0, 0] *= 1.01;
    rt4::solve(p);
  });
  expect_throw("asymmetric absorption", [&] {
    auto p                              = good();
    p.optics[0].absorption[down, 2, 0] *= 1.01;
    rt4::solve(p);
  });
  expect_throw("asymmetric phase [down, up] vs [up, down]", [&] {
    auto p                                   = good();
    p.optics[0].phase[down, up, 1, 2, 0, 0] *= 1.01;
    rt4::solve(p);
  });
  expect_throw("asymmetric phase [down, down] vs [up, up]", [&] {
    auto p                                     = good();
    p.optics[0].phase[down, down, 1, 2, 0, 0] *= 1.01;
    rt4::solve(p);
  });
  expect_throw("temperature has nlay values", [&] {
    auto p        = good();
    p.temperature = Vector{220.0, 240.0, 265.0};
    rt4::solve(p);
  });
  expect_throw("optics without the extra angle", [&] {
    auto p                 = good();
    p.optics[0].extinction = Tensor4(2, 8, 2, 2, 0.0);
    rt4::solve(p);
  });
  expect_throw("layer_optics_index out of range", [&] {
    auto p                  = good();
    p.layer_optics_index[2] = 1;
    rt4::solve(p);
  });
  expect_throw("discrete_surface of wrong shape", [&] {
    auto p   = good();
    p.ground = rt4::discrete_surface{.reflection = Tensor4(9, 9, 1, 1, 0.0), .emission = Matrix(9, 2, 0.0)};
    rt4::solve(p);
  });
  expect_throw("extra_mu = 0", [&] {
    auto p     = good();
    p.extra_mu = Vector{0.0};
    rt4::solve(p);
  });
  expect_throw("max_delta_tau = 0", [&] {
    auto p          = good();
    p.max_delta_tau = 0.0;
    rt4::solve(p);
  });
  expect_throw("max_delta_tau < 0", [&] {
    auto p          = good();
    p.max_delta_tau = -1e-6;
    rt4::solve(p);
  });
  expect_throw("temperature = 0 K", [&] {
    auto p           = good();
    p.temperature[1] = 0.0;
    rt4::solve(p);
  });
}
}  // namespace

int main() try {
  if (not rt4::available()) throw std::runtime_error("rt4-test needs ENABLE_RT4=ON");
  test_quadrature();
  test_planck();
  test_gas_only();
  test_doubling_convergence();
  test_layout();
  test_surfaces();
  test_kirchhoff(true);
  test_kirchhoff(false);
  test_errors();
  std::cout << "All RT4 tests passed\n";
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << "rt4-test failed: " << e.what() << '\n';
  return EXIT_FAILURE;
}
