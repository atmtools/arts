/* VDISORT against Evans' RT4 polarized doubling-adding solver (src/core/rt4).

   Where the expected answer comes from.  RT4 builds layer reflection,
   transmission and source operators by doubling a thin initial layer and
   combines them by adding.  That is a framework external to VDISORT, which
   eigen-decomposes the discrete-ordinate system of each layer and matches
   boundary and continuity conditions.  Both solvers get the same discrete
   problem:

   - RT4's double-Gauss rule with nmu nodes per hemisphere, which is also
     VDISORT's rule with NQuad = 2 nmu streams (asserted);
   - the same azimuthally averaged phase matrix on those streams;
   - the same scalar extinction and absorption per layer;
   - the same Planck function, linear in optical depth within a layer;
   - the same sky and surface.

   With identical angular discretisation, both solve the same linear system
   of ODEs in optical depth.  VDISORT solves it to round-off.  RT4 solves it
   to round-off in gas-only layers, which it integrates analytically.  In
   every other layer, RT4's error is that of its first-order initial layer,
   so it is first order in that layer's thickness delta <= max_delta_tau.
   The difference must therefore vanish linearly with max_delta_tau (C11),
   and a Richardson extrapolation of RT4 in max_delta_tau must agree with
   VDISORT to second order and round-off.  A convention mismatch (stream
   order, hemisphere, Stokes or stream transpose, normalisation or source)
   would not converge.

   Mapping, layer l top-down, with thickness dz, gas extinction kg, particle
   extinction kp, scattering coefficient sigma <= kp, and Zbar the m = 0
   azimuthal mean of the phase matrix per unit length and steradian:

   - RT4: extinction kp 1, absorption [kp - sigma, 0], phase the [I, Q] block
     of Zbar, gas_extinction kg.
   - VDISORT:
     - tau is the cumulative (kg + kp) dz and omega = sigma / (kg + kp);
     - phase_matrix[cosine, 0, l, i, j] holds 4 pi Zbar / sigma in the [I, Q]
       block, and phase_matrix[sine, 0, l, i, j] the [U, V] block;
     - stream i < N is RT4 (up, mu_i) and stream N + i is RT4 (down, mu_i);
     - the source is c0 + c1 tau in global tau;
     - NFourier = 1.

   The tolerances are derived above direct_tolerance() below. */
#include <arts_constants.h>
#include <physics_funcs.h>
#include <rt4.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "lab-frame.h"

namespace {
constexpr Numeric pi        = Constant::pi;
constexpr Numeric frequency = 89e9;

using rt4::down;
using rt4::up;

Index ssize(const auto& v) { return static_cast<Index>(v.size()); }

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

//! The signed direction cosine of RT4 stream mu in hemisphere h (VDISORT convention: > 0 upward)
Numeric signed_mu(Index h, Numeric mu) { return h == up ? mu : -mu; }

//! RT4's double-Gauss nodes (ascending), followed by the zero-weight extra angles
struct streams {
  Index  nmu{};
  Vector mu{};
  Vector w{};

  Index total() const { return ssize(mu); }
};

streams make_streams(Index nmu, const Vector& extra_mu = {}) {
  const auto  q = rt4::get_quadrature(nmu, rt4::quadrature_type::double_gauss);
  const Index n = nmu + ssize(extra_mu);
  streams     s{.nmu = nmu, .mu = Vector(n, 0.0), .w = Vector(n, 0.0)};
  for (Index i = 0; i < nmu; i++) {
    s.mu[i] = q.mu[i];
    s.w[i]  = q.weights[i];
  }
  for (Index e = 0; e < ssize(extra_mu); e++) s.mu[nmu + e] = extra_mu[e];
  return s;
}

/* The m = 0 azimuthal mean (1 / 2 pi) int Z dDelta-phi of the [I, Q, U, V]
   phase matrix on the streams of both hemispheres, indexed [h_out, h_in,
   mu_out, mu_in] with rt4::down / rt4::up.  It is per unit length and
   steradian, so it includes the scattering coefficient. */
using stream_phase = rtepack::muelmat_tensor4;

//! A normalised ((1 / 4 pi) int P11 dOmega = 1) m = 0 phase matrix P(mu_out <- mu_in), signed cosines
using phase_function = std::function<rtepack::muelmat(Numeric, Numeric)>;

stream_phase on_streams(const streams& s, Numeric sigma, const phase_function& P) {
  const Index  n = s.total();
  stream_phase Z(2, 2, n, n, rtepack::muelmat{0.0});
  for (Index ho = 0; ho < 2; ho++)
    for (Index hi = 0; hi < 2; hi++)
      for (Index io = 0; io < n; io++)
        for (Index ii = 0; ii < n; ii++)
          Z[ho, hi, io, ii] = sigma / (4 * pi) * P(signed_mu(ho, s.mu[io]), signed_mu(hi, s.mu[ii]));
  return Z;
}

/* Rayleigh, m = 0, Q = I_v - I_h, mo outgoing, mi incident (signed):
     P_II = 3/8 (3 - mo^2 - mi^2 + 3 mo^2 mi^2),  P_IQ = 3/8 (1 - 3 mo^2)(1 - mi^2),
     P_QI = 3/8 (1 - mo^2)(1 - 3 mi^2),           P_QQ = 9/8 (1 - mo^2)(1 - mi^2),
     P_UU = 0,  P_VV = 3/2 mo mi.
   The numerical builder below reproduces it to round-off. */
rtepack::muelmat rayleigh(Numeric mo, Numeric mi) {
  const Numeric    a = mo * mo, b = mi * mi;
  rtepack::muelmat P{0.0};
  P[0, 0] = 3.0 / 8.0 * (3 - a - b + 3 * a * b);
  P[0, 1] = 3.0 / 8.0 * (1 - 3 * a) * (1 - b);
  P[1, 0] = 3.0 / 8.0 * (1 - a) * (1 - 3 * b);
  P[1, 1] = 9.0 / 8.0 * (1 - a) * (1 - b);
  P[3, 3] = 1.5 * mo * mi;
  return P;
}

rtepack::muelmat rayleigh_intensity_only(Numeric mo, Numeric mi) {
  rtepack::muelmat P{0.0};
  P[0, 0] = rayleigh(mo, mi)[0, 0];
  return P;
}

/* A constructed, non-reciprocal and Stokes-asymmetric m = 0 phase matrix.
   x = |mo|, y = |mi|, and q = +1 within a hemisphere, -1 between them:
     P_II = 1 + 0.3 (x - 1/2)(1 - y^2) + 0.25 q x (1 - y)^2
     P_IQ = 0.1 (1 - 2 y^2)(1/2 + x)
     P_QI = 0.2 (1 - x^2)(1 + y/2)
     P_QQ = 0.3 (0.2 + x y^2)(1 + 0.2 q)
   It depends on the hemispheres only through q, so it is mirror symmetric,
   as RT4 requires.  The per-incident I normalisation
   (1/2) sum_{both hemispheres} W_o P_II = 1 is exact for any double-Gauss
   rule: sum W (x - 1/2) = 0, and the q terms cancel between the
   hemispheres.  Every transpose changes it:
   - stream: P_II(x, y) != P_II(y, x);
   - Stokes: P_IQ != P_QI;
   - both: P_II and P_QQ are not symmetric;
   - quadrant: q -> -q. */
rtepack::muelmat constructed(Numeric mo, Numeric mi) {
  const Numeric    x = std::abs(mo), y = std::abs(mi), q = mo * mi > 0 ? 1.0 : -1.0;
  rtepack::muelmat P{0.0};
  P[0, 0] = 1 + 0.3 * (x - 0.5) * (1 - y * y) + 0.25 * q * x * (1 - y) * (1 - y);
  P[0, 1] = 0.1 * (1 - 2 * y * y) * (0.5 + x);
  P[1, 0] = 0.2 * (1 - x * x) * (1 + 0.5 * y);
  P[1, 1] = 0.3 * (0.2 + x * y * y) * (1 + 0.2 * q);
  return P;
}

using vdisort_test::lab_frame;
using vdisort_test::tro_elements;
using vdisort_test::tro_matrix;

tro_elements rayleigh_scattering_matrix(Numeric c) {
  return {.F11 = 0.75 * (1 + c * c),
          .F12 = -0.75 * (1 - c * c),
          .F22 = 0.75 * (1 + c * c),
          .F33 = 1.5 * c,
          .F34 = 0.0,
          .F44 = 1.5 * c};
}

//! Henyey-Greenstein g = 0.7 for F11 = F22, with a Rayleigh-like polarisation, F12 = -0.8 F11 sin^2 / (1 + cos^2)
tro_elements polarized_henyey_greenstein(Numeric c) {
  constexpr Numeric g   = 0.7;
  const Numeric     F11 = (1 - g * g) / std::pow(1 + g * g - 2 * g * c, 1.5);
  const Numeric     r   = 1 + c * c;
  return {.F11 = F11,
          .F12 = -0.8 * F11 * (1 - c * c) / r,
          .F22 = F11,
          .F33 = 2 * c * F11 / r,
          .F34 = 0.0,
          .F44 = 2 * c * F11 / r};
}

/* Zbar by the periodic trapezoid rule in Delta-phi, at the nphi midpoints
   (k + 1/2) 2 pi / nphi.  This is spectrally accurate for the smooth
   periodic integrand, and exact for Rayleigh, a trigonometric polynomial of
   degree 2.  The midpoints avoid Delta-phi = 0 and pi, so the forward and
   backward branch of lab_frame (lab-frame.h) is only taken for mu = 1
   pairs. */
stream_phase numeric_azimuthal_mean(const streams& s, Numeric sigma, const tro_matrix& F, Index nphi) {
  const Index  n = s.total();
  stream_phase Z(2, 2, n, n, rtepack::muelmat{0.0});
  for (Index ho = 0; ho < 2; ho++) {
    for (Index hi = 0; hi < 2; hi++) {
      for (Index io = 0; io < n; io++) {
        for (Index ii = 0; ii < n; ii++) {
          rtepack::muelmat sum{0.0};
          for (Index k = 0; k < nphi; k++)
            sum += lab_frame(F,
                             signed_mu(hi, s.mu[ii]),
                             0.0,
                             signed_mu(ho, s.mu[io]),
                             (static_cast<Numeric>(k) + 0.5) * 2 * pi / static_cast<Numeric>(nphi));
          Z[ho, hi, io, ii] = sigma / (4 * pi * static_cast<Numeric>(nphi)) * sum;
        }
      }
    }
  }
  return Z;
}

/* Rescale each incident column (h, j) of Z, all Stokes elements, as ARTS 2
   did for RT4, so that discrete energy conservation holds on the streams:
     2 pi sum_i w_i [Z11(up <- h) + Z11(down <- h)](i, j) = sigma.
   Returns the largest |sigma_j / sigma - 1| before rescaling. */
Numeric renormalize(stream_phase& Z, const streams& s, Numeric sigma) {
  Numeric worst = 0.0;
  for (Index hi = 0; hi < 2; hi++) {
    for (Index j = 0; j < s.nmu; j++) {
      Numeric sj = 0.0;
      for (Index i = 0; i < s.nmu; i++) sj += 2 * pi * s.w[i] * (Z[up, hi, i, j][0, 0] + Z[down, hi, i, j][0, 0]);
      worst = std::max(worst, std::abs(sj / sigma - 1));
      for (Index ho = 0; ho < 2; ho++)
        for (Index i = 0; i < s.total(); i++) Z[ho, hi, i, j] *= sigma / sj;
    }
  }
  return worst;
}

//! Particle optics shared by the layers that use them.  Particle absorption is kp - sigma.
struct optics_set {
  Numeric      kp{};
  Numeric      sigma{};
  stream_phase Z{};
};

/* A non-specular surface given by its m = 0 BRDF rho(mu_out, mu_in) per
   steradian, [I, Q] block, mu >= 0.  This is VDISORT's R^0: VDISORT adds
   pi W_j mu_j R^0_ij.  RT4's discrete operator is pi w_j mu_j rho(mu_i, mu_j). */
struct discrete_reflection {
  std::function<rtepack::muelmat(Numeric, Numeric)> rho;
};

using ground = std::variant<rt4::lambertian_surface, rt4::fresnel_surface, discrete_reflection>;

struct setup {
  Index                   nstokes{2};
  streams                 s{};
  Numeric                 max_delta_tau{1e-7};
  Vector                  height{};        // [nlay + 1], top-down; extinctions are per unit of it
  Vector                  temperature{};   // [nlay + 1] at the interfaces
  Vector                  gas{};           // [nlay]
  ArrayOfIndex            optics_index{};  // [nlay], < 0 for gas-only
  std::vector<optics_set> optics{};
  Numeric                 sky{Constant::cosmic_microwave_background_temperature};
  Numeric                 surface{};
  ground                  g{rt4::lambertian_surface{.albedo = 0.0}};

  Index   nlay() const { return ssize(height) - 1; }
  Numeric dz(Index l) const { return std::abs(height[l] - height[l + 1]); }
  Numeric kp(Index l) const { return optics_index[l] < 0 ? 0.0 : optics[optics_index[l]].kp; }
  Numeric sigma(Index l) const { return optics_index[l] < 0 ? 0.0 : optics[optics_index[l]].sigma; }
};

//! The Kirchhoff emission ([1, 0] - sum_j R(i, j)[:, 0]) B_s of a discrete surface into stream i, [I, Q]
Vector2 discrete_emission(const discrete_reflection& d, const streams& s, Index i, Numeric Bs) {
  Vector2 e{Bs, 0.0};
  for (Index j = 0; j < s.nmu; j++) {
    const auto R  = pi * s.w[j] * s.mu[j] * d.rho(s.mu[i], s.mu[j]);
    e[0]         -= R[0, 0] * Bs;
    e[1]         -= R[1, 0] * Bs;
  }
  return e;
}

rt4::problem rt4_problem(const setup& c) {
  const Index  n = c.s.total(), ns = c.nstokes;
  rt4::problem p;
  p.nstokes  = ns;
  p.nmu      = c.s.nmu;
  p.quad     = rt4::quadrature_type::double_gauss;
  p.extra_mu = Vector(n - c.s.nmu, 0.0);
  for (Index e = 0; e < n - c.s.nmu; e++) p.extra_mu[e] = c.s.mu[c.s.nmu + e];
  p.max_delta_tau       = c.max_delta_tau;
  p.frequency           = frequency;
  p.height              = c.height;
  p.temperature         = c.temperature;
  p.gas_extinction      = c.gas;
  p.layer_optics_index  = c.optics_index;
  p.sky_temperature     = c.sky;
  p.surface_temperature = c.surface;

  for (const auto& o : c.optics) {
    rt4::layer_optics lo{.extinction = Tensor4(2, n, ns, ns, 0.0),
                         .absorption = Tensor3(2, n, ns, 0.0),
                         .phase      = Tensor6(2, 2, n, n, ns, ns, 0.0)};
    for (Index h = 0; h < 2; h++) {
      for (Index i = 0; i < n; i++) {
        for (Index s = 0; s < ns; s++) lo.extinction[h, i, s, s] = o.kp;
        lo.absorption[h, i, 0] = o.kp - o.sigma;
      }
    }
    for (Index ho = 0; ho < 2; ho++)
      for (Index hi = 0; hi < 2; hi++)
        for (Index io = 0; io < n; io++)
          for (Index ii = 0; ii < n; ii++)
            for (Index so = 0; so < ns; so++)
              for (Index si = 0; si < ns; si++) lo.phase[ho, hi, io, ii, so, si] = o.Z[ho, hi, io, ii][so, si];
    p.optics.push_back(std::move(lo));
  }

  if (const auto* d = std::get_if<discrete_reflection>(&c.g)) {
    const Numeric         Bs = planck(frequency, c.surface);
    rt4::discrete_surface ds{.reflection = Tensor4(n, n, ns, ns, 0.0), .emission = Matrix(n, ns, 0.0)};
    for (Index io = 0; io < n; io++) {
      const auto e = discrete_emission(*d, c.s, io, Bs);
      for (Index so = 0; so < ns; so++) ds.emission[io, so] = e[so];
      // The extra-angle columns are applied without weights by RT4, so they stay zero
      for (Index ii = 0; ii < c.s.nmu; ii++) {
        const auto R = pi * c.s.w[ii] * c.s.mu[ii] * d->rho(c.s.mu[io], c.s.mu[ii]);
        for (Index so = 0; so < ns; so++)
          for (Index si = 0; si < ns; si++) ds.reflection[io, ii, so, si] = R[so, si];
      }
    }
    p.ground = ds;
  } else if (const auto* f = std::get_if<rt4::fresnel_surface>(&c.g)) {
    p.ground = *f;
  } else {
    p.ground = std::get<rt4::lambertian_surface>(c.g);
  }
  return p;
}

//! Deliberate VDISORT input errors, to show that the comparison detects them
enum class mistake {
  none,
  stream_transpose,   // P(i <- j) taken from Z(j <- i)
  stokes_transpose,   // P[s, t] taken from Z[t, s]
  full_transpose,     // both: a reciprocal, mirror-symmetric Z is invariant under it
  quadrant_exchange,  // same- and opposite-hemisphere quadrants exchanged
};

std::string_view name(mistake m) {
  switch (m) {
    case mistake::none:              return "none";
    case mistake::stream_transpose:  return "stream transpose";
    case mistake::stokes_transpose:  return "Stokes transpose";
    case mistake::full_transpose:    return "stream and Stokes transpose";
    case mistake::quadrant_exchange: return "same/opposite-hemisphere quadrants exchanged";
  }
  return "";
}

//! The VDISORT phase-matrix block for streams io <- ii, from Z (with a deliberate mistake if asked)
rtepack::muelmat vdisort_block(const stream_phase& Z, Index N, Index io, Index ii, mistake m) {
  Index ho = io < N ? up : down, mo = io % N;
  Index hi = ii < N ? up : down, mi = ii % N;
  if (m == mistake::stream_transpose or m == mistake::full_transpose) {
    std::swap(ho, hi);
    std::swap(mo, mi);
  }
  if (m == mistake::quadrant_exchange) ho = 1 - ho;
  rtepack::muelmat z = Z[ho, hi, mo, mi];
  if (m == mistake::stokes_transpose or m == mistake::full_transpose) matpack::inplace_transpose(z);
  return z;
}

//! The cosine-mode ([I, Q]) and sine-mode ([U, V]) m = 0 blocks of 4 pi Z / sigma
void set_vdisort_phase(
    rtepack::muelmat& cosine, rtepack::muelmat& sine, const rtepack::muelmat& z, Numeric sigma, Index nstokes) {
  const Numeric scale = 4 * pi / sigma;
  for (Index s = 0; s < nstokes; s++)
    for (Index t = 0; t < nstokes; t++) cosine[s, t] = scale * z[s, t];
  for (Index s = 2; s < 4; s++)
    for (Index t = 2; t < 4; t++) sine[s, t] = scale * z[s, t];
}

/* A custom VDISORT surface from rho, optionally deliberately transposed
   (stream or Stokes).  rho is a smooth function, so the callback evaluates
   it at whatever cosines VDISORT passes.  At the nodes this is what matching
   them would give (as fresnel_fourier_modes does), and it also defines the
   reflection into off-node user angles.  The sine (U, V) system is not
   excited, so its callback is zero. */
vdisort::BDRF discrete_bdrf(const discrete_reflection& d, mistake m) {
  const auto cosine =
      [rho = d.rho, m](rtepack::muelmat_matrix_view out, const ConstVectorView& mu_out, const ConstVectorView& mu_in) {
        out = rtepack::muelmat{0.0};
        for (Index i = 0; i < out.nrows(); i++) {
          for (Index j = 0; j < out.ncols(); j++) {
            const Numeric x = std::abs(mu_out[i]), y = std::abs(mu_in[j]);
            const auto    r = m == mistake::stream_transpose ? rho(y, x) : rho(x, y);
            for (Index s = 0; s < 2; s++)
              for (Index t = 0; t < 2; t++) out[i, j][s, t] = m == mistake::stokes_transpose ? r[t, s] : r[s, t];
          }
        }
      };
  const auto zero = [](rtepack::muelmat_matrix_view out, const ConstVectorView&, const ConstVectorView&) {
    out = rtepack::muelmat{0.0};
  };
  return {.cosine      = vdisort::BDRF::func_t{cosine},
          .sine        = vdisort::BDRF::func_t{zero},
          .beam_cosine = vdisort::BDRF::func_t{zero},
          .beam_sine   = vdisort::BDRF::func_t{zero}};
}

vdisort::main_data vdisort_solver(const setup& c,
                                  mistake      phase_mistake   = mistake::none,
                                  mistake      surface_mistake = mistake::none) {
  const Index N = c.s.nmu, NQuad = 2 * N, NL = c.nlay(), ns = c.nstokes;

  Vector  tau(NL), omega(NL);
  Numeric t = 0.0;
  for (Index l = 0; l < NL; l++) {
    const Numeric k  = c.gas[l] + c.kp(l);
    t               += k * c.dz(l);
    tau[l]           = t;
    omega[l]         = c.sigma(l) / k;
  }

  vdisort::phase_matrix_data P(2, 1, NL, NQuad, NQuad, rtepack::muelmat{0.0});
  for (Index l = 0; l < NL; l++) {
    if (c.sigma(l) == 0.0) continue;
    const auto& o = c.optics[c.optics_index[l]];
    for (Index io = 0; io < NQuad; io++)
      for (Index ii = 0; ii < NQuad; ii++)
        set_vdisort_phase(P[vdisort::cosine_mode, 0, l, io, ii],
                          P[vdisort::sine_mode, 0, l, io, ii],
                          vdisort_block(o.Z, N, io, ii, phase_mistake),
                          o.sigma,
                          ns);
  }

  const Numeric            Bs = planck(frequency, c.surface);
  rtepack::stokvec_tensor3 bottom(2, 1, N), top(2, 1, N);
  bottom = rtepack::stokvec{};
  top    = rtepack::stokvec{};
  std::vector<vdisort::BDRF> brdf;
  for (Index i = 0; i < N; i++) {
    top[vdisort::cosine_mode, 0, i] = {planck(frequency, c.sky), 0.0, 0.0, 0.0};
    Vector2 e{};
    if (const auto* d = std::get_if<discrete_reflection>(&c.g)) {
      e = discrete_emission(*d, c.s, i, Bs);
    } else if (const auto* f = std::get_if<rt4::fresnel_surface>(&c.g)) {
      // rtepack::fresnel_reflectance, independent of RT4's Fresnel code
      const auto F = vdisort::brdf::Fresnel{f->refractive_index}(c.s.mu[i]);
      e            = {(1 - F[0, 0]) * Bs, -F[1, 0] * Bs};
    } else {
      e = {(1 - std::get<rt4::lambertian_surface>(c.g).albedo) * Bs, 0.0};
    }
    bottom[vdisort::cosine_mode, 0, i] = {e[0], ns > 1 ? e[1] : 0.0, 0.0, 0.0};
  }
  if (const auto* d = std::get_if<discrete_reflection>(&c.g))
    brdf.push_back(discrete_bdrf(*d, surface_mistake));
  else if (const auto* f = std::get_if<rt4::fresnel_surface>(&c.g))
    brdf = vdisort::brdf::fresnel_fourier_modes(f->refractive_index, 1);
  else
    brdf = vdisort::brdf::lambertian_fourier_modes(std::get<rt4::lambertian_surface>(c.g).albedo, 1);

  // B(tau) = c0 + c1 tau in the global optical depth, linear within each layer
  rtepack::stokvec_matrix source(NL, 2);
  Numeric                 tau_top = 0.0;
  for (Index l = 0; l < NL; l++) {
    const Numeric B0 = planck(frequency, c.temperature[l]), B1 = planck(frequency, c.temperature[l + 1]);
    const Numeric slope = (B1 - B0) / (tau[l] - tau_top);
    source[l, 0]        = {B0 - slope * tau_top, 0.0, 0.0, 0.0};
    source[l, 1]        = {slope, 0.0, 0.0, 0.0};
    tau_top             = tau[l];
  }

  vdisort::main_data v(NQuad,
                       1,
                       AscendingGrid{std::move(tau)},
                       std::move(omega),
                       std::move(P),
                       std::move(bottom),
                       std::move(top),
                       std::move(source),
                       std::move(brdf),
                       0.5,
                       rtepack::stokvec{},
                       0.0);

  Numeric node = 0.0;
  for (Index i = 0; i < N; i++)
    node = std::max({node,
                     std::abs(v.mu()[i] - c.s.mu[i]),
                     std::abs(v.mu()[N + i] + c.s.mu[i]),
                     std::abs(v.weights()[i] - c.s.w[i])});
  require(node < 1e-14,
          std::format("VDISORT's double-Gauss streams and weights must equal RT4's for nmu = {}; they differ by up "
                      "to {:.1e}",
                      N,
                      node));
  return v;
}

//! [level (0 = top), stream, I/Q] for the up- and downward directions
struct field {
  Tensor3 up;
  Tensor3 down;
};

field vdisort_streams(const vdisort::main_data& v, Index nstokes) {
  const Index N = ssize(v.weights()), NL = ssize(v.tau());
  field       f{.up = Tensor3(NL + 1, N, nstokes, 0.0), .down = Tensor3(NL + 1, N, nstokes, 0.0)};
  for (Index l = 0; l <= NL; l++) {
    vdisort::u_data d;
    v.u(d, l == 0 ? 0.0 : v.tau()[l - 1], 0.0);
    for (Index i = 0; i < N; i++) {
      for (Index s = 0; s < nstokes; s++) {
        f.up[l, i, s]   = d.intensities[i][s];
        f.down[l, i, s] = d.intensities[N + i][s];
      }
    }
  }
  return f;
}

//! VDISORT's user-angle formal solution at the setup's extra angles, at every level, both directions
/* VDISORT at the extra angles by its user-angle formal solution, up and down or (upward = false) down only.
   With exact_boundary the boundary radiances where the user rays start are given (the Kirchhoff emission of a
   Lambertian or Fresnel surface, the sky); otherwise VDISORT interpolates them from the streams. */
field vdisort_extra_angles(const vdisort::main_data& v,
                           const setup&              c,
                           bool                      upward         = true,
                           bool                      exact_boundary = false) {
  const Index N = c.s.nmu, NQuad = 2 * N, NL = c.nlay(), ne = c.s.total() - N;
  const Index first = upward ? 0 : ne;  // the user directions are the extra angles up, then down
  Vector      user_mu(2 * ne - first);
  for (Index e = 0; e < ne; e++) {
    if (upward) user_mu[e] = c.s.mu[N + e];
    user_mu[ne + e - first] = -c.s.mu[N + e];
  }
  vdisort::phase_matrix_data user_phase(2, 1, NL, 2 * ne - first, NQuad, rtepack::muelmat{0.0});
  for (Index l = 0; l < NL; l++) {
    if (c.sigma(l) == 0.0) continue;
    const auto& o = c.optics[c.optics_index[l]];
    for (Index u = first; u < 2 * ne; u++) {
      for (Index j = 0; j < NQuad; j++) {
        set_vdisort_phase(user_phase[vdisort::cosine_mode, 0, l, u - first, j],
                          user_phase[vdisort::sine_mode, 0, l, u - first, j],
                          o.Z[u < ne ? up : down, j < N ? up : down, N + u % ne, j % N],
                          o.sigma,
                          c.nstokes);
      }
    }
  }
  Vector levels(NL + 1, 0.0);
  for (Index l = 0; l < NL; l++) levels[l + 1] = v.tau()[l];
  const AscendingGrid      tau{std::move(levels)};
  rtepack::stokvec_tensor3 boundary;
  if (exact_boundary) {
    const Numeric Bs = planck(frequency, c.surface), Bsky = planck(frequency, c.sky);
    boundary.resize(2, 1, 2 * ne - first);
    boundary = rtepack::stokvec{};
    for (Index u = 0; u < 2 * ne - first; u++) {
      if (user_mu[u] < 0.0) {
        boundary[vdisort::cosine_mode, 0, u] = {Bsky, 0.0, 0.0, 0.0};
      } else if (const auto* f = std::get_if<rt4::fresnel_surface>(&c.g)) {
        const auto R                         = vdisort::brdf::Fresnel{f->refractive_index}(user_mu[u]);
        boundary[vdisort::cosine_mode, 0, u] = {(1 - R[0, 0]) * Bs, -R[1, 0] * Bs, -R[2, 0] * Bs, -R[3, 0] * Bs};
      } else {
        boundary[vdisort::cosine_mode, 0, u] = {
            (1 - std::get<rt4::lambertian_surface>(c.g).albedo) * Bs, 0.0, 0.0, 0.0};
      }
    }
  }
  rtepack::stokvec_tensor3 out(NL + 1, 1, 2 * ne - first);
  v.ungridded_u_user(out, tau, Vector{0.0}, user_mu, user_phase, {}, boundary);

  field f{.up = Tensor3(NL + 1, ne, c.nstokes, 0.0), .down = Tensor3(NL + 1, ne, c.nstokes, 0.0)};
  for (Index l = 0; l <= NL; l++) {
    for (Index e = 0; e < ne; e++) {
      for (Index s = 0; s < c.nstokes; s++) {
        if (upward) f.up[l, e, s] = out[l, 0, e][s];
        f.down[l, e, s] = out[l, 0, ne + e - first][s];
      }
    }
  }
  return f;
}

//! max |VDISORT - RT4| over levels, directions and streams, relative to max |I_RT4| there
struct deviation {
  Numeric I{}, Q{};

  Numeric max() const { return std::max(I, Q); }
};

enum class directions { both, up_only, down_only };

//! v holds RT4 streams first .. first + v.up.extent(1) - 1
deviation compare(const rt4::result& r, const field& v, Index first, directions which = directions::both) {
  const Index count = v.up.extent(1), ns = v.up.extent(2), nlev = v.up.extent(0);
  const bool  use_up = which != directions::down_only, use_down = which != directions::up_only;
  Numeric     scale = 0.0;
  for (Index l = 0; l < nlev; l++)
    for (Index k = 0; k < count; k++)
      scale = std::max({scale, std::abs(r.up[l, first + k, 0]), std::abs(r.down[l, first + k, 0])});
  deviation d;
  for (Index l = 0; l < nlev; l++) {
    for (Index k = 0; k < count; k++) {
      for (Index s = 0; s < ns; s++) {
        Numeric& worst = s == 0 ? d.I : d.Q;
        if (use_up) worst = std::max(worst, std::abs(v.up[l, k, s] - r.up[l, first + k, s]) / scale);
        if (use_down) worst = std::max(worst, std::abs(v.down[l, k, s] - r.down[l, first + k, s]) / scale);
      }
    }
  }
  return d;
}

//! max |Q| / max |I| of RT4 over levels, directions and streams first .. first + count - 1, to show Q is exercised
Numeric polarization(const rt4::result& r, Index first, Index count) {
  if (r.up.extent(2) < 2) return 0.0;
  Numeric I = 0.0, Q = 0.0;
  for (const auto* t : {&r.up, &r.down}) {
    for (Index l = 0; l < t->extent(0); l++) {
      for (Index i = first; i < first + count; i++) {
        I = std::max(I, std::abs((*t)[l, i, 0]));
        Q = std::max(Q, std::abs((*t)[l, i, 1]));
      }
    }
  }
  return Q / I;
}

/* RT4's doubling, as in radtran4.f.  A layer with optics is doubled
   n = int(F) + 1 times, with F = log2(max(tau, 1e-7) / max_delta_tau), or
   n = 0 if F <= 0.  The initial layer then has a vertical optical thickness
   delta = tau / 2^n.  tau is the extinction at the first stream plus the
   gas, times the thickness. */
Numeric doublings(Numeric tau, Numeric max_delta_tau) {
  const Numeric F = std::log2(std::max(tau, 1e-7) / max_delta_tau);
  return F > 0 ? std::floor(F) + 1 : 0.0;
}

struct doubling {
  bool    any{false};
  Numeric delta{0.0};  // the largest initial-layer thickness of the doubled layers
  Numeric n{0.0};      // the largest number of doublings
};

doubling rt4_doubling(const setup& c, Numeric max_delta_tau) {
  doubling d;
  for (Index l = 0; l < c.nlay(); l++) {
    if (c.optics_index[l] < 0) continue;
    const Numeric tau = (c.kp(l) + c.gas[l]) * c.dz(l);
    const Numeric n   = doublings(tau, max_delta_tau);
    d.any             = true;
    d.delta           = std::max(d.delta, tau / std::exp2(n));
    d.n               = std::max(d.n, n);
  }
  return d;
}

/* Tolerances, relative to max |I|.  delta is the largest initial-layer
   thickness and n the largest number of doublings.
   - No doubled layer.  Both solvers are exact up to round-off: 1e-12.
   - Direct, at max_delta_tau.  RT4's error is first order in delta.  The
     difference halves exactly when max_delta_tau halves (C11).  Over all
     cases here it is 0.1 to 0.7 delta, printed as "c delta".  The tolerance
     is 10 delta plus RT4's doubling round-off, which grows like 2^n eps (as
     in RT4's isothermal Kirchhoff test).
   - Richardson, 2 RT4(max_delta_tau / 2) - RT4(max_delta_tau) at
     max_delta_tau = 1e-5.  Every doubled layer gets exactly one more
     doubling (asserted), so this cancels the first-order term.  What
     remains is second order in delta plus three times the round-off.  The
     second-order term dominates at max_delta_tau = 1e-4 and 3e-5.  Measured
     there, it is 0.1 to 7 delta^2; 7 is for the forward-peaked C7.  At
     1e-5 and below the round-off dominates.  The tolerance is
     50 delta^2 + 3 2^n eps. */
constexpr Numeric eps                      = std::numeric_limits<Numeric>::epsilon();
constexpr Numeric richardson_max_delta_tau = 1e-5;

Numeric direct_tolerance(const doubling& d) { return d.any ? 10 * d.delta + std::exp2(d.n) * eps : 1e-12; }

Numeric richardson_tolerance(const setup& c) {
  const auto coarse = rt4_doubling(c, richardson_max_delta_tau);
  const auto fine   = rt4_doubling(c, richardson_max_delta_tau / 2);
  return 50 * coarse.delta * coarse.delta + 3 * std::exp2(fine.n) * eps;
}

//! RT4 extrapolated to max_delta_tau = 0 from richardson_max_delta_tau and half of it
rt4::result richardson(const setup& c) {
  for (Index l = 0; l < c.nlay(); l++) {
    if (c.optics_index[l] < 0) continue;
    const Numeric tau = (c.kp(l) + c.gas[l]) * c.dz(l);
    require(doublings(tau, richardson_max_delta_tau / 2) == doublings(tau, richardson_max_delta_tau) + 1,
            std::format("The Richardson extrapolation needs exactly one more doubling at max_delta_tau = {:.1e} "
                        "than at {:.1e}; a layer of optical thickness {} does not get it",
                        richardson_max_delta_tau / 2,
                        richardson_max_delta_tau,
                        tau));
  }
  auto a             = c;
  a.max_delta_tau    = richardson_max_delta_tau;
  const auto coarse  = rt4::solve(rt4_problem(a));
  a.max_delta_tau   /= 2;
  auto fine          = rt4::solve(rt4_problem(a));
  for (Index l = 0; l < fine.up.extent(0); l++) {
    for (Index i = 0; i < fine.up.extent(1); i++) {
      for (Index s = 0; s < fine.up.extent(2); s++) {
        fine.up[l, i, s]   = 2 * fine.up[l, i, s] - coarse.up[l, i, s];
        fine.down[l, i, s] = 2 * fine.down[l, i, s] - coarse.down[l, i, s];
      }
    }
  }
  return fine;
}

void report(std::string_view what, deviation d, Numeric tol, std::string_view note = {}) {
  std::cout << std::format("{:<66} I {:9.3e}  Q {:9.3e}  tolerance {:.1e}{}\n", what, d.I, d.Q, tol, note);
  require(d.max() <= tol,
          std::format("{}: VDISORT and RT4 must agree to max |difference| / max |I| <= {:.1e}, got {:.3e} for I "
                      "and {:.3e} for Q",
                      what,
                      tol,
                      d.I,
                      d.Q));
}

struct comparison {
  rt4::result        r;
  vdisort::main_data v;
  field              f;
};

comparison run(const setup& c) {
  auto r = rt4::solve(rt4_problem(c));
  auto v = vdisort_solver(c);
  auto f = vdisort_streams(v, c.nstokes);
  return {.r = std::move(r), .v = std::move(v), .f = std::move(f)};
}

//! Compare RT4 streams first .. first + v.up.extent(1) - 1 at max_delta_tau and, if doubled, Richardson-extrapolated
deviation check(std::string_view what, const setup& c, const rt4::result& r, const field& v, Index first) {
  const auto    dd = rt4_doubling(c, c.max_delta_tau);
  const auto    d  = compare(r, v, first);
  const Numeric qi = polarization(r, first, v.up.extent(1));
  report(what,
         d,
         direct_tolerance(dd),
         dd.any ? std::format(" = {:.2f} delta; max|Q|/max|I| {:.1e}", d.max() / dd.delta, qi)
                : std::format(", no doubling; max|Q|/max|I| {:.1e}", qi));
  if (dd.any)
    report("    RT4 Richardson-extrapolated to max_delta_tau = 0",
           compare(richardson(c), v, first),
           richardson_tolerance(c));
  return d;
}

deviation check(std::string_view what, const setup& c) {
  const auto x = run(c);
  return check(std::format("{} (nmu {})", what, c.s.nmu), c, x.r, x.f, 0);
}

//////////////////////////////////////////////////////////////////////////////
// Cases
//////////////////////////////////////////////////////////////////////////////

//! Three gas-only layers with a temperature gradient (heights in km, extinction per km)
setup gas_atmosphere(ground g) {
  setup c;
  c.s            = make_streams(8);
  c.height       = Vector{3.0, 2.0, 1.0, 0.0};
  c.temperature  = Vector{220.0, 240.0, 265.0, 285.0};
  c.gas          = Vector{0.1, 0.3, 0.5};
  c.optics_index = ArrayOfIndex{-1, -1, -1};
  c.surface      = 290.0;
  c.g            = g;
  return c;
}

//! C1-C3: RT4 integrates gas-only layers analytically, so both solvers are exact
void test_gas_only() {
  check("C1 gas-only, black surface", gas_atmosphere(rt4::lambertian_surface{.albedo = 0.0}));
  check("C2 gas-only, Fresnel n = 1.5", gas_atmosphere(rt4::fresnel_surface{.refractive_index = Complex{1.5, 0.0}}));
  check("C2 gas-only, Fresnel n = 3+0.2i", gas_atmosphere(rt4::fresnel_surface{.refractive_index = Complex{3.0, 0.2}}));
  check("C3 gas-only, Lambertian A = 0.3", gas_atmosphere(rt4::lambertian_surface{.albedo = 0.3}));
}

//! C4: one Rayleigh layer, omega 0.9, tau 1, over a black surface; Q comes from scattering only
void test_rayleigh_layer() {
  setup c;
  c.s            = make_streams(8);
  c.height       = Vector{1.0, 0.0};
  c.temperature  = Vector{220.0, 280.0};
  c.gas          = Vector{0.0};
  c.optics_index = ArrayOfIndex{0};
  c.optics       = {optics_set{.kp = 1.0, .sigma = 0.9, .Z = on_streams(c.s, 0.9, rayleigh)}};
  c.surface      = 285.0;
  check("C4 Rayleigh omega 0.9 tau 1, black surface", c);
}

/* C5: top-down, a gas-only layer (tau 0.3), then four scattering layers:
   - omega 0.5, tau 2, with gas;
   - omega 1 (exactly conservative), tau 5, no gas;
   - omega 0, tau 0.1, with gas (a doubled, non-scattering layer);
   - omega 0.95, tau 0.5, with gas.
   Rayleigh phase, over a Fresnel 3+0.2i surface. */
setup multilayer(Index nmu, const phase_function& P = rayleigh, Index nstokes = 2) {
  setup c;
  c.nstokes      = nstokes;
  c.s            = make_streams(nmu);
  c.height       = Vector{5.0, 4.0, 3.0, 2.0, 1.0, 0.0};
  c.temperature  = Vector{200.0, 225.0, 250.0, 262.0, 270.0, 285.0};
  c.gas          = Vector{0.3, 0.4, 0.0, 0.05, 0.01};
  c.optics_index = ArrayOfIndex{-1, 0, 1, 2, 3};
  const std::array<std::array<Numeric, 2>, 4> kp_sigma{{{1.6, 1.0}, {5.0, 5.0}, {0.05, 0.0}, {0.49, 0.475}}};
  for (const auto& [kp, sigma] : kp_sigma)
    c.optics.push_back({.kp = kp, .sigma = sigma, .Z = on_streams(c.s, sigma, P)});
  c.surface = 295.0;
  c.g       = rt4::fresnel_surface{.refractive_index = Complex{3.0, 0.2}};
  return c;
}

void test_multilayer() { check("C5 multilayer omega {0.5, 1, 0, 0.95}, Fresnel 3+0.2i", multilayer(8)); }

//! C6: a thick, exactly conservative Rayleigh layer over a warm Fresnel surface under a cold sky
void test_thick_conservative() {
  setup c;
  c.s            = make_streams(8);
  c.height       = Vector{1.0, 0.0};
  c.temperature  = Vector{210.0, 230.0};
  c.gas          = Vector{0.0};
  c.optics_index = ArrayOfIndex{0};
  c.optics       = {optics_set{.kp = 20.0, .sigma = 20.0, .Z = on_streams(c.s, 20.0, rayleigh)}};
  c.surface      = 300.0;
  c.g            = rt4::fresnel_surface{.refractive_index = Complex{1.5, 0.0}};
  check("C6 conservative Rayleigh tau 20, Fresnel 1.5, cold sky", c);
}

/* C7: a forward-peaked, polarizing phase matrix by numerical azimuthal
   averaging.  The builder is first validated against the Rayleigh closed
   form, all 16 elements, on the C7 streams plus two off-node angles,
   mu = 1 included. */
void test_numerical_phase_matrix() {
  constexpr Index nphi = 720;
  {
    const streams s    = make_streams(16, Vector{0.35, 1.0});
    const auto    num  = numeric_azimuthal_mean(s, 4 * pi, rayleigh_scattering_matrix, nphi);
    const auto    ref  = on_streams(s, 4 * pi, rayleigh);
    Numeric       diff = 0.0;
    for (Index ho = 0; ho < 2; ho++)
      for (Index hi = 0; hi < 2; hi++)
        for (Index io = 0; io < s.total(); io++)
          for (Index ii = 0; ii < s.total(); ii++)
            for (Index a = 0; a < 4; a++)
              for (Index b = 0; b < 4; b++)
                diff = std::max(diff, std::abs(num[ho, hi, io, ii][a, b] - ref[ho, hi, io, ii][a, b]));
    std::cout << std::format("{:<66} max |P_num - P_closed| {:9.3e}  tolerance 1.0e-12\n",
                             "C7 numerical azimuthal mean of Rayleigh, all 4 x 4 elements",
                             diff);
    require(diff <= 1e-12,
            std::format("The numerical azimuthal mean of the Rayleigh lab-frame matrix must equal the m = 0 closed "
                        "form to 1e-12 (nphi = {}), got {:.3e}",
                        nphi,
                        diff));
  }

  setup             c;
  constexpr Numeric sigma = 1.2;
  c.s                     = make_streams(16);
  auto Z                  = numeric_azimuthal_mean(c.s, sigma, polarized_henyey_greenstein, nphi);

  // RT4 requires mirror symmetry; VDISORT's m = 0 combined form drops the [I, Q] <-> [U, V] blocks
  Numeric asym = 0.0, cross_block = 0.0, zmax = 0.0;
  for (Index io = 0; io < c.s.total(); io++) {
    for (Index ii = 0; ii < c.s.total(); ii++) {
      for (Index a = 0; a < 4; a++) {
        for (Index b = 0; b < 4; b++) {
          zmax = std::max(zmax, std::abs(Z[down, down, io, ii][a, b]));
          asym = std::max({asym,
                           std::abs(Z[down, down, io, ii][a, b] - Z[up, up, io, ii][a, b]),
                           std::abs(Z[down, up, io, ii][a, b] - Z[up, down, io, ii][a, b])});
          if ((a < 2) != (b < 2))
            for (Index ho = 0; ho < 2; ho++)
              for (Index hi = 0; hi < 2; hi++) cross_block = std::max(cross_block, std::abs(Z[ho, hi, io, ii][a, b]));
        }
      }
    }
  }
  const Numeric renorm = renormalize(Z, c.s, sigma);
  std::cout << std::format(
      "    polarized HG g = 0.7: mirror asymmetry {:.1e} and [I, Q]-[U, V] coupling {:.1e} (relative to max Z); "
      "per-stream renormalisation up to {:.1e}\n",
      asym / zmax,
      cross_block / zmax,
      renorm);
  require(asym <= 1e-13 * zmax and cross_block <= 1e-13 * zmax,
          std::format("The numerical m = 0 polarized HG matrix must be mirror symmetric and have no [I, Q] - [U, V] "
                      "coupling to 1e-13 relative; got {:.1e} and {:.1e}",
                      asym / zmax,
                      cross_block / zmax));

  c.height       = Vector{2.0, 1.0, 0.0};
  c.temperature  = Vector{210.0, 240.0, 280.0};
  c.gas          = Vector{0.2, 0.15};
  c.optics_index = ArrayOfIndex{-1, 0};
  c.optics       = {optics_set{.kp = 1.35, .sigma = sigma, .Z = std::move(Z)}};
  c.surface      = 290.0;
  c.g            = rt4::fresnel_surface{.refractive_index = Complex{1.5, 0.0}};
  check("C7 polarized HG g 0.7 omega 0.8 tau 1.5, Fresnel 1.5", c);
}

/* C8: the constructed non-reciprocal, Stokes-asymmetric phase matrix.  The
   case is not blind: each deliberately wrong VDISORT phase input must miss
   RT4 by more than 100 times the tolerance. */
void test_nonreciprocal() {
  setup c;
  c.s            = make_streams(8);
  c.height       = Vector{2.0, 1.0, 0.0};
  c.temperature  = Vector{215.0, 245.0, 275.0};
  c.gas          = Vector{0.2, 0.1};
  c.optics_index = ArrayOfIndex{-1, 0};
  c.optics       = {optics_set{.kp = 0.9, .sigma = 0.8, .Z = on_streams(c.s, 0.8, constructed)}};
  c.surface      = 290.0;
  c.g            = rt4::fresnel_surface{.refractive_index = Complex{3.0, 0.2}};
  check("C8 non-reciprocal Stokes-asymmetric phase, Fresnel 3+0.2i", c);

  const auto    r   = rt4::solve(rt4_problem(c));
  const Numeric tol = direct_tolerance(rt4_doubling(c, c.max_delta_tau));
  for (auto m :
       {mistake::stream_transpose, mistake::stokes_transpose, mistake::full_transpose, mistake::quadrant_exchange}) {
    const auto d = compare(r, vdisort_streams(vdisort_solver(c, m), c.nstokes), 0);
    std::cout << std::format(
        "    deliberately wrong VDISORT phase: {:<44} I {:9.3e}  Q {:9.3e}  ({:.0f} x tolerance)\n",
        name(m),
        d.I,
        d.Q,
        d.max() / tol);
    require(d.max() > 100 * tol,
            std::format("C8 must detect a {} of the VDISORT phase matrix by more than 100 times the tolerance {:.1e}, "
                        "got {:.3e}",
                        name(m),
                        tol,
                        d.max()));
  }
}

/* C9: a non-specular, polarizing, non-reciprocal discrete surface under a
   Rayleigh layer.  rho is per steradian (VDISORT's R^0); the total
   reflectance is at most 0.5 (1.4)(0.4) = 0.28. */
rtepack::muelmat surface_rho(Numeric x, Numeric y) {
  rtepack::muelmat R{0.0};
  R[0, 0] = 0.5 / pi * (1 + 0.4 * x) * (1 - 0.3 * y);
  R[0, 1] = -0.05 / pi * x * (1 - y * y);
  R[1, 0] = 0.08 / pi * (1 - x * x) * (1 + y);
  R[1, 1] = 0.3 / pi * (0.5 + x * y * y);
  return R;
}

void test_discrete_surface() {
  setup c;
  c.s            = make_streams(8);
  c.height       = Vector{2.0, 1.0, 0.0};
  c.temperature  = Vector{220.0, 250.0, 270.0};
  c.gas          = Vector{0.15, 0.0};
  c.optics_index = ArrayOfIndex{-1, 0};
  c.optics       = {optics_set{.kp = 0.8, .sigma = 0.72, .Z = on_streams(c.s, 0.72, rayleigh)}};
  c.surface      = 285.0;
  c.g            = discrete_reflection{.rho = surface_rho};
  check("C9 discrete polarizing non-reciprocal surface, Rayleigh", c);

  const auto    r   = rt4::solve(rt4_problem(c));
  const Numeric tol = direct_tolerance(rt4_doubling(c, c.max_delta_tau));
  for (auto m : {mistake::stream_transpose, mistake::stokes_transpose}) {
    const auto d = compare(r, vdisort_streams(vdisort_solver(c, mistake::none, m), c.nstokes), 0);
    std::cout << std::format(
        "    deliberately wrong VDISORT surface: {:<42} I {:9.3e}  Q {:9.3e}  ({:.0f} x tolerance)\n",
        name(m),
        d.I,
        d.Q,
        d.max() / tol);
    require(d.max() > 100 * tol,
            std::format("C9 must detect a {} of the VDISORT surface by more than 100 times the tolerance {:.1e}, got "
                        "{:.3e}",
                        name(m),
                        tol,
                        d.max()));
  }
}

//! C10: C5 for nmu = 2 ... 32 (RT4 needs nstokes nmu <= 64)
void test_stream_sweep() {
  for (Index nmu : {2, 4, 8, 16, 32}) check("C10 C5 stream sweep", multilayer(nmu));
}

/* C11: the C5 difference must shrink with max_delta_tau like RT4's
   first-order doubling error.  The number of doublings is
   int(log2(tau / max_delta_tau)) + 1:
   - halving max_delta_tau adds exactly one doubling to every layer, so the
     difference must halve;
   - a decade changes delta by a factor of 8 or 16. */
void test_convergence() {
  auto                 c = multilayer(8);
  const auto           v = vdisort_streams(vdisort_solver(c), c.nstokes);
  std::vector<Numeric> err;
  for (Numeric mdt : {1e-5, 5e-6, 1e-6, 1e-7}) {
    c.max_delta_tau = mdt;
    const auto dd   = rt4_doubling(c, mdt);
    const auto d    = compare(rt4::solve(rt4_problem(c)), v, 0);
    err.push_back(d.max());
    report(std::format("C11 C5 convergence, max_delta_tau {:.0e}", mdt),
           d,
           direct_tolerance(dd),
           std::format(" = {:.2f} delta", d.max() / dd.delta));
  }
  const Numeric halving = err[0] / err[1];
  std::cout << std::format("    difference ratio for max_delta_tau 1e-5 / 5e-6: {:.3f}\n", halving);
  require(halving > 1.9 and halving < 2.1,
          std::format("Halving max_delta_tau adds one doubling to every layer, so the VDISORT-RT4 difference must "
                      "halve (ratio 1.9 to 2.1); got {:.3f}",
                      halving));
  for (std::size_t i : {std::size_t{0}, std::size_t{2}}) {
    const Numeric ratio = err[i] / err[i + (i == 0 ? 2 : 1)];
    std::cout << std::format("    difference ratio per decade of max_delta_tau: {:.2f}\n", ratio);
    require(ratio > 4.0 and ratio < 25.0,
            std::format("The VDISORT-RT4 difference must fall like RT4's first-order doubling error, by a factor "
                        "4 to 25 per decade of max_delta_tau; got {:.2f}",
                        ratio));
  }
}

//! C12: nstokes = 1 on the C5 geometry with the intensity-only Rayleigh phase function
void test_scalar() { check("C12 C5 geometry, nstokes 1, Rayleigh P_II", multilayer(8, rayleigh_intensity_only, 1)); }

/* Extra angles: RT4's zero-weight streams against VDISORT's user-angle
   formal solution, at mu = 0.35 and 1, which are not nodes.  Black and
   Lambertian surfaces are asserted.  For the Fresnel surface only the
   downward radiances are asserted; the upward ones are reported, see
   below. */
setup extra_angle_setup(ground g) {
  setup c;
  c.s            = make_streams(8, Vector{0.35, 1.0});
  c.height       = Vector{2.0, 1.0, 0.0};
  c.temperature  = Vector{210.0, 230.0, 280.0};
  c.gas          = Vector{0.2, 0.0};
  c.optics_index = ArrayOfIndex{-1, 0};
  c.optics       = {optics_set{.kp = 1.0, .sigma = 0.9, .Z = on_streams(c.s, 0.9, rayleigh)}};
  c.surface      = 290.0;
  c.g            = g;
  return c;
}

void test_extra_angles() {
  for (const auto& [what, g] : {std::pair<std::string_view, ground>{"black", rt4::lambertian_surface{.albedo = 0.0}},
                                {"Lambertian A = 0.3", rt4::lambertian_surface{.albedo = 0.3}}}) {
    const auto c = extra_angle_setup(g);
    const auto x = run(c);
    check(std::format("Extra-angle setup, {}: quadrature streams", what), c, x.r, x.f, 0);
    check(std::format("Extra-angle setup, {}: mu = 0.35, 1 (user angles)", what),
          c,
          x.r,
          vdisort_extra_angles(x.v, c),
          c.s.nmu);
  }

  /* The Fresnel surface.  VDISORT reflects the downward user-angle radiance
     at -mu into mu (BDRF::specular), so it needs both directions; without
     the downward partner it must refuse an upward user angle.  The surface
     emission at a user angle is interpolated from the streams (barycentric
     in mu, disort_common::barycentric_interpolate), so VDISORT's upward
     radiance carries the interpolation error d(mu) = interpolated - exact
     emission [(1 - R11) B_s, -R21 B_s], attenuated by exp(-(tau_s - tau) / mu).
     With that term removed, every level must agree with RT4 to its doubling
     error: the reflection R(mu) I_down(mu) is then exact. */
  const auto c = extra_angle_setup(rt4::fresnel_surface{.refractive_index = Complex{3.0, 0.2}});
  const auto x = run(c);
  check("Extra-angle setup, Fresnel 3+0.2i: quadrature streams", c, x.r, x.f, 0);
  const auto    vu  = vdisort_extra_angles(x.v, c);
  const auto    tol = direct_tolerance(rt4_doubling(c, c.max_delta_tau));
  const Index   N = c.s.nmu, L = c.nlay(), ne = c.s.total() - N;
  const Numeric Bs = planck(frequency, c.surface);
  report("Extra-angle setup, Fresnel 3+0.2i: mu = 0.35, 1, downward", compare(x.r, vu, N, directions::down_only), tol);

  const Vector nodes{x.v.mu()[Range(0, N)]};
  Vector       weights(N), emission_i(N), emission_q(N);
  disort_common::barycentric_weights(weights, nodes);
  const vdisort::brdf::Fresnel fresnel{.refractive_index = Complex{3.0, 0.2}};
  for (Index i = 0; i < N; i++) {
    const auto R  = fresnel(nodes[i]);
    emission_i[i] = (1 - R[0, 0]) * Bs;
    emission_q[i] = -R[1, 0] * Bs;
  }
  Numeric scale = 0.0, raw = 0.0, corrected = 0.0, interpolation = 0.0;
  for (Index l = 0; l <= L; l++)
    for (Index e = 0; e < ne; e++)
      scale = std::max({scale, std::abs(x.r.up[l, N + e, 0]), std::abs(x.r.down[l, N + e, 0])});
  for (Index e = 0; e < ne; e++) {
    const Numeric mu  = c.s.mu[N + e];
    const auto    R   = fresnel(mu);
    const Numeric d_i = disort_common::barycentric_interpolate(nodes, weights, emission_i, mu) - (1 - R[0, 0]) * Bs;
    const Numeric d_q = disort_common::barycentric_interpolate(nodes, weights, emission_q, mu) + R[1, 0] * Bs;
    interpolation     = std::max({interpolation, std::abs(d_i) / Bs, std::abs(d_q) / Bs});
    for (Index l = 0; l <= L; l++) {
      const Numeric tau = l == 0 ? 0.0 : x.v.tau()[l - 1];
      const Numeric att = std::exp(-(x.v.tau()[L - 1] - tau) / mu);
      for (Index st = 0; st < std::min<Index>(c.nstokes, 2); st++) {
        const Numeric diff = vu.up[l, e, st] - x.r.up[l, N + e, st];
        raw                = std::max(raw, std::abs(diff) / scale);
        corrected          = std::max(corrected, std::abs(diff - (st == 0 ? d_i : d_q) * att) / scale);
      }
    }
  }
  std::cout << std::format(
      "{:<66} I, Q {:9.3e} after removing the emission interpolation ({:.1e} of B_s; {:.1e} of max I before); "
      "tolerance {:.1e}\n",
      "Extra-angle setup, Fresnel 3+0.2i: mu = 0.35, 1, upward",
      corrected,
      interpolation,
      raw,
      tol);
  require(corrected <= tol,
          std::format("VDISORT's upward user-angle radiance over a Fresnel surface must be RT4's up to its emission "
                      "interpolation, to {:.1e}; got {:.2e}",
                      tol,
                      corrected));

  // With the exact boundary radiances given there is no interpolation: VDISORT is RT4 at every level
  const auto exact = vdisort_extra_angles(x.v, c, true, true);
  report("Extra-angle setup, Fresnel 3+0.2i: mu = 0.35, 1, exact boundaries given", compare(x.r, exact, N), tol);

  // Without the downward partner an upward user angle must be refused
  bool refused = false;
  try {
    rtepack::stokvec_tensor3 out(1, 1, 1);
    x.v.ungridded_u_user(out,
                         AscendingGrid{Vector{0.0}},
                         Vector{0.0},
                         Vector{c.s.mu[N]},
                         vdisort::phase_matrix_data(2, 1, L, 1, 2 * N, rtepack::muelmat{0.0}));
  } catch (const std::exception&) { refused = true; }
  std::cout << std::format("{:<66} {}\n",
                           "Extra-angle setup, Fresnel 3+0.2i: upward without its downward partner",
                           refused ? "refused, as required" : "NOT REFUSED");
  require(refused, "VDISORT must refuse an upward user angle over a Fresnel surface without its downward partner");
}

}  // namespace

int main() try {
  require(rt4::available(), "This test requires ENABLE_RT4=ON");
  test_gas_only();
  test_rayleigh_layer();
  test_multilayer();
  test_thick_conservative();
  test_numerical_phase_matrix();
  test_nonreciprocal();
  test_discrete_surface();
  test_stream_sweep();
  test_convergence();
  test_scalar();
  test_extra_angles();
  std::cout << "vdisort-rt4 comparison passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
