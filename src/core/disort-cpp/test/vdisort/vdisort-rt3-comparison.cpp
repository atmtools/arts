/* VDISORT against Evans' RT3 polarized doubling-adding solver (src/core/rt3).

   Where the expected answer comes from.  For every Fourier azimuth mode,
   RT3 builds the reflection, transmission and source operators of each
   layer by doubling a thin initial layer, and combines the layers by
   adding (Evans and Stephens 1991).  It makes its own Fourier modes of the
   phase matrix from the Legendre series of the six scattering-plane
   elements, with its own rotation formulas (ROTATE_PHASE_MATRIX) and FFT,
   and its own direct-beam pseudo-source.  That is a framework external to
   VDISORT, which eigen-decomposes the combined cosine and sine
   discrete-ordinate systems of each mode and matches boundary and
   continuity conditions.  Both solvers get the same discrete problem:

   - RT3's double-Gauss rule with nmu nodes per hemisphere, which is also
     VDISORT's rule with NQuad = 2 nmu streams (asserted);
   - the same Legendre series of F11, F12, F22, F33, F34 and F44.  This test
     builds VDISORT's Fourier coefficients from it by vector geometry
     (lab-frame.h), at RT3's azimuth samples;
   - the same scalar extinction and single-scattering albedo per layer;
   - the same solar beam, Planck source linear in optical depth, sky and
     Lambertian surface.

   With identical angular discretisation, both solve the same linear system
   of ODEs in optical depth for every mode.  VDISORT solves it to round-off.
   RT3 solves it to round-off in gas-only layers, which it integrates
   analytically.  In every other layer, RT3's error is that of its
   first-order initial layer, so it is first order in that layer's
   thickness delta <= max_delta_tau.  The difference must therefore vanish
   linearly with max_delta_tau (R7), and a Richardson extrapolation of RT3
   in max_delta_tau must agree with VDISORT to second order and round-off.
   A convention mismatch (azimuth origin or sense, U or V sign, Fourier
   normalisation, cosine/sine system, beam normalisation or source) would
   not converge.

   Mapping, layer l top-down, with thickness dz, gas extinction kg and the
   scattering set (k, sigma, legendre) of the layer:

   - RT3: the set as given, gas_extinction kg, double_gauss, aziorder =
     NFourier - 1, direct_flux F on the horizontal at direct_mu = mu0,
     Lambertian albedo A.
   - VDISORT:
     - tau is the cumulative (kg + k) dz and omega = sigma / (kg + k);
     - the ordinary Fourier coefficients, without epsilon_m, on the signed
       streams (> 0 upward),
         C^m(mu_o, mu_i) = (1 / 2 pi) int Z(mu_o, 0; mu_i, phi') cos(m phi') dphi',
         S^m(mu_o, mu_i) = (1 / 2 pi) int Z(mu_o, 0; mu_i, phi') sin(m phi') dphi',
       with Z the lab-frame phase matrix of the series; the diffuse
       operator is vdisort::combine_phase_matrices(C, S);
     - the beam operator is vdisort::combine_beam_phase_matrices of
       C^m(mu_i, -mu0) and S^m(mu_i, -mu0), checked against combine_beam()
       below;
     - beam_stokes = [F / mu0, 0, 0, 0], the irradiance normal to the beam;
     - stream i < N is RT3 (up, mu_i), stream N + i is RT3 (down, mu_i);
     - the radiance at azimuth phi0 + psi is RT3's at psi: both azimuths
       are those of the propagation direction, and VDISORT's beam
       propagates toward phi0, RT3's toward 0;
     - the thermal source is c0 + c1 tau in global tau, the sky is
       planck(T_sky), the surface emits (1 - A) planck(T_s), and
       brdf::lambertian_fourier_modes(A, NFourier).
   - With delta-M (R9), RT3 scales every set (GET_SCAT_SET) and VDISORT is
     given the scaled set, computed here with RT3's algebra.  That is not
     VDISORT's own IMS/TMS-corrected delta-M, which is out of scope.

   The tolerances are derived above direct_tolerance() below. */
#include <arts_constants.h>
#include <legendre.h>
#include <physics_funcs.h>
#include <rt3.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "evans-scripts.h"
#include "lab-frame.h"

namespace {
constexpr Numeric pi  = Constant::pi;
constexpr Numeric eps = std::numeric_limits<Numeric>::epsilon();

using vdisort_test::lab_frame;
using vdisort_test::tro_elements;
using vdisort_test::tro_matrix;

Index isize(const auto& v) { return static_cast<Index>(v.size()); }

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

//! The RT3 test azimuths psi [deg], principal plane and off it
constexpr std::array<Numeric, 7> azimuths_deg{0.0, 30.0, 75.0, 90.0, 135.0, 180.0, 250.0};

Vector test_azimuths() {
  Vector psi(isize(azimuths_deg));
  for (Index k = 0; k < isize(azimuths_deg); k++) psi[k] = azimuths_deg[k] * pi / 180.0;
  return psi;
}

//! An azimuth in [0, 2 pi), as VDISORT requires
Numeric wrapped(Numeric phi) {
  Numeric x = std::fmod(phi, 2 * pi);
  if (x < 0) x += 2 * pi;
  return x >= 2 * pi ? 0.0 : x;
}

//////////////////////////////////////////////////////////////////////////////
// Scattering sets (columns F11, F12, F33, F34, F22, F44, as RT3's files)
//////////////////////////////////////////////////////////////////////////////

Matrix legendre(std::initializer_list<std::array<Numeric, 6>> rows) {
  Matrix m(isize(rows), 6);
  Index  l = 0;
  for (const auto& r : rows) {
    for (Index k = 0; k < 6; k++) m[l, k] = r[k];
    l++;
  }
  return m;
}

//! mietest.sca of Evans' runmietest and runtesta (3rdparty/polradtran)
Matrix mie_legendre() {
  return legendre({
      {1.00000000, -.32071711, .71206342, -.01882245, 1.00000000, .71206342},
      {1.45529318, -.20350675, 1.76014119, -.04725108, 1.45529318, 1.76014119},
      {1.05402631, .24638948, 1.06682431, .00894436, 1.05402631, 1.06682431},
      {.39758994, .18605748, .39651104, .04505815, .39758994, .39651104},
      {.11659302, .07124848, .09576412, .00958275, .11659302, .09576412},
      {.02387477, .01700757, .01765088, .00215761, .02387477, .01765088},
      {.00395010, .00302534, .00261549, .00029195, .00395010, .00261549},
      {.00053888, .00043592, .00032713, .00003502, .00053888, .00032713},
      {.00006372, .00005326, .00003583, .00000337, .00006372, .00003583},
      {.00000667, .00000572, .00000351, .00000029, .00000667, .00000351},
      {.00000063, .00000055, .00000031, .00000002, .00000063, .00000031},
      {.00000006, .00000005, .00000003, .00000000, .00000006, .00000003},
  });
}

//! rayleigh.sca of runtesta: F11 = 3/4 (1 + x^2), F12 = -3/4 (1 - x^2), F33 = F44 = 3/2 x
Matrix rayleigh_legendre() {
  return legendre({{1.0, -0.5, 0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 1.5, 0.0, 0.0, 1.5}, {0.5, 0.5, 0.0, 0.0, 0.5, 0.0}});
}

//! The first degree + 1 rows of a series
Matrix truncated(const Matrix& coef, Index degree) {
  const Index n = std::min(degree + 1, coef.nrows());
  Matrix      m(n, 6);
  for (Index l = 0; l < n; l++)
    for (Index k = 0; k < 6; k++) m[l, k] = coef[l, k];
  return m;
}

//! The highest row with a non-zero coefficient (0 if none), as rt3::solve strips the series
Index stripped_degree(const Matrix& coef) {
  for (Index l = coef.nrows() - 1; l > 0; l--)
    for (Index k = 0; k < 6; k++)
      if (coef[l, k] != 0.0) return l;
  return 0;
}

//! F(cos Theta) of a series, each element a plain Legendre series (RT3's SUM_LEGENDRE)
tro_matrix legendre_matrix(const Matrix& coef) {
  return [coef](Numeric x) {
    Vector p(coef.nrows()), sum(6, 0.0);
    Legendre::legendre_polynomials(p, x);
    for (Index l = 0; l < coef.nrows(); l++)
      for (Index k = 0; k < 6; k++) sum[k] += coef[l, k] * p[l];
    return tro_elements{.F11 = sum[0], .F12 = sum[1], .F22 = sum[4], .F33 = sum[2], .F34 = sum[3], .F44 = sum[5]};
  };
}

//////////////////////////////////////////////////////////////////////////////
// The problem
//////////////////////////////////////////////////////////////////////////////

constexpr Numeric solar_frequency = Constant::c / 0.5e-6;

struct setup {
  Index                            nstokes{4};
  Index                            nmu{8};
  Index                            aziorder{3};
  Numeric                          max_delta_tau{1e-7};
  bool                             delta_m{false};  // RT3's delta-M; VDISORT gets the scaled sets
  Numeric                          frequency{solar_frequency};
  Vector                           height{};       // [nlay + 1], top-down; extinctions are per unit of it
  Vector                           temperature{};  // [nlay + 1] at the interfaces
  Vector                           gas{};          // [nlay]
  ArrayOfIndex                     set_index{};    // [nlay], < 0 for gas-only
  std::vector<rt3::scattering_set> sets{};
  Numeric                          direct_flux{0.0};  // on the horizontal at the top
  Numeric                          mu0{0.5};
  bool                             thermal{false};
  Numeric                          sky{0.0};
  Numeric                          surface{0.0};
  Numeric                          albedo{0.0};
  std::optional<Complex>           fresnel{};  // a Fresnel surface of this index instead of the Lambertian albedo
  Numeric                          phi0{0.0};  // VDISORT's beam azimuth

  Index   nlay() const { return isize(height) - 1; }
  Index   nfourier() const { return aziorder + 1; }
  Numeric dz(Index l) const { return std::abs(height[l] - height[l + 1]); }
};

/* A scattering set as RT3 transports it: trailing zero rows dropped (as
   rt3::solve does), and with delta-M scaled exactly as GET_SCAT_SET
   (radscat3.f), M = 2 nmu:
     f = legendre[M, 0] / (2 M + 1), k' = (1 - w f) k with w = sigma / k,
     w' = (1 - f) w / (1 - w f), sigma' = w' k',
     diagonal (F11, F33, F22, F44): (2 l + 1) (c_l / (2 l + 1) - f) / (1 - f),
     off-diagonal (F12, F34): (2 l + 1) (c_l / (2 l + 1)) / (1 - f),
   for l <= M - 1, then truncated to RT3's NLEGLIM.  rt3::solve rejects a
   problem in which that truncation would drop a non-zero coefficient. */
rt3::scattering_set transport_set(const rt3::scattering_set& s, const setup& c) {
  const Index nleg = stripped_degree(s.legendre);
  if (not c.delta_m)
    return {.extinction = s.extinction, .scattering = s.scattering, .legendre = truncated(s.legendre, nleg)};

  const Index   M      = 2 * c.nmu;
  const Numeric f      = M <= nleg ? s.legendre[M, 0] / static_cast<Numeric>(2 * M + 1) : 0.0;
  const Numeric w      = s.scattering / s.extinction;
  const Numeric ext    = (1 - w * f) * s.extinction;
  const Numeric sca    = (1 - f) * w / (1 - w * f) * ext;
  const Index   degree = std::min(M - 1, rt3::max_legendre_degree(c.nmu, rt3::quadrature_type::double_gauss));
  Matrix        L(degree + 1, 6, 0.0);
  for (Index l = 0; l <= degree; l++) {
    const auto m = static_cast<Numeric>(2 * l + 1);
    for (Index k = 0; k < 6; k++) {
      const Numeric x    = l <= nleg ? s.legendre[l, k] : 0.0;
      const bool    diag = k == 0 or k == 2 or k == 4 or k == 5;
      L[l, k]            = diag ? m * (x / m - f) / (1 - f) : m * (x / m) / (1 - f);
    }
  }
  return {.extinction = ext, .scattering = sca, .legendre = std::move(L)};
}

rt3::problem rt3_problem(const setup& c) {
  rt3::problem p;
  p.nstokes                = c.nstokes;
  p.nmu                    = c.nmu;
  p.quad                   = rt3::quadrature_type::double_gauss;
  p.aziorder               = c.aziorder;
  p.max_delta_tau          = c.max_delta_tau;
  p.delta_m                = c.delta_m;
  p.direct_flux            = c.direct_flux;
  p.direct_mu              = c.mu0;
  p.thermal                = c.thermal;
  p.frequency              = c.frequency;
  p.height                 = c.height;
  p.temperature            = c.temperature;
  p.gas_extinction         = c.gas;
  p.scattering_sets        = c.sets;
  p.layer_scattering_index = c.set_index;
  p.sky_temperature        = c.sky;
  p.surface_temperature    = c.surface;
  if (c.fresnel)
    p.ground = rt3::fresnel_surface{.refractive_index = *c.fresnel};
  else
    p.ground = rt3::lambertian_surface{.albedo = c.albedo};
  return p;
}

//////////////////////////////////////////////////////////////////////////////
// VDISORT's Fourier coefficients
//////////////////////////////////////////////////////////////////////////////

/* RT3's number of azimuth samples of the phase matrix for a series summed
   to degree L (SCATTERING and DIRECT_SCATTERING in radscat3.f), with the
   same floating-point expression:
     aziorder > 0: 2 * 2^int(ln(L + 4) / ln(2) + 1),
     aziorder = 0: 2 int((L + 1) / 2) + 4. */
Index rt3_azimuth_samples(Index degree, Index aziorder) {
  if (aziorder == 0) return 2 * ((degree + 1) / 2) + 4;
  const auto e = static_cast<int>(std::log(static_cast<Numeric>(degree + 4)) / std::log(2.0) + 1.0);
  return 2 * (Index{1} << e);
}

struct ordinary_modes {
  rtepack::muelmat_tensor3 C;  // [m, out, in]
  rtepack::muelmat_tensor3 S;
};

/* C^m and S^m on signed cosines (> 0 upward), by the trapezoid rule at the
   nphi equidistant phi' = 2 pi k / nphi, k = 0 .. nphi - 1.  These are
   RT3's samples, including phi' = 0 and pi, and RT3's FFT is the same
   trapezoid rule; RT3 keeps the modes m < nphi / 2.  For a band-limited Z
   (Rayleigh, or any series that is regular at Theta = 0 and pi) the
   samples do not matter; Evans' Mie series is regular only to 1e-8, which
   changes the coefficients by 3e-10 relative between sample sets (printed
   by V0).
   Rows and columns from nstokes on are zeroed before the transform: RT3
   with nstokes < 4 transports the leading nstokes x nstokes block.
   midpoints shifts the samples by half a step (for V0), and
   azimuth_sign = -1 takes Z at -phi' (a deliberate mistake). */
ordinary_modes fourier_modes(const tro_matrix& F,
                             const Vector&     mu_out,
                             const Vector&     mu_in,
                             Index             nfourier,
                             Index             nphi,
                             Index             nstokes,
                             bool              midpoints    = false,
                             Numeric           azimuth_sign = 1.0) {
  const Index    no = isize(mu_out), ni = isize(mu_in);
  ordinary_modes r{.C = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0}),
                   .S = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0})};
  const Index    nmode = std::min(nfourier, nphi / 2);
  const auto     n     = static_cast<Numeric>(nphi);
  for (Index io = 0; io < no; io++) {
    for (Index ii = 0; ii < ni; ii++) {
      for (Index k = 0; k < nphi; k++) {
        const Numeric phi = (2 * pi * (static_cast<Numeric>(k) + (midpoints ? 0.5 : 0.0))) / n;
        auto          Z   = lab_frame(F, mu_in[ii], azimuth_sign * phi, mu_out[io], 0.0);
        for (Index a = 0; a < 4; a++)
          for (Index b = 0; b < 4; b++)
            if (a >= nstokes or b >= nstokes) Z[a, b] = 0.0;
        for (Index m = 0; m < nmode; m++) {
          r.C[m, io, ii] += (std::cos(static_cast<Numeric>(m) * phi) / n) * Z;
          r.S[m, io, ii] += (std::sin(static_cast<Numeric>(m) * phi) / n) * Z;
        }
      }
    }
  }
  return r;
}

/* The combined beam operator.  The beam is a delta in azimuth at phi0, with
   the expansion (1 / 2 pi) sum_m eps_m cos m(phi0 - phi): cosine terms
   only, for every Stokes component of the beam.  Its source is
     Z(mu, phi; -mu0, phi0) S_b = sum_m eps_m [C^m cos m(phi0 - phi) + S^m sin m(phi0 - phi)] S_b,
   so the cosine system (I^c, Q^c, U^s, V^s) gets rows I, Q of C^m and rows
   U, V of S^m, and the sine system (I^s, Q^s, U^c, V^c) rows I, Q of S^m
   and rows U, V of C^m, in all four columns.  VDISORT applies eps_m
   itself.  vdisort::combine_beam_phase_matrices, which the solver set-up
   below uses, must equal this exactly (V0).  It used to reuse the diffuse
   combination of Lin et al. Eq. 81, which this comparison exposed: 7 to 58
   per cent errors for m >= 1. */
vdisort::beam_phase_matrix_data combine_beam(const rtepack::muelmat_tensor3& C, const rtepack::muelmat_tensor3& S) {
  const auto [NF, NL, NQ] = C.shape();
  vdisort::beam_phase_matrix_data B(2, NF, NL, NQ, rtepack::muelmat{0.0});
  for (Index m = 0; m < NF; m++) {
    for (Index l = 0; l < NL; l++) {
      for (Index i = 0; i < NQ; i++) {
        for (Index a = 0; a < 4; a++) {
          for (Index b = 0; b < 4; b++) {
            B[vdisort::cosine_mode, m, l, i][a, b] = a < 2 ? C[m, l, i][a, b] : S[m, l, i][a, b];
            B[vdisort::sine_mode, m, l, i][a, b]   = a < 2 ? S[m, l, i][a, b] : C[m, l, i][a, b];
          }
        }
      }
    }
  }
  return B;
}

//! Deliberate VDISORT input errors, to show that the comparison detects them
enum class mistake {
  none,
  sine_sign,         // S^m -> -S^m, diffuse and beam: a U-sign error
  mirrored_azimuth,  // Z(mu_o, 0; mu_i, -phi') in the transform, diffuse and beam: an azimuth-sense error
  diffuse_epsilon,   // diffuse C^m and S^m (m > 0) multiplied by eps_m = 2
  beam_epsilon,      // beam C^m and S^m (m > 0) multiplied by eps_m = 2, which VDISORT applies itself
  swapped_systems,   // the cosine (I^c, Q^c, U^s, V^s) and sine (I^s, Q^s, U^c, V^c) systems exchanged
};

std::string_view name(mistake m) {
  switch (m) {
    case mistake::none:             return "none";
    case mistake::sine_sign:        return "S^m sign flipped (diffuse and beam)";
    case mistake::mirrored_azimuth: return "Z at -phi' in the transform";
    case mistake::diffuse_epsilon:  return "diffuse C^m, S^m (m > 0) times eps = 2";
    case mistake::beam_epsilon:     return "beam C^m, S^m (m > 0) times eps = 2";
    case mistake::swapped_systems:  return "cosine and sine systems swapped";
  }
  return "";
}

//! The signed VDISORT streams on RT3's double-Gauss nodes: mu_i, then -mu_i
Vector signed_streams(Index nmu) {
  const auto q = rt3::get_quadrature(nmu, rt3::quadrature_type::double_gauss);
  Vector     mu(2 * nmu);
  for (Index i = 0; i < nmu; i++) {
    mu[i]       = q.mu[i];
    mu[nmu + i] = -q.mu[i];
  }
  return mu;
}

vdisort::main_data vdisort_solver(const setup& c, mistake mk = mistake::none) {
  const Index  N = c.nmu, NQuad = 2 * N, NL = c.nlay(), NF = c.nfourier(), ns = c.nstokes;
  const Vector mu = signed_streams(N);

  std::vector<rt3::scattering_set> T;
  for (const auto& s : c.sets) T.push_back(transport_set(s, c));

  Vector  tau(NL), omega(NL);
  Numeric t = 0.0;
  for (Index l = 0; l < NL; l++) {
    const Index   s  = c.set_index[l];
    const Numeric k  = c.gas[l] + (s < 0 ? 0.0 : T[s].extinction);
    t               += k * c.dz(l);
    tau[l]           = t;
    omega[l]         = s < 0 ? 0.0 : T[s].scattering / k;
  }

  std::vector<ordinary_modes> diffuse, beam;
  for (const auto& s : T) {
    const auto    F    = legendre_matrix(s.legendre);
    const Index   nphi = rt3_azimuth_samples(s.legendre.nrows() - 1, c.aziorder);
    const Numeric sign = mk == mistake::mirrored_azimuth ? -1.0 : 1.0;
    diffuse.push_back(fourier_modes(F, mu, mu, NF, nphi, ns, false, sign));
    beam.push_back(fourier_modes(F, mu, Vector{-c.mu0}, NF, nphi, ns, false, sign));
  }

  // The ordinary coefficients per layer, [m, layer, out, in] and [m, layer, out] for the beam
  rtepack::muelmat_tensor4 Cd(NF, NL, NQuad, NQuad, rtepack::muelmat{0.0});
  rtepack::muelmat_tensor4 Sd(NF, NL, NQuad, NQuad, rtepack::muelmat{0.0});
  rtepack::muelmat_tensor3 Cb(NF, NL, NQuad, rtepack::muelmat{0.0});
  rtepack::muelmat_tensor3 Sb(NF, NL, NQuad, rtepack::muelmat{0.0});
  for (Index l = 0; l < NL; l++) {
    const Index s = c.set_index[l];
    if (s < 0) continue;
    for (Index m = 0; m < NF; m++) {
      const Numeric diffuse_factor = m > 0 and mk == mistake::diffuse_epsilon ? 2.0 : 1.0;
      const Numeric beam_factor    = m > 0 and mk == mistake::beam_epsilon ? 2.0 : 1.0;
      const Numeric sine_factor    = mk == mistake::sine_sign ? -1.0 : 1.0;
      for (Index i = 0; i < NQuad; i++) {
        for (Index j = 0; j < NQuad; j++) {
          Cd[m, l, i, j] = diffuse_factor * diffuse[s].C[m, i, j];
          Sd[m, l, i, j] = (sine_factor * diffuse_factor) * diffuse[s].S[m, i, j];
        }
        Cb[m, l, i] = beam_factor * beam[s].C[m, i, 0];
        Sb[m, l, i] = (sine_factor * beam_factor) * beam[s].S[m, i, 0];
      }
    }
  }
  auto P = vdisort::combine_phase_matrices(Cd, Sd);
  auto B = vdisort::combine_beam_phase_matrices(Cb, Sb);
  if (mk == mistake::swapped_systems) {
    for (Index m = 0; m < NF; m++) {
      for (Index l = 0; l < NL; l++) {
        for (Index i = 0; i < NQuad; i++) {
          for (Index j = 0; j < NQuad; j++)
            std::swap(P[vdisort::cosine_mode, m, l, i, j], P[vdisort::sine_mode, m, l, i, j]);
          std::swap(B[vdisort::cosine_mode, m, l, i], B[vdisort::sine_mode, m, l, i]);
        }
      }
    }
  }

  rtepack::stokvec_tensor3 bottom(2, NF, N), top(2, NF, N);
  bottom              = rtepack::stokvec{};
  top                 = rtepack::stokvec{};
  const Numeric Bsky  = c.sky > 0 ? planck(c.frequency, c.sky) : 0.0;
  const Numeric Bs    = c.surface > 0 ? planck(c.frequency, c.surface) : 0.0;
  const Numeric Bsurf = c.thermal ? (1 - c.albedo) * Bs : 0.0;
  for (Index i = 0; i < N; i++) {
    top[vdisort::cosine_mode, 0, i] = {Bsky, 0.0, 0.0, 0.0};
    if (c.fresnel) {
      // Kirchhoff emission B ([1, 0, 0, 0] - R(mu)[:, 0]), always on, as RT3's Fresnel surface
      const auto R                       = vdisort::brdf::Fresnel{.refractive_index = *c.fresnel}(mu[i]);
      bottom[vdisort::cosine_mode, 0, i] = {(1 - R[0, 0]) * Bs, -R[1, 0] * Bs, -R[2, 0] * Bs, -R[3, 0] * Bs};
    } else {
      bottom[vdisort::cosine_mode, 0, i] = {Bsurf, 0.0, 0.0, 0.0};
    }
  }

  // B(tau) = c0 + c1 tau in the global optical depth, linear within each layer
  rtepack::stokvec_matrix source(NL, 2);
  source          = rtepack::stokvec{};
  Numeric tau_top = 0.0;
  for (Index l = 0; l < NL and c.thermal; l++) {
    const Numeric B0 = planck(c.frequency, c.temperature[l]), B1 = planck(c.frequency, c.temperature[l + 1]);
    const Numeric slope = (B1 - B0) / (tau[l] - tau_top);
    source[l, 0]        = {B0 - slope * tau_top, 0.0, 0.0, 0.0};
    source[l, 1]        = {slope, 0.0, 0.0, 0.0};
    tau_top             = tau[l];
  }

  const rtepack::stokvec beam_stokes{c.direct_flux > 0 ? c.direct_flux / c.mu0 : 0.0, 0.0, 0.0, 0.0};
  vdisort::main_data     v(NQuad,
                           NF,
                           AscendingGrid{std::move(tau)},
                           std::move(omega),
                           std::move(P),
                           std::move(bottom),
                           std::move(top),
                           std::move(source),
                           c.fresnel ? vdisort::brdf::fresnel_fourier_modes(*c.fresnel, NF)
                                     : vdisort::brdf::lambertian_fourier_modes(c.albedo, NF),
                           c.mu0,
                           beam_stokes,
                           c.phi0,
                           std::move(B));

  const auto q    = rt3::get_quadrature(N, rt3::quadrature_type::double_gauss);
  Numeric    node = 0.0;
  for (Index i = 0; i < N; i++)
    node = std::max({node,
                     std::abs(v.mu()[i] - q.mu[i]),
                     std::abs(v.mu()[N + i] + q.mu[i]),
                     std::abs(v.weights()[i] - q.weights[i])});
  require(node < 1e-14,
          std::format("VDISORT's double-Gauss streams and weights must equal RT3's for nmu = {}; they differ by up "
                      "to {:.1e}",
                      N,
                      node));
  return v;
}

//////////////////////////////////////////////////////////////////////////////
// Comparison
//////////////////////////////////////////////////////////////////////////////

//! How the VDISORT azimuth phi is chosen for RT3's psi (only correct is right)
enum class azimuth_map {
  correct,    // phi = phi0 + psi
  sense,      // phi = phi0 - psi: the azimuth sense reversed
  beam_sign,  // phi = psi - phi0: the beam azimuth with the wrong sign
};

//! max |VDISORT - RT3| / max |I_RT3| per Stokes component, and for the fluxes / max |F_RT3|
struct deviation {
  Vector4 stokes{};
  Numeric flux{};
  Numeric beyond{};  // VDISORT's Stokes components from nstokes on, which must be 0

  Numeric max() const { return std::max({stokes[0], stokes[1], stokes[2], stokes[3], flux}); }
};

/* Every level (top, interfaces, bottom), every stream, both directions and
   the test azimuths, for the Stokes components RT3 computes, and the up-
   and downward fluxes of I (down including the direct beam) and of Q. */
deviation compare(const setup&              c,
                  const rt3::result&        r,
                  const vdisort::main_data& v,
                  azimuth_map               map = azimuth_map::correct) {
  const Index  N = c.nmu, NL = c.nlay(), ns = c.nstokes;
  const Vector psi = test_azimuths();
  const auto   up = rt3::azimuth_radiance(r.up, psi), down = rt3::azimuth_radiance(r.down, psi);

  Numeric scale = 0.0;
  for (const auto* t : {&up, &down})
    for (Index l = 0; l <= NL; l++)
      for (Index k = 0; k < isize(psi); k++)
        for (Index i = 0; i < N; i++) scale = std::max(scale, std::abs((*t)[l, k, i, 0]));

  deviation d;
  for (Index l = 0; l <= NL; l++) {
    const Numeric tau = l == 0 ? 0.0 : v.tau()[l - 1];
    for (Index k = 0; k < isize(psi); k++) {
      const Numeric   phi = map == azimuth_map::correct ? c.phi0 + psi[k]
                            : map == azimuth_map::sense ? c.phi0 - psi[k]
                                                        : psi[k] - c.phi0;
      vdisort::u_data u;
      v.u(u, tau, wrapped(phi));
      for (Index i = 0; i < N; i++) {
        for (Index s = 0; s < 4; s++) {
          if (s < ns) {
            d.stokes[s] = std::max({d.stokes[s],
                                    std::abs(u.intensities[i][s] - up[l, k, i, s]) / scale,
                                    std::abs(u.intensities[N + i][s] - down[l, k, i, s]) / scale});
          } else {
            d.beyond =
                std::max({d.beyond, std::abs(u.intensities[i][s]) / scale, std::abs(u.intensities[N + i][s]) / scale});
          }
        }
      }
    }
  }

  Numeric fscale = 0.0;
  for (Index l = 0; l <= NL; l++) fscale = std::max({fscale, std::abs(r.up_flux[l, 0]), std::abs(r.down_flux[l, 0])});
  for (Index l = 0; l <= NL; l++) {
    const Numeric      tau = l == 0 ? 0.0 : v.tau()[l - 1];
    vdisort::flux_data fd;
    const auto         f = v.flux(fd, tau);
    d.flux               = std::max({d.flux,
                                     std::abs(f.up - r.up_flux[l, 0]) / fscale,
                                     std::abs(f.down_diffuse + f.down_direct - r.down_flux[l, 0]) / fscale});
    if (ns > 1) {
      Numeric qup = 0.0, qdown = 0.0;
      for (Index i = 0; i < N; i++) {
        qup   += 2 * pi * v.weights()[i] * v.mu()[i] * fd.u0[i][1];
        qdown += 2 * pi * v.weights()[i] * v.mu()[i] * fd.u0[N + i][1];
      }
      d.flux =
          std::max({d.flux, std::abs(qup - r.up_flux[l, 1]) / fscale, std::abs(qdown - r.down_flux[l, 1]) / fscale});
    }
  }
  return d;
}

//! max |Q|, |U|, |V| / max |I| of RT3 at the test azimuths, to show that each is exercised
Vector3 polarization(const rt3::result& r) {
  const Vector psi = test_azimuths();
  Vector4      top{};
  for (const auto* t : {&r.up, &r.down}) {
    const auto x = rt3::azimuth_radiance(*t, psi);
    for (Index l = 0; l < x.extent(0); l++)
      for (Index k = 0; k < x.extent(1); k++)
        for (Index i = 0; i < x.extent(2); i++)
          for (Index s = 0; s < x.extent(3); s++) top[s] = std::max(top[s], std::abs(x[l, k, i, s]));
  }
  return {top[1] / top[0], top[2] / top[0], top[3] / top[0]};
}

/* RT3's doubling, as in radtran3.f.  A layer with a scattering set (albedo
   > 0) is doubled n = int(F) + 1 times, with
   F = ln(max(tau, 1e-7) / max_delta_tau) / ln(2), or n = 0 if F <= 0.  The
   initial layer then has a vertical optical thickness delta = tau / 2^n.
   tau is the (delta-M scaled) particle plus gas extinction times the
   thickness. */
Numeric doublings(Numeric tau, Numeric max_delta_tau) {
  const Numeric F = std::log(std::max(tau, 1e-7) / max_delta_tau) / std::log(2.0);
  return F > 0 ? std::floor(F) + 1 : 0.0;
}

struct doubling {
  bool    any{false};
  Numeric delta{0.0};  // the largest initial-layer thickness of the doubled layers
  Numeric n{0.0};      // the largest number of doublings
};

doubling rt3_doubling(const setup& c, Numeric max_delta_tau) {
  doubling d;
  for (Index l = 0; l < c.nlay(); l++) {
    if (c.set_index[l] < 0) continue;
    const auto s = transport_set(c.sets[c.set_index[l]], c);
    if (s.scattering == 0.0) continue;
    const Numeric tau = (s.extinction + c.gas[l]) * c.dz(l);
    const Numeric n   = doublings(tau, max_delta_tau);
    d.any             = true;
    d.delta           = std::max(d.delta, tau / std::exp2(n));
    d.n               = std::max(d.n, n);
  }
  return d;
}

/* Tolerances, relative to max |I| (radiances) and max |F| (fluxes).  delta
   is the largest initial-layer thickness and n the largest number of
   doublings.
   - Direct, at max_delta_tau.  RT3's error is first order in delta.  The
     difference halves exactly when max_delta_tau halves (R7).  Its scale
     is the initial layer's slant thickness along the beam, delta / mu0,
     with a beam, and delta without one: RT3's initial-layer beam source
     (INITIAL_SOURCE) takes exp(-tau / mu0) at the top of the sublayer, a
     relative error of delta / (2 mu0).  Measured over all cases here and
     over mu0 = 0.2, 0.4 and 0.8 for R2, it is 0.41 to 0.68 delta / mu0
     with a beam, and 0.2 (R3) to 0.7 (R2) delta with thermal sources only
     (RT4: 0.1 to 0.7 delta).  It is printed as "c delta / mu0" or
     "c delta".  The tolerance is 10 times that scale plus RT3's doubling
     round-off, which grows like 2^n eps.
   - Richardson, 2 RT3(max_delta_tau / 2) - RT3(max_delta_tau) at
     max_delta_tau = 1e-5.  Every doubled layer gets exactly one more
     doubling (asserted), so this cancels the first-order term.  What
     remains is second order in delta plus three times the round-off.  The
     second-order term dominates from max_delta_tau = 1e-4 down to 1e-5
     (the residual / delta^2 is constant there) and grows with the number
     of streams like 1 / mu_min, the smallest stream cosine.  Measured, it
     is 0.04 to 0.12 delta^2 / mu_min for nmu = 8 and 16, with and without
     a beam, and at most 0.23 for nmu = 2 (R6).  At 3e-6 the round-off
     takes over.  The tolerance is delta^2 / mu_min + 3 2^n eps. */
constexpr Numeric direct_factor            = 10.0;
constexpr Numeric richardson_max_delta_tau = 1e-5;

//! The scale of RT3's first-order error: delta / mu0 with a beam, delta without one
Numeric error_scale(const setup& c, Numeric delta) { return c.direct_flux > 0 ? delta / std::min(c.mu0, 1.0) : delta; }

std::string_view error_scale_name(const setup& c) { return c.direct_flux > 0 ? "delta / mu0" : "delta"; }

Numeric direct_tolerance(const setup& c, const doubling& d) {
  return d.any ? direct_factor * error_scale(c, d.delta) + std::exp2(d.n) * eps : 1e-12;
}

Numeric richardson_tolerance(const setup& c) {
  const auto    coarse = rt3_doubling(c, richardson_max_delta_tau);
  const auto    fine   = rt3_doubling(c, richardson_max_delta_tau / 2);
  const Numeric mu_min = rt3::get_quadrature(c.nmu, rt3::quadrature_type::double_gauss).mu[0];
  return coarse.delta * coarse.delta / mu_min + 3 * std::exp2(fine.n) * eps;
}

//! RT3 extrapolated to max_delta_tau = 0 from richardson_max_delta_tau and half of it
rt3::result richardson(const setup& c) {
  for (Index l = 0; l < c.nlay(); l++) {
    if (c.set_index[l] < 0) continue;
    const Numeric tau = (transport_set(c.sets[c.set_index[l]], c).extinction + c.gas[l]) * c.dz(l);
    require(doublings(tau, richardson_max_delta_tau / 2) == doublings(tau, richardson_max_delta_tau) + 1,
            std::format("The Richardson extrapolation needs exactly one more doubling at max_delta_tau = {:.1e} "
                        "than at {:.1e}; a layer of optical thickness {} does not get it",
                        richardson_max_delta_tau / 2,
                        richardson_max_delta_tau,
                        tau));
  }
  auto a             = c;
  a.max_delta_tau    = richardson_max_delta_tau;
  const auto coarse  = rt3::solve(rt3_problem(a));
  a.max_delta_tau   /= 2;
  auto fine          = rt3::solve(rt3_problem(a));
  for (Index x = 0; x < static_cast<Index>(fine.up.size()); x++) {
    fine.up.data_handle()[x]   = 2 * fine.up.data_handle()[x] - coarse.up.data_handle()[x];
    fine.down.data_handle()[x] = 2 * fine.down.data_handle()[x] - coarse.down.data_handle()[x];
  }
  for (Index x = 0; x < static_cast<Index>(fine.up_flux.size()); x++) {
    fine.up_flux.data_handle()[x]   = 2 * fine.up_flux.data_handle()[x] - coarse.up_flux.data_handle()[x];
    fine.down_flux.data_handle()[x] = 2 * fine.down_flux.data_handle()[x] - coarse.down_flux.data_handle()[x];
  }
  return fine;
}

std::string format(const deviation& d) {
  return std::format("I {:9.3e}  Q {:9.3e}  U {:9.3e}  V {:9.3e}  F {:9.3e}",
                     d.stokes[0],
                     d.stokes[1],
                     d.stokes[2],
                     d.stokes[3],
                     d.flux);
}

void report(std::string_view what, const deviation& d, Numeric tol, std::string_view note = {}) {
  std::cout << std::format("{:<58} {}  tolerance {:.1e}{}\n", what, format(d), tol, note);
  require(d.max() <= tol,
          std::format("{}: VDISORT and RT3 must agree to max |difference| / max |I| (radiances) and / max |F| "
                      "(fluxes) <= {:.1e}; got {:.3e} for I, {:.3e} for Q, {:.3e} for U, {:.3e} for V and {:.3e} "
                      "for the fluxes",
                      what,
                      tol,
                      d.stokes[0],
                      d.stokes[1],
                      d.stokes[2],
                      d.stokes[3],
                      d.flux));
  require(d.beyond <= 1e-14,
          std::format("{}: VDISORT's Stokes components beyond RT3's nstokes must stay 0 (<= 1e-14 of max |I|), got "
                      "{:.3e}",
                      what,
                      d.beyond));
}

struct comparison {
  rt3::result        r;
  vdisort::main_data v;
};

comparison run(const setup& c) { return {.r = rt3::solve(rt3_problem(c)), .v = vdisort_solver(c)}; }

//! Compare at max_delta_tau and, if doubled, against RT3 Richardson-extrapolated
deviation check(std::string_view what, const setup& c, const comparison& x) {
  const auto dd = rt3_doubling(c, c.max_delta_tau);
  const auto d  = compare(c, x.r, x.v);
  const auto p  = polarization(x.r);
  report(what,
         d,
         direct_tolerance(c, dd),
         std::format("{}; |Q|,|U|,|V| / |I| = {:.1e}, {:.1e}, {:.1e}",
                     dd.any ? std::format(" = {:.2f} {}", d.max() / error_scale(c, dd.delta), error_scale_name(c))
                            : std::string{", no doubling"},
                     p[0],
                     p[1],
                     p[2]));
  if (dd.any)
    report("    RT3 Richardson-extrapolated to max_delta_tau = 0",
           compare(c, richardson(c), x.v),
           richardson_tolerance(c));
  return d;
}

deviation check(std::string_view what, const setup& c) {
  return check(std::format("{} (nmu {}, nstokes {})", what, c.nmu, c.nstokes), c, run(c));
}

//////////////////////////////////////////////////////////////////////////////
// Cases
//////////////////////////////////////////////////////////////////////////////

/* Rayleigh, m = 0, Q = I_v - I_h, mo outgoing, mi incident (signed):
     P_II = 3/8 (3 - mo^2 - mi^2 + 3 mo^2 mi^2),  P_IQ = 3/8 (1 - 3 mo^2)(1 - mi^2),
     P_QI = 3/8 (1 - mo^2)(1 - 3 mi^2),           P_QQ = 9/8 (1 - mo^2)(1 - mi^2),
     P_UU = 0,  P_VV = 3/2 mo mi, and no [I, Q] - [U, V] coupling. */
rtepack::muelmat rayleigh_m0(Numeric mo, Numeric mi) {
  const Numeric    a = mo * mo, b = mi * mi;
  rtepack::muelmat P{0.0};
  P[0, 0] = 3.0 / 8.0 * (3 - a - b + 3 * a * b);
  P[0, 1] = 3.0 / 8.0 * (1 - 3 * a) * (1 - b);
  P[1, 0] = 3.0 / 8.0 * (1 - a) * (1 - 3 * b);
  P[1, 1] = 9.0 / 8.0 * (1 - a) * (1 - b);
  P[3, 3] = 1.5 * mo * mi;
  return P;
}

/* Rayleigh (dipole) scattering of unpolarized light from k_in into k_out,
   independent of the rotation formulas of lab-frame.h: unpolarized light
   is the incoherent sum of two orthogonal linear polarizations e (the
   meridional e_v, e_h of k_in), the dipole radiates the part of e normal
   to k_out, and its components on the e_v, e_h of k_out are
   E_v = e_v . e, E_h = e_h . e, so
     I = sum (E_v^2 + E_h^2), Q = sum (E_v^2 - E_h^2), U = sum 2 E_v E_h, V = 0,
   scaled by 3/4 so that I = 3/4 (1 + cos^2 Theta). */
rtepack::stokvec rayleigh_column(Numeric mu_out, Numeric phi_out, Numeric mu_in, Numeric phi_in) {
  const auto basis = [](Numeric mu, Numeric phi) {
    const Numeric s = std::sqrt(1 - mu * mu);
    return std::pair<Vector3, Vector3>{{mu * std::cos(phi), mu * std::sin(phi), -s},
                                       {-std::sin(phi), std::cos(phi), 0.0}};
  };
  const auto [vi, hi] = basis(mu_in, phi_in);
  const auto [vo, ho] = basis(mu_out, phi_out);
  rtepack::stokvec z{};
  for (const auto& e : {vi, hi}) {
    const Numeric ev = dot(vo, e), eh = dot(ho, e);
    z[0] += 0.75 * (ev * ev + eh * eh);
    z[1] += 0.75 * (ev * ev - eh * eh);
    z[2] += 0.75 * 2.0 * ev * eh;
  }
  return z;
}

/* The single-scattering beam column sum_m eps_m (...) synthesised from the
   combined beam operator B [2, NF, 1, NQuad] as VDISORT's u() sums the
   field: I, Q from the cosine system with cos m(phi0 - phi) and the sine
   system with sin m(phi0 - phi); U, V the other way around. */
rtepack::stokvec synthesised(const vdisort::beam_phase_matrix_data& B, Index i, Numeric phi, Numeric phi0) {
  rtepack::stokvec z{};
  for (Index m = 0; m < B.extent(1); m++) {
    const Numeric e = m == 0 ? 1.0 : 2.0;
    const Numeric c = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const Numeric s = std::sin(static_cast<Numeric>(m) * (phi0 - phi));
    for (Index a = 0; a < 4; a++) {
      const Numeric cos_system  = B[vdisort::cosine_mode, m, 0, i][a, 0];
      const Numeric sin_system  = B[vdisort::sine_mode, m, 0, i][a, 0];
      z[a]                     += e * (a < 2 ? cos_system * c + sin_system * s : sin_system * c + cos_system * s);
    }
  }
  return z;
}

/* V0: the Fourier builder.
   (a) m = 0 Rayleigh, all 16 elements, against the closed form, on the
       8-node streams plus mu = 0.35 and 1, both hemispheres, at RT3's 16
       azimuth samples for degree 2.
   (b) The single-scattering beam column synthesised from combine_beam()
       against the independent dipole column, for mu0 = 0.6 and
       phi0 = 0 and 1.1, at the 16 streams and the test azimuths.  This
       pins C^m, S^m, the beam combination and VDISORT's field convention
       together.  vdisort::combine_beam_phase_matrices must pass the same
       check and equal combine_beam() exactly.
   (c) Evans' Mie series at RT3's 32 samples against 1024 midpoints:
       printed, it is the effect of the series' 1e-8 irregularity at
       Theta = 0 and 180 deg (F12 = 1e-8 there, F22 - F33 = 2e-8 forward). */
void test_fourier_builder() {
  const auto ray = legendre_matrix(rayleigh_legendre());
  {
    const Vector q = signed_streams(8);
    Vector       mu(isize(q) + 4);
    for (Index i = 0; i < isize(q); i++) mu[i] = q[i];
    mu[isize(q)]     = 0.35;
    mu[isize(q) + 1] = 1.0;
    mu[isize(q) + 2] = -0.35;
    mu[isize(q) + 3] = -1.0;
    const auto m     = fourier_modes(ray, mu, mu, 3, rt3_azimuth_samples(2, 2), 4);
    Numeric    diff  = 0.0;
    for (Index io = 0; io < isize(mu); io++)
      for (Index ii = 0; ii < isize(mu); ii++)
        for (Index a = 0; a < 4; a++)
          for (Index b = 0; b < 4; b++)
            diff = std::max(diff, std::abs(m.C[0, io, ii][a, b] - rayleigh_m0(mu[io], mu[ii])[a, b]));
    std::cout << std::format(
        "{:<58} max |C^0 - P_closed| {:9.3e}  tolerance 1.0e-13\n", "V0 (a) m = 0 Rayleigh, all 4 x 4 elements", diff);
    require(diff <= 1e-13,
            std::format("The m = 0 Fourier coefficient of the Rayleigh lab-frame matrix must equal the closed form to "
                        "1e-13 (16 azimuth samples), got {:.3e}",
                        diff));
  }

  {
    const Vector             mu  = signed_streams(8);
    const Vector             psi = test_azimuths();
    constexpr Numeric        mu0 = 0.6;
    const auto               b   = fourier_modes(ray, mu, Vector{-mu0}, 4, rt3_azimuth_samples(2, 3), 4);
    rtepack::muelmat_tensor3 C(4, 1, isize(mu), rtepack::muelmat{0.0}), S = C;
    for (Index m = 0; m < 4; m++) {
      for (Index i = 0; i < isize(mu); i++) {
        C[m, 0, i] = b.C[m, i, 0];
        S[m, 0, i] = b.S[m, i, 0];
      }
    }
    const auto derived = combine_beam(C, S), library = vdisort::combine_beam_phase_matrices(C, S);
    for (Numeric phi0 : {0.0, 1.1}) {
      Numeric scale = 0.0, d_derived = 0.0, d_library = 0.0, u_scale = 0.0;
      for (Index i = 0; i < isize(mu); i++) {
        for (Index k = 0; k < isize(psi); k++) {
          const Numeric phi = wrapped(phi0 + psi[k]);
          const auto    ref = rayleigh_column(mu[i], phi, -mu0, phi0);
          const auto    zd = synthesised(derived, i, phi, phi0), zl = synthesised(library, i, phi, phi0);
          scale   = std::max(scale, std::abs(ref[0]));
          u_scale = std::max(u_scale, std::abs(ref[2]));
          for (Index a = 0; a < 4; a++) {
            d_derived = std::max(d_derived, std::abs(zd[a] - ref[a]));
            d_library = std::max(d_library, std::abs(zl[a] - ref[a]));
          }
        }
      }
      std::cout << std::format(
          "{:<58} max |dev| / max |I| {:9.3e}  tolerance 1.0e-13; max |U| / max |I| {:.2f}\n"
          "    the same from vdisort::combine_beam_phase_matrices: {:.3e}  tolerance 1.0e-13\n",
          std::format("V0 (b) Rayleigh beam column vs dipole, phi0 = {}", phi0),
          d_derived / scale,
          u_scale / scale,
          d_library / scale);
      require(u_scale > 0.1 * scale, "V0 (b): the single-scattering U must be large enough to pin its sign");
      require(d_derived <= 1e-13 * scale,
              std::format("The single-scattering beam column synthesised from the combined beam operator must equal "
                          "the dipole column to 1e-13 of max I (mu0 = 0.6, phi0 = {}), got {:.3e}",
                          phi0,
                          d_derived / scale));
      require(d_library <= 1e-13 * scale,
              std::format("vdisort::combine_beam_phase_matrices must give the dipole column to 1e-13 of max I "
                          "(mu0 = 0.6, phi0 = {}), got {:.3e}",
                          phi0,
                          d_library / scale));
    }
    Numeric d_operator = 0.0;
    for (Index x = 0; x < static_cast<Index>(derived.size()); x++)
      for (Index a = 0; a < 4; a++)
        for (Index c = 0; c < 4; c++)
          d_operator = std::max(d_operator, std::abs(derived.data_handle()[x][a, c] - library.data_handle()[x][a, c]));
    require(d_operator == 0.0,
            std::format("vdisort::combine_beam_phase_matrices must equal the combination derived here exactly, they "
                        "differ by {:.3e}",
                        d_operator));
  }

  {
    const auto   mie = legendre_matrix(mie_legendre());
    const Vector mu  = signed_streams(8);
    const auto   a   = fourier_modes(mie, mu, mu, 12, rt3_azimuth_samples(11, 11), 4);
    const auto   b   = fourier_modes(mie, mu, mu, 12, 1024, 4, true);
    Numeric      d = 0.0, top = 0.0;
    for (Index x = 0; x < static_cast<Index>(a.C.size()); x++) {
      for (Index s = 0; s < 4; s++) {
        for (Index t = 0; t < 4; t++) {
          top = std::max(top, std::abs(b.C.data_handle()[x][s, t]));
          d   = std::max({d,
                          std::abs(a.C.data_handle()[x][s, t] - b.C.data_handle()[x][s, t]),
                          std::abs(a.S.data_handle()[x][s, t] - b.S.data_handle()[x][s, t])});
        }
      }
    }
    std::cout << std::format("{:<58} max |dC|, |dS| / max |C| {:9.3e}  (printed)\n",
                             "V0 (c) Mie C^m, S^m: RT3's 32 samples vs 1024 midpoints",
                             d / top);
  }
}

//! R1: a solar Rayleigh layer, tau 0.5, omega 0.95, mu0 0.6, Lambertian 0.1
setup rayleigh_layer() {
  setup c;
  c.nmu         = 8;
  c.aziorder    = 3;
  c.height      = Vector{1.0, 0.0};
  c.temperature = Vector{0.0, 0.0};
  c.gas         = Vector{0.0};
  c.set_index   = ArrayOfIndex{0};
  c.sets        = {{.extinction = 0.5, .scattering = 0.475, .legendre = rayleigh_legendre()}};
  c.direct_flux = 1.7;
  c.mu0         = 0.6;
  c.albedo      = 0.1;
  return c;
}

/* R1.  NFourier = 4: the m = 3 mode of Rayleigh scattering vanishes, in
   RT3's output and in VDISORT's input. */
void test_rayleigh() {
  const auto c = rayleigh_layer();
  const auto x = run(c);
  check(std::format("R1 Rayleigh tau 0.5 omega 0.95 mu0 0.6 A 0.1 (nmu {})", c.nmu), c, x);

  Numeric top = 0.0, m3 = 0.0;
  for (const auto* t : {&x.r.up, &x.r.down}) {
    for (Index l = 0; l < t->extent(0); l++) {
      for (Index i = 0; i < t->extent(2); i++) {
        for (Index s = 0; s < t->extent(3); s++) {
          top = std::max(top, std::abs((*t)[l, 0, i, s]));
          m3  = std::max(m3, std::abs((*t)[l, 3, i, s]));
        }
      }
    }
  }
  const auto F     = legendre_matrix(rayleigh_legendre());
  const auto modes = fourier_modes(F, signed_streams(c.nmu), signed_streams(c.nmu), 4, rt3_azimuth_samples(2, 3), 4);
  Numeric    in3 = 0.0, in0 = 0.0;
  for (Index i = 0; i < modes.C.extent(1); i++) {
    for (Index j = 0; j < modes.C.extent(2); j++) {
      for (Index a = 0; a < 4; a++) {
        for (Index b = 0; b < 4; b++) {
          in0 = std::max(in0, std::abs(modes.C[0, i, j][a, b]));
          in3 = std::max({in3, std::abs(modes.C[3, i, j][a, b]), std::abs(modes.S[3, i, j][a, b])});
        }
      }
    }
  }
  std::cout << std::format(
      "    m = 3: RT3 max |c_3| / max |c_0| {:.1e}; VDISORT input max |C^3|, |S^3| / max |C^0| {:.1e}\n",
      m3 / top,
      in3 / in0);
  require(m3 <= 1e-12 * top and in3 <= 1e-13 * in0,
          std::format("The m = 3 mode of Rayleigh scattering must vanish: RT3 output {:.1e} (<= 1e-12) and VDISORT "
                      "input {:.1e} (<= 1e-13), relative to m = 0",
                      m3 / top,
                      in3 / in0));
}

/* R2: Evans' mietest (runmietest, Evans and Stephens 1991): tau 1,
   omega 0.99, mu0 0.2, Lambertian 0.1, flux 0.2 pi, NFourier 12.  F34 makes
   V. */
setup mie_layer(Index nmu) {
  setup c;
  c.nmu         = nmu;
  c.aziorder    = 11;
  c.height      = Vector{1.0, 0.0};
  c.temperature = Vector{0.0, 0.0};
  c.gas         = Vector{0.0};
  c.set_index   = ArrayOfIndex{0};
  c.sets        = {{.extinction = 1.0, .scattering = 0.99, .legendre = mie_legendre()}};
  c.direct_flux = 0.2 * pi;
  c.mu0         = 0.2;
  c.albedo      = 0.1;
  return c;
}

void test_mie() {
  for (Index nmu : {8, 12}) check("R2 Evans mietest tau 1 omega 0.99 mu0 0.2 A 0.1", mie_layer(nmu));
}

/* R3: solar and thermal at 3 um, top-down a Rayleigh layer (tau 0.5, gas
   0.2), a Mie layer (tau 1.05, omega_p 0.99, gas 0.15) and a gas-only
   layer (tau 0.2), 200 to 295 K, over a Lambertian A = 0.25 surface at
   300 K, under a 250 K sky, with a flux of 0.5 W m-2 um-1 (a tenth of
   runtesta's) at mu0 = 0.5, so that solar and thermal radiances are of the
   same order. */
constexpr Numeric thermal_wavelength = 3e-6;

setup multilayer(Index nmu) {
  setup c;
  c.nmu         = nmu;
  c.aziorder    = 11;
  c.frequency   = Constant::c / thermal_wavelength;
  c.height      = Vector{15.0, 5.0, 2.0, 0.0};
  c.temperature = Vector{200.0, 260.0, 285.0, 295.0};
  c.gas         = Vector{0.02, 0.05, 0.1};
  c.set_index   = ArrayOfIndex{0, 1, -1};
  c.sets        = {{.extinction = 0.05, .scattering = 0.05, .legendre = rayleigh_legendre()},
                   {.extinction = 0.35, .scattering = 0.3465, .legendre = mie_legendre()}};
  c.direct_flux = 0.5 * (1e6 * thermal_wavelength) / c.frequency;
  c.mu0         = 0.5;
  c.thermal     = true;
  c.sky         = 250.0;
  c.surface     = 300.0;
  c.albedo      = 0.25;
  return c;
}

void test_multilayer() {
  const auto c = multilayer(8);
  const auto x = run(c);
  check(std::format("R3 Rayleigh/Mie/gas, solar + thermal 3 um, A 0.25 (nmu {})", c.nmu), c, x);

  // The thermal source alone, which must be a sizeable part of the solar + thermal radiance
  auto t        = c;
  t.direct_flux = 0.0;
  const auto y  = run(t);
  check(std::format("R3 thermal source only (nmu {})", t.nmu), t, y);
  Numeric th = 0.0, all = 0.0;
  for (Index k = 0; k < static_cast<Index>(y.r.up.size()); k++) {
    th  = std::max(th, std::abs(y.r.up.data_handle()[k]));
    all = std::max(all, std::abs(x.r.up.data_handle()[k]));
  }
  std::cout << std::format("    thermal-only max |I_up| / solar + thermal max |I_up|: {:.2f}\n", th / all);
  require(th > 0.1 * all and th < 0.9 * all,
          std::format("R3: both sources must matter, i.e. the thermal-only max |I_up| must be 0.1 to 0.9 of the "
                      "solar + thermal one; got {:.2f}",
                      th / all));
}

/* R4: VDISORT's beam azimuth phi0 = 1.1 on R2.  VDISORT's radiance at
   phi0 + psi must be RT3's at psi.  Taking VDISORT at phi0 - psi (the
   azimuth sense reversed; it flips U and V) or at psi - phi0 (phi0 with
   the wrong sign) must miss RT3 by more than 100 times the tolerance. */
void test_beam_azimuth() {
  auto c       = mie_layer(8);
  c.phi0       = 1.1;
  const auto x = run(c);
  check("R4 R2 with VDISORT phi0 = 1.1 rad (nmu 8)", c, x);
  const Numeric tol = direct_tolerance(c, rt3_doubling(c, c.max_delta_tau));
  for (auto [map, what] : {std::pair{azimuth_map::sense, "VDISORT at phi0 - psi (azimuth sense)"},
                           std::pair{azimuth_map::beam_sign, "VDISORT at psi - phi0 (phi0 sign)"}}) {
    const auto d = compare(c, x.r, x.v, map);
    std::cout << std::format(
        "    deliberately wrong: {:<40} {}  ({:.0f} x tolerance)\n", what, format(d), d.max() / tol);
    require(
        d.max() > 100 * tol,
        std::format("R4 must detect {} by more than 100 times the tolerance {:.1e}, got {:.3e}", what, tol, d.max()));
  }
}

//! R5: nstokes 1, 2 and 3 in RT3, against VDISORT with the leading nstokes x nstokes block of Z
void test_nstokes() {
  for (Index ns : {1, 2, 3}) {
    auto a    = mie_layer(8);
    a.nstokes = ns;
    check("R5 R2", a);
    auto b    = multilayer(8);
    b.nstokes = ns;
    check("R5 R3", b);
  }
}

/* R6: R2 and R3 for nmu = 2 ... 16 (RT3 needs nstokes nmu <= 64).  RT3's
   double-Gauss rule takes Legendre degrees up to 2 nmu - 3, so for nmu 2
   and 4 both solvers get the series truncated there, and NFourier =
   degree + 1.  The truncated matrices are not regular at Theta = 0 and
   180 deg, and with RT3's azimuth samples both still get the same discrete
   problem. */
void test_stream_sweep() {
  for (Index nmu : {2, 4, 8, 16}) {
    const Index degree = rt3::max_legendre_degree(nmu, rt3::quadrature_type::double_gauss);
    for (auto c : {mie_layer(nmu), multilayer(nmu)}) {
      Index top = 0;
      for (auto& s : c.sets) {
        s.legendre = truncated(s.legendre, degree);
        top        = std::max(top, stripped_degree(s.legendre));
      }
      c.aziorder = std::min(c.aziorder, top);
      check(std::format("R6 {} stream sweep, degree <= {}", c.nlay() == 1 ? "R2" : "R3", degree), c);
    }
  }
}

/* R7: the R3 difference must shrink with max_delta_tau like RT3's
   first-order doubling error.  The number of doublings is
   int(log2(tau / max_delta_tau)) + 1:
   - halving max_delta_tau adds exactly one doubling to every layer, so the
     difference must halve;
   - a decade changes delta by a factor of 8 or 16. */
void test_convergence() {
  auto                 c = multilayer(8);
  const auto           v = vdisort_solver(c);
  std::vector<Numeric> err;
  for (Numeric mdt : {1e-5, 5e-6, 1e-6, 1e-7}) {
    c.max_delta_tau = mdt;
    const auto dd   = rt3_doubling(c, mdt);
    const auto d    = compare(c, rt3::solve(rt3_problem(c)), v);
    err.push_back(d.max());
    report(std::format("R7 R3 convergence, max_delta_tau {:.0e}", mdt),
           d,
           direct_tolerance(c, dd),
           std::format(" = {:.2f} {}", d.max() / error_scale(c, dd.delta), error_scale_name(c)));
  }
  const Numeric halving = err[0] / err[1];
  std::cout << std::format("    difference ratio for max_delta_tau 1e-5 / 5e-6: {:.3f}\n", halving);
  require(halving > 1.9 and halving < 2.1,
          std::format("Halving max_delta_tau adds one doubling to every layer, so the VDISORT-RT3 difference must "
                      "halve (ratio 1.9 to 2.1); got {:.3f}",
                      halving));
  for (std::size_t i : {std::size_t{0}, std::size_t{2}}) {
    const Numeric ratio = err[i] / err[i + (i == 0 ? 2 : 1)];
    std::cout << std::format("    difference ratio per decade of max_delta_tau: {:.2f}\n", ratio);
    require(ratio > 4.0 and ratio < 25.0,
            std::format("The VDISORT-RT3 difference must fall like RT3's first-order doubling error, by a factor "
                        "4 to 25 per decade of max_delta_tau; got {:.2f}",
                        ratio));
  }
}

/* R8: the comparison is not blind.  Each deliberately wrong VDISORT input
   must miss RT3 (R2, nmu 8) by more than 100 times the tolerance.  A sign
   flip of S^m flips U and V only (it is the similarity transform
   diag(1, 1, -1, -1) of the combined systems), so U and V carry that
   detection.  Z at -phi' in the transform is the same error (it gives
   -S^m and the same C^m), and so is the reversed azimuth sense of R4. */
void test_mistakes() {
  const auto    c   = mie_layer(8);
  const auto    r   = rt3::solve(rt3_problem(c));
  const Numeric tol = direct_tolerance(c, rt3_doubling(c, c.max_delta_tau));
  for (auto m : {mistake::sine_sign,
                 mistake::mirrored_azimuth,
                 mistake::diffuse_epsilon,
                 mistake::beam_epsilon,
                 mistake::swapped_systems}) {
    const auto d = compare(c, r, vdisort_solver(c, m));
    std::cout << std::format(
        "    deliberately wrong VDISORT input: {:<40} {}  ({:.0f} x tolerance)\n", name(m), format(d), d.max() / tol);
    require(
        d.max() > 100 * tol,
        std::format("R8 must detect the VDISORT input error \"{}\" by more than 100 times the tolerance {:.1e}, got "
                    "{:.3e}",
                    name(m),
                    tol,
                    d.max()));
  }
}

/* R9: delta-M.  RT3 with delta_m against VDISORT given the identical
   scaled problem (transport_set()), at every level.  RT3's double-Gauss
   rule takes degrees up to 2 nmu - 3, while delta-M raises the degree to
   2 nmu - 1, so the scaled diagonal coefficients at 2 nmu - 2 and
   2 nmu - 1 must vanish.  The set is therefore (1 - f) Mie plus f times
   the forward peak truncated at degree M = 2 nmu, i.e. coefficients
   (1 - f) c_l + f (2 l + 1) on the diagonal, with the dyadic f = 1/4 so
   that c_l / (2 l + 1) = f exactly at l >= 12.  RT3 then finds f = 1/4,
   scales k by 1 - omega f, attenuates the beam with the scaled tau, and
   recovers the Mie series.  The Rayleigh set, of degree 2 < M, has f = 0. */
void test_delta_m() {
  auto          c   = multilayer(8);
  const Index   M   = 2 * c.nmu;
  const Numeric f   = 0.25;
  const Matrix  mie = mie_legendre();
  Matrix        peaked(M + 1, 6, 0.0);
  for (Index l = 0; l <= M; l++) {
    const auto w = static_cast<Numeric>(2 * l + 1);
    for (Index k = 0; k < 6; k++) {
      const Numeric x    = l < mie.nrows() ? mie[l, k] : 0.0;
      const bool    diag = k == 0 or k == 2 or k == 4 or k == 5;
      peaked[l, k]       = (1 - f) * x + (diag ? f * w : 0.0);
    }
  }
  c.sets[1].legendre = std::move(peaked);
  c.delta_m          = true;

  const auto s   = transport_set(c.sets[1], c);
  Numeric    dev = std::abs(s.extinction - (1 - 0.99 * f) * 0.35) / 0.35;
  for (Index l = 0; l < s.legendre.nrows(); l++)
    for (Index k = 0; k < 6; k++) dev = std::max(dev, std::abs(s.legendre[l, k] - (l < mie.nrows() ? mie[l, k] : 0.0)));
  std::cout << std::format("    R9 scaled set: f = {}, k' = {:.6f}, omega' = {:.6f}; max |scaled - Mie| {:.1e}\n",
                           f,
                           s.extinction,
                           s.scattering / s.extinction,
                           dev);
  require(dev < 1e-14, std::format("R9: RT3's delta-M algebra must recover the Mie series, got {:.1e}", dev));
  check("R9 R3 with a forward peak f = 1/4, RT3 delta-M", c);
}
//////////////////////////////////////////////////////////////////////////////
// Evans' benchmark settings (E)
//////////////////////////////////////////////////////////////////////////////

//! Exact-SI radiation constants in RT3's units, 2 h c^2 [W m-2 sr-1 um^4] and h c / k [um K]
constexpr Numeric planck_c1 = 2.0 * Constant::h * Constant::c * Constant::c * 1e24;
constexpr Numeric planck_c2 = Constant::h * Constant::c / Constant::k * 1e6;

/* A problem of Evans' scripts (3rdparty/polradtran/runmietest and
   runtesta), read from the script itself: its answers to rt3.f, its layer
   and scattering files and its expected output (the table).  The table is
   per micrometre and was made with Evans' 5-digit Planck constants;
   per_um converts it to W m-2 Hz-1 sr-1, and the temperatures are those at
   which the exact Planck function equals Evans' (see evans-scripts.h). */
struct evans_case {
  std::string             name;
  evans::rt3_settings     settings;
  std::vector<evans::row> table;
  setup                   c;                  // c.nmu is VDISORT's, to be set
  Numeric                 per_um;             // W m-2 Hz-1 sr-1 per W m-2 um-1 sr-1
  Vector                  mu;                 // Evans' streams, the table's |MU|
  bool                    brightness{false};  // the table holds V, H brightness temperatures (rt4.f, UNITS T)
};

evans_case read_evans(const std::string& script) {
  const auto sc = evans::read_script(std::filesystem::path(POLRADTRAN_DIR) / script);
  evans_case e;
  e.name  = script;
  e.table = evans::read_output(sc.files.at(sc.check));
  std::map<std::string, std::string> legendre_of;  // an RT4 scattering file made by scatcnv -> its Legendre input
  if (sc.solver().program == "rt3") {
    e.settings = evans::read_rt3_settings(sc);
  } else {
    // An RT4 script as the equivalent RT3 problem: thermal, m = 0, randomly oriented particles
    const auto r4 = evans::read_rt4_settings(sc);
    require(r4.units == 'T' and r4.polarization == "VH" and r4.nstokes == 2,
            std::format("{}: only EBB temperatures in V and H with nstokes 2 are handled", script));
    e.brightness = true;
    e.settings   = {.nstokes            = r4.nstokes,
                    .nmu                = r4.nmu,
                    .quad               = r4.quad,
                    .aziorder           = 0,
                    .layer_file         = r4.layer_file,
                    .delta_m            = false,
                    .src_code           = 2,
                    .ground_temperature = r4.ground_temperature,
                    .ground_type        = r4.ground_type,
                    .albedo             = r4.albedo,
                    .ground_index       = r4.ground_index,
                    .sky_temperature    = r4.sky_temperature,
                    .wavelength         = r4.wavelength,
                    .output             = r4.output};
    for (const auto& r : sc.runs)
      if (r.program == "scatcnv") {
        evans::answers a(r);
        const auto     in = a.first(), out = a.first();
        legendre_of[out] = in;
      }
  }
  const auto& st = e.settings;
  require(
      st.quad == 'G' and st.ground_type != 'S' and not st.delta_m,
      std::format("{}: only Gauss quadrature, a Lambertian or Fresnel ground and no delta-M are handled here", script));

  const auto levels = evans::read_layers(sc.files.at(st.layer_file));
  const auto T = [&](Numeric t) { return evans::exact_temperature_of_5digit(st.wavelength, t, planck_c1, planck_c2); };
  setup&     c = e.c;
  c.nstokes    = st.nstokes;
  c.aziorder   = st.aziorder;
  c.frequency  = Constant::c / (st.wavelength * 1e-6);
  e.per_um     = st.wavelength / c.frequency;
  const Index nlay = isize(levels) - 1;
  c.height         = Vector(nlay + 1);
  c.temperature    = Vector(nlay + 1);
  c.gas            = Vector(nlay);
  c.set_index      = ArrayOfIndex(nlay, -1);
  std::map<std::string, Index> set_of;
  for (Index l = 0; l <= nlay; l++) {
    c.height[l]      = levels[l].height;
    c.temperature[l] = T(levels[l].temperature);
    if (l == nlay) break;
    c.gas[l]         = levels[l].gas;
    const auto& file = levels[l].scattering_file;
    if (file.empty()) continue;
    if (not set_of.contains(file)) {
      require(sc.files.contains(legendre_of.contains(file) ? legendre_of[file] : file),
              std::format("{}: the scattering file {} is not a Legendre series (RT3 format)", script, file));
      const auto f = evans::read_scattering(sc.files.at(legendre_of.contains(file) ? legendre_of[file] : file));
      Matrix     L(isize(f.legendre), 6);
      for (Index i = 0; i < L.nrows(); i++)
        for (Index k = 0; k < 6; k++) L[i, k] = f.legendre[i][k];
      set_of[file] = isize(c.sets);
      c.sets.push_back({.extinction = f.extinction, .scattering = f.scattering, .legendre = std::move(L)});
    }
    c.set_index[l] = set_of[file];
  }
  c.direct_flux = st.direct_flux * e.per_um;
  c.mu0         = st.direct_mu;
  c.thermal     = st.src_code >= 2;
  c.sky         = T(st.sky_temperature);
  c.surface     = T(st.ground_temperature);
  c.albedo      = st.albedo;
  if (st.ground_type == 'F') c.fresnel = st.ground_index;
  e.mu = rt3::get_quadrature(st.nmu, rt3::quadrature_type::gauss).mu;
  return e;
}

//! A solution at the table's rows, per micrometre: [I, Q, U, V], fluxes (MU = -2, 2) in I and Q
using row_solution = std::function<Vector4(const evans::row&)>;

//! The table level of a height and Evans' stream of a |MU| (printed with 5 decimals)
Index level_of(const evans_case& e, Numeric z) {
  for (Index l = 0; l < isize(e.c.height); l++)
    if (std::abs(e.c.height[l] - z) < 1e-9) return l;
  throw std::runtime_error(std::format("{}: no level at Z = {}", e.name, z));
}

Index stream_of(const Vector& mu, Numeric m) {
  for (Index i = 0; i < isize(mu); i++)
    if (std::abs(mu[i] - std::abs(m)) < 6e-6) return i;
  throw std::runtime_error(std::format("No stream at MU = {}", m));
}

//! VDISORT at Evans' streams (user angles, formal solution) and fluxes, for a solved problem
row_solution vdisort_at_table(const vdisort::main_data& v, const evans_case& e, const setup& c) {
  const Index  N = c.nmu, NQuad = 2 * N, NL = c.nlay(), NF = c.nfourier(), ns = c.nstokes, n = isize(e.mu);
  const Vector mu = signed_streams(N);
  Vector       user(2 * n);  // up first, as VDISORT's streams
  for (Index i = 0; i < n; i++) {
    user[i]     = e.mu[i];
    user[n + i] = -e.mu[i];
  }

  rtepack::muelmat_tensor4 Cd(NF, NL, 2 * n, NQuad, rtepack::muelmat{0.0}), Sd = Cd;
  rtepack::muelmat_tensor3 Cb(NF, NL, 2 * n, rtepack::muelmat{0.0}), Sb        = Cb;
  for (Index l = 0; l < NL; l++) {
    const Index s = c.set_index[l];
    if (s < 0) continue;
    const auto  T    = transport_set(c.sets[s], c);
    const auto  F    = legendre_matrix(T.legendre);
    const Index nphi = rt3_azimuth_samples(T.legendre.nrows() - 1, c.aziorder);
    const auto  d    = fourier_modes(F, user, mu, NF, nphi, ns);
    const auto  b    = fourier_modes(F, user, Vector{-c.mu0}, NF, nphi, ns);
    for (Index m = 0; m < NF; m++) {
      for (Index u = 0; u < 2 * n; u++) {
        for (Index j = 0; j < NQuad; j++) {
          Cd[m, l, u, j] = d.C[m, u, j];
          Sd[m, l, u, j] = d.S[m, u, j];
        }
        Cb[m, l, u] = b.C[m, u, 0];
        Sb[m, l, u] = b.S[m, u, 0];
      }
    }
  }
  Vector levels(NL + 1, 0.0);
  for (Index l = 0; l < NL; l++) levels[l + 1] = v.tau()[l];
  const AscendingGrid tau{std::move(levels)};
  const Vector        psi{0.0, pi / 2, pi};  // the table's PHI = 0, 90, 180 deg
  Vector              phi(3);
  for (Index k = 0; k < 3; k++) phi[k] = wrapped(c.phi0 + psi[k]);
  // Every Evans angle in both directions: over a Fresnel surface the upward one reflects the downward one
  auto out = std::make_shared<rtepack::stokvec_tensor3>(NL + 1, 3, 2 * n);
  v.ungridded_u_user(
      *out, tau, phi, user, vdisort::combine_phase_matrices(Cd, Sd), vdisort::combine_beam_phase_matrices(Cb, Sb));

  auto flux = std::make_shared<Matrix>(4, NL + 1);
  v.ungridded_flux((*flux)[0], (*flux)[1], (*flux)[2], (*flux)[3], tau);

  return [out, flux, n, &e](const evans::row& r) -> Vector4 {
    const Index l = level_of(e, r.z);
    if (std::abs(r.mu) == 2.0) {
      const Numeric f = r.mu < 0 ? (*flux)[0, l] : (*flux)[1, l] + (*flux)[2, l];
      return {f / e.per_um, 0.0, 0.0, 0.0};
    }
    const Index k = static_cast<Index>(std::lround(r.phi / 90.0)), i = stream_of(e.mu, r.mu);
    const auto& x = (*out)[l, k, r.mu < 0 ? i : n + i];
    return {x[0] / e.per_um, x[1] / e.per_um, x[2] / e.per_um, x[3] / e.per_um};
  };
}

//! RT3 (the ARTS wrapper) on its own nmu nodes of the quadrature, which must be e.mu
row_solution rt3_at_table(const evans_case& e, Index nmu, rt3::quadrature_type quad, Numeric max_delta_tau) {
  setup c         = e.c;
  c.nmu           = nmu;
  c.max_delta_tau = max_delta_tau;
  auto p          = rt3_problem(c);
  p.quad          = quad;
  const auto   r  = std::make_shared<rt3::result>(rt3::solve(p));
  const Vector psi{0.0, pi / 2, pi};
  auto         up = std::make_shared<Tensor4>(rt3::azimuth_radiance(r->up, psi));
  auto         dn = std::make_shared<Tensor4>(rt3::azimuth_radiance(r->down, psi));
  require(isize(r->mu) == isize(e.mu), "RT3's streams must be the evaluation streams");
  for (Index i = 0; i < isize(e.mu); i++)
    require(std::abs(r->mu[i] - e.mu[i]) < 1e-15, "RT3's streams must be the evaluation streams");
  return [r, up, dn, &e](const evans::row& row) -> Vector4 {
    const Index l = level_of(e, row.z);
    if (std::abs(row.mu) == 2.0) {
      const auto& f = row.mu < 0 ? r->up_flux : r->down_flux;
      return {f[l, 0] / e.per_um, f.ncols() > 1 ? f[l, 1] / e.per_um : 0.0, 0.0, 0.0};
    }
    const Index k = static_cast<Index>(std::lround(row.phi / 90.0)), i = stream_of(e.mu, row.mu);
    const auto& t = row.mu < 0 ? *up : *dn;
    Vector4     x{};
    for (Index s = 0; s < t.extent(3); s++) x[s] = t[l, k, i, s] / e.per_um;
    return x;
  };
}

//! max |a - b| over the table's radiance rows per Stokes component / max I, and of the I fluxes / max flux
struct table_deviation {
  Vector4 stokes{};
  Numeric flux{};
};

//! The table's rows (fluxes and radiances at PHI = 0, 90, 180 deg and every level) on the streams mu
std::vector<evans::row> rows_on(const evans_case& e, const Vector& mu) {
  std::vector<evans::row> rows;
  for (Index l = 0; l < isize(e.c.height); l++) {
    rows.push_back({.z = e.c.height[l], .phi = 0.0, .mu = -2.0, .n = 4, .iquv = {}, .unit = {}});
    rows.push_back({.z = e.c.height[l], .phi = 0.0, .mu = 2.0, .n = 4, .iquv = {}, .unit = {}});
    for (Numeric phi : {0.0, 90.0, 180.0})
      for (Numeric sign : {-1.0, 1.0})
        for (Index i = 0; i < isize(mu); i++)
          rows.push_back({.z = e.c.height[l], .phi = phi, .mu = sign * mu[i], .n = 4, .iquv = {}, .unit = {}});
  }
  return rows;
}

//! max |a - b| over rows, per Stokes component / max I of b, and of the I fluxes / max flux of b
table_deviation deviation_from(const evans_case&              e,
                               const row_solution&            a,
                               const row_solution&            b,
                               const std::vector<evans::row>& rows) {
  Numeric max_i = 0.0, max_f = 0.0;
  for (const auto& r : rows) {
    const Numeric x = std::abs(b(r)[0]);
    if (std::abs(r.mu) == 2.0)
      max_f = std::max(max_f, x);
    else
      max_i = std::max(max_i, x);
  }
  table_deviation d;
  for (const auto& r : rows) {
    const auto x = a(r), y = b(r);
    if (std::abs(r.mu) == 2.0) {
      d.flux = std::max(d.flux, std::abs(x[0] - y[0]) / max_f);
    } else {
      for (Index s = 0; s < e.c.nstokes; s++) d.stokes[s] = std::max(d.stokes[s], std::abs(x[s] - y[s]) / max_i);
    }
  }
  return d;
}

//! Evans' table as a solution
row_solution table_of() {
  return [](const evans::row& r) -> Vector4 { return {r.iquv[0], r.iquv[1], r.iquv[2], r.iquv[3]}; };
}

Numeric max_of(const table_deviation& d) {
  return std::max({d.stokes[0], d.stokes[1], d.stokes[2], d.stokes[3], d.flux});
}

table_deviation print_deviation(std::string_view what, const table_deviation& d, Numeric tol = -1.0) {
  std::cout << std::format("    {:<62} I {:8.2e}  Q {:8.2e}  U {:8.2e}  V {:8.2e}  F {:8.2e}{}\n",
                           what,
                           d.stokes[0],
                           d.stokes[1],
                           d.stokes[2],
                           d.stokes[3],
                           d.flux,
                           tol < 0 ? "" : std::format("  tolerance {:.0e}", tol));
  require(tol < 0 or max_of(d) <= tol, std::format("{}: {:.2e} exceeds the tolerance {:.0e}", what, max_of(d), tol));
  return d;
}

/* E: VDISORT on the settings of Evans' scripts, runmietest (the Mie case
   of Evans and Stephens 1991) and runtesta (Rayleigh over Mie, gas, solar
   and thermal), read from the scripts.  Evans' tables are RT3 solutions with
   Gauss quadrature (the positive half of a 2 nmu-point Gauss-Legendre rule
   on [-1, 1]; nmu 8 and 4).  VDISORT has double-Gauss streams only, so it
   solves the same physical problem on its own streams and is evaluated at
   Evans' angles by its formal solution (ungridded_u_user).  The fluxes are
   compared in I.
   - RT3 (the ARTS wrapper) at Evans' settings must reproduce his tables to
     print precision.  That checks the reading of the scripts, the units
     (per micrometre, Evans' 5-digit Planck function) and mu0.
   - VDISORT must be converged in its streams at Evans' angles (nmu 16 and
     32 agree).
   - RT3 with double-Gauss nodes must agree with VDISORT on the same
     streams to RT3's doubling error: the same discrete problem (as R6).
   - Evans' tables differ from VDISORT by up to 1.3e-2 of max I, largest at
     grazing upwelling angles at the top.  That is the error of Gauss
     quadrature, which handles the discontinuity of the radiance at the
     horizon poorly: RT3 with Gauss quadrature approaches VDISORT steadily
     (about like 1 / nmu) as nmu grows from 4 to 16.  At Evans' nmu, RT3 is
     his table (above), so its distance from VDISORT is the table's. */
void test_evans_settings() {
  for (const std::string script : {"runmietest", "runtesta"}) {
    const auto e = read_evans(script);
    std::cout << std::format(
        "E {} (Evans' benchmark): {} layer(s), nstokes {}, Gauss nmu {}, aziorder {}, source code {}, "
        "mu0 {:.6f}\n",
        script,
        e.c.nlay(),
        e.settings.nstokes,
        e.settings.nmu,
        e.settings.aziorder,
        e.settings.src_code,
        e.c.mu0);
    const auto table = table_of();
    const auto evans = rt3_at_table(e, e.settings.nmu, rt3::quadrature_type::gauss, 1e-6);  // rt3.f's MAX_DELTA_TAU
    print_deviation("RT3 at Evans' settings vs his table", deviation_from(e, evans, table, e.table), 2e-6);

    setup c16 = e.c, c32 = e.c;
    c16.nmu        = 16;
    c32.nmu        = 32;
    const auto v16 = vdisort_solver(c16), v32 = vdisort_solver(c32);
    const auto v_at = [&](const evans_case& g) { return vdisort_at_table(v32, g, c32); };
    const auto v    = v_at(e);
    print_deviation(
        "VDISORT, nmu 16 vs 32, at Evans' angles", deviation_from(e, vdisort_at_table(v16, e, c16), v, e.table), 1e-6);
    print_deviation("VDISORT (nmu 32) at Evans' angles vs his table", deviation_from(e, v, table, e.table));

    // RT3 with Gauss quadrature approaches VDISORT (32 streams) at its nodes
    Numeric previous = 1.0;
    for (Index nmu : {4, 6, 8, 12, 16}) {
      evans_case g = e;
      g.mu         = rt3::get_quadrature(nmu, rt3::quadrature_type::gauss).mu;
      const auto d = print_deviation(
          std::format("RT3 Gauss nmu {:2} vs VDISORT, at RT3's nodes", nmu),
          deviation_from(g, rt3_at_table(g, nmu, rt3::quadrature_type::gauss, 1e-7), v_at(g), rows_on(g, g.mu)));
      require(max_of(d) < previous, "RT3 with Gauss quadrature must approach VDISORT as nmu grows");
      previous = max_of(d);
    }

    // RT3 and VDISORT on the same double-Gauss streams: the same discrete problem
    for (Index nmu : {8, 16}) {
      evans_case g = e;
      g.mu         = rt3::get_quadrature(nmu, rt3::quadrature_type::double_gauss).mu;
      setup c      = e.c;
      c.nmu        = nmu;
      print_deviation(std::format("RT3 double-Gauss nmu {:2} vs VDISORT on the same streams", nmu),
                      deviation_from(g,
                                     rt3_at_table(g, nmu, rt3::quadrature_type::double_gauss, 1e-7),
                                     vdisort_at_table(vdisort_solver(c), g, c),
                                     rows_on(g, g.mu)),
                      2e-6);
    }
  }
}

//! Evans' CONVERT_OUTPUT for UNITS 'T': [I, Q] per micrometre to the effective blackbody temperatures of V and H
std::array<Numeric, 2> brightness_vh(Numeric i, Numeric q, Numeric lambda, bool flux) {
  std::array<Numeric, 2> t{};
  for (Index k = 0; k < 2; k++) {
    Numeric rad = 2.0 * 0.5 * (i + (k == 0 ? q : -q));
    if (flux) rad /= pi;
    const Numeric sign = rad < 0 ? -1.0 : 1.0;
    t[k] = rad == 0 ? 0.0 : sign * 1.4388e4 / (lambda * std::log(1.0 + 1.1911e8 / (sign * rad * std::pow(lambda, 5))));
  }
  return t;
}

//! max |T_V, T_H of a - the table| over rows [K]
Numeric brightness_deviation(const evans_case& e, const row_solution& a, const std::vector<evans::row>& rows) {
  Numeric d = 0.0;
  for (const auto& r : rows) {
    const auto x = a(r);
    const auto t = brightness_vh(x[0], x[1], e.settings.wavelength, std::abs(r.mu) == 2.0);
    for (Index k = 0; k < 2; k++) d = std::max(d, std::abs(t[k] - r.iquv[k]));
  }
  return d;
}

/* E for Evans' RT4 script runtestr, read from the script: a 2 mm/h rain
   layer of spherical drops at 85 GHz over water (Fresnel, n = 3.17 -
   1.75i), 8 Gauss streams, I and Q, brightness temperatures in V and H.
   The drops are randomly oriented, so RT3 with aziorder 0 and VDISORT solve
   the same problem from the Mie Legendre series that Evans' scatcnv turns
   into RT4's scattering file.  For two Stokes components the Fresnel
   reflectivities depend only on |r_v|^2 and |r_h|^2, so the sign of the
   imaginary part of n does not matter.
   - RT3 at Evans' settings must reproduce his RT4 table to 0.01 K (its
     last printed digit): RT3's m = 0 doubling is RT4's.
   - VDISORT is evaluated at Evans' angles by its formal solution; over the
     Fresnel surface the upward radiance at mu reflects the downward one at
     -mu (BDRF::specular).  It must be converged in the streams (16 against
     32) to 0.01 K.  It differs from the table by the Gauss quadrature of
     the table, as in E above: RT3 with Gauss quadrature approaches VDISORT
     as nmu grows.
   - On the same double-Gauss streams VDISORT and RT3 solve the same
     discrete problem, with the reflection, and must agree to RT3's
     doubling error at every level, stream and direction. */
void test_evans_rt4_settings() {
  const auto e = read_evans("runtestr");
  std::cout << std::format(
      "E runtestr (Evans' RT4 benchmark): {} layer(s), nstokes {}, Gauss nmu {}, Fresnel n = {}{:+}i, "
      "brightness temperatures in V and H\n",
      e.c.nlay(),
      e.settings.nstokes,
      e.settings.nmu,
      e.c.fresnel->real(),
      e.c.fresnel->imag());
  const auto    evans = rt3_at_table(e, e.settings.nmu, rt3::quadrature_type::gauss, 1e-6);
  const Numeric rt3_k = brightness_deviation(e, evans, e.table);
  std::cout << std::format(
      "    {:<62} {:.4f} K (tolerance 0.01 K)\n", "RT3 at Evans' settings vs his RT4 table", rt3_k);
  require(rt3_k <= 0.01 + 1e-9, "RT3 at Evans' settings must reproduce his RT4 table to 0.01 K");

  // VDISORT at Evans' angles, 16 and 32 streams
  setup coarse = e.c, fine = e.c;
  coarse.nmu                       = 16;
  fine.nmu                         = 32;
  const auto              v_coarse = vdisort_solver(coarse), v_fine = vdisort_solver(fine);
  const auto              v = vdisort_at_table(v_fine, e, fine);
  std::vector<evans::row> radiances;
  for (const auto& r : e.table)
    if (std::abs(r.mu) != 2.0) radiances.push_back(r);
  const Numeric converged = [&] {
    const auto a = vdisort_at_table(v_coarse, e, coarse);
    Numeric    d = 0.0;
    for (const auto& r : radiances) {
      const auto x = brightness_vh(a(r)[0], a(r)[1], e.settings.wavelength, false);
      const auto y = brightness_vh(v(r)[0], v(r)[1], e.settings.wavelength, false);
      d            = std::max({d, std::abs(x[0] - y[0]), std::abs(x[1] - y[1])});
    }
    return d;
  }();
  std::cout << std::format(
      "    {:<62} {:.4f} K (tolerance 0.01 K)\n", "VDISORT, 16 vs 32 streams, at Evans' angles", converged);
  require(converged <= 0.01, "VDISORT's radiances at Evans' angles must be converged in the streams to 0.01 K");
  std::cout << std::format("    {:<62} {:.4f} K\n",
                           "VDISORT (32 streams) at Evans' angles vs his table",
                           brightness_deviation(e, v, radiances));

  // RT3 with Gauss quadrature approaches VDISORT at RT3's nodes
  Numeric previous = 1e9;
  for (Index nmu : {4, 8, 16}) {
    evans_case g  = e;
    g.mu          = rt3::get_quadrature(nmu, rt3::quadrature_type::gauss).mu;
    const auto r3 = rt3_at_table(g, nmu, rt3::quadrature_type::gauss, 1e-7);
    const auto vd = vdisort_at_table(v_fine, g, fine);
    Numeric    d  = 0.0;
    for (const auto& r : rows_on(g, g.mu)) {
      if (std::abs(r.mu) == 2.0 or r.phi != 0.0) continue;
      const auto a = brightness_vh(r3(r)[0], r3(r)[1], e.settings.wavelength, false);
      const auto b = brightness_vh(vd(r)[0], vd(r)[1], e.settings.wavelength, false);
      d            = std::max({d, std::abs(a[0] - b[0]), std::abs(a[1] - b[1])});
    }
    std::cout << std::format(
        "    {:<62} {:.4f} K\n", std::format("RT3 Gauss nmu {:2} vs VDISORT, at RT3's nodes", nmu), d);
    require(d < previous, "RT3 with Gauss quadrature must approach VDISORT as nmu grows");
    previous = d;
  }
  for (Index nmu : {8, 16}) {
    evans_case g    = e;
    g.mu            = rt3::get_quadrature(nmu, rt3::quadrature_type::double_gauss).mu;
    setup c         = e.c;
    c.nmu           = nmu;
    const auto rows = [&] {
      std::vector<evans::row> r;
      for (const auto& x : rows_on(g, g.mu))
        if (std::abs(x.mu) != 2.0 and x.phi == 0.0) r.push_back(x);  // m = 0: one azimuth
      return r;
    }();
    print_deviation(std::format("RT3 double-Gauss nmu {:2} vs VDISORT on the same streams", nmu),
                    deviation_from(g,
                                   rt3_at_table(g, nmu, rt3::quadrature_type::double_gauss, 1e-7),
                                   vdisort_at_table(vdisort_solver(c), g, c),
                                   rows),
                    2e-6);
  }
}

}  // namespace

int main() try {
  require(rt3::available(), "This test requires ENABLE_RT3=ON");
  test_fourier_builder();
  test_rayleigh();
  test_mie();
  test_multilayer();
  test_beam_azimuth();
  test_nstokes();
  test_stream_sweep();
  test_convergence();
  test_mistakes();
  test_delta_m();
  test_evans_settings();
  test_evans_rt4_settings();
  std::cout << "vdisort-rt3 comparison passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
