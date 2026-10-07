#include "rt4.h"

#include <arts_constants.h>
#include <debug.h>
#include <integration.h>

#include <algorithm>
#include <cmath>

#ifdef ARTS_HAS_RT4
#include "radtran4.h"
#include "radutil4.h"
#endif

namespace rt4 {
namespace {
#ifdef ARTS_HAS_RT4
//! Fixed sizes in radtran4.f (MAXV, MAXLAY, MAXLM), kept by radtrano().
//! MAXM = MAXV^2 and the MINVERT limit of 256 are implied by MAXV.
constexpr Index max_vector       = 64;
constexpr Index max_layers       = 400;
constexpr Index max_layer_matrix = 301 * 4096;

void check_mirror_symmetry(const layer_optics& o, Index iset, Index nmu, Index ns) {
  constexpr Numeric rel = 1e-10;

  const auto mismatch = [](Numeric a, Numeric b, Numeric scale) { return std::abs(a - b) > rel * scale; };

  Numeric kscale = 0.0, ascale = 0.0, zscale = 0.0;
  for (auto x : o.extinction | by_elem) kscale = std::max(kscale, std::abs(x));
  for (auto x : o.absorption | by_elem) ascale = std::max(ascale, std::abs(x));
  for (auto x : o.phase | by_elem) zscale = std::max(zscale, std::abs(x));

  bool k_ok = true, a_ok = true, zt_ok = true, zr_ok = true;
  for (Index i = 0; i < nmu; i++) {
    for (Index s = 0; s < ns; s++) {
      a_ok = a_ok and not mismatch(o.absorption[down, i, s], o.absorption[up, i, s], ascale);
      for (Index t = 0; t < ns; t++)
        k_ok = k_ok and not mismatch(o.extinction[down, i, s, t], o.extinction[up, i, s, t], kscale);
    }
    for (Index j = 0; j < nmu; j++) {
      for (Index s = 0; s < ns; s++) {
        for (Index t = 0; t < ns; t++) {
          zt_ok = zt_ok and not mismatch(o.phase[down, down, i, j, s, t], o.phase[up, up, i, j, s, t], zscale);
          zr_ok = zr_ok and not mismatch(o.phase[down, up, i, j, s, t], o.phase[up, down, i, j, s, t], zscale);
        }
      }
    }
  }

  ARTS_USER_ERROR_IF(not(k_ok and a_ok and zt_ok and zr_ok),
                     "RT4 doubling is hard-coded for media that are mirror symmetric between the hemispheres "
                     "(SYMMETRIC in radtran4.f), to {} relative to the largest magnitude of each quantity.  Optics "
                     "set {} breaks this: "
                     "extinction[down] == extinction[up]: {}, absorption[down] == absorption[up]: {}, "
                     "phase[down, down] == phase[up, up]: {}, phase[down, up] == phase[up, down]: {}",
                     rel,
                     iset,
                     k_ok,
                     a_ok,
                     zt_ok,
                     zr_ok);
}

/* RT4 conserves energy only if every incident quadrature stream scatters
   K11 - a1 into the quadrature streams of both hemispheres, 2 pi sum_i w_i
   (phase[down, h, i, j] + phase[up, h, i, j])[I, I].  RT4's own check
   (CHECK_NORM in radscat4.f) only printed a warning and was disabled in
   ARTS 2, whose interface checked and renormalised the phase matrices
   instead.  Here it is an error.  If sampled phase matrices (forward peaks
   between the streams) make this a real problem, an explicit
   renormalisation step could be added, as ARTS 2 had. */
void check_normalisation(const layer_optics& o, Index iset, const quadrature& q, Index nquad, Numeric tolerance) {
  if (std::isinf(tolerance)) return;  // the user accepts any optics
  for (Index h = 0; h < 2; h++) {
    for (Index j = 0; j < nquad; j++) {
      Numeric scattered = 0.0;
      for (Index ho = 0; ho < 2; ho++)
        for (Index i = 0; i < nquad; i++) scattered += 2.0 * Constant::pi * q.weights[i] * o.phase[ho, h, i, j, 0, 0];
      const Numeric k11 = o.extinction[h, j, 0, 0], expected = k11 - o.absorption[h, j, 0];
      ARTS_USER_ERROR_IF(
          not(std::abs(scattered - expected) <= tolerance * std::abs(k11)),
          "RT4 optics set {} does not conserve energy on the streams: the {} quadrature stream {} (mu = {}) "
          "scatters 2 pi sum_i w_i phase[I, I] = {} into the quadrature streams, but K11 - a1 = {} - {} = {} (a "
          "difference of {:.3e} of K11, above normalisation_tolerance = {:.1e}).  Sample the phase matrix more "
          "finely or use more streams",
          iset,
          h == down ? "downward" : "upward",
          j,
          q.mu[j],
          scattered,
          k11,
          o.absorption[h, j, 0],
          expected,
          std::abs(scattered - expected) / std::abs(k11),
          tolerance);
    }
  }
}
#endif
}  // namespace

bool available() {
#ifdef ARTS_HAS_RT4
  return true;
#else
  return false;
#endif
}

quadrature get_quadrature(Index nmu, quadrature_type type) {
  ARTS_USER_ERROR_IF(nmu < 1, "RT4 needs at least one quadrature node per hemisphere, got nmu = {}", nmu);

  // The positive half of ARTS's 2 nmu-point rule on [-1, 1]
  const auto positive_half = [nmu](const auto& rule) {
    return quadrature{.mu      = Vector{rule.get_nodes()[Range{nmu, nmu}]},
                      .weights = Vector{rule.get_weights()[Range{nmu, nmu}]}};
  };
  switch (type) {
    case quadrature_type::double_gauss: return positive_half(scattering::DoubleGaussQuadrature(2 * nmu));
    case quadrature_type::gauss:        return positive_half(scattering::GaussLegendreQuadrature(2 * nmu));
    case quadrature_type::lobatto:      return positive_half(scattering::LobattoQuadrature(2 * nmu));
  }
  ARTS_USER_ERROR("Unknown RT4 quadrature type {}", static_cast<int>(type));
}

result solve(const problem& p) {
  ARTS_USER_ERROR_IF(not available(), "RT4 requires ENABLE_RT4=ON");

#ifdef ARTS_HAS_RT4
  const Index ns     = p.nstokes;
  const Index nquad  = p.nmu;
  const Index nextra = static_cast<Index>(p.extra_mu.size());
  const Index nmu    = nquad + nextra;
  const Index n      = ns * nmu;
  const Index nlay   = static_cast<Index>(p.height.size()) - 1;
  const Index nsl    = static_cast<Index>(p.optics.size());

  ARTS_USER_ERROR_IF(ns != 1 and ns != 2, "RT4 supports nstokes 1 ([I]) or 2 ([I, Q]), got {}", ns);
  ARTS_USER_ERROR_IF(nquad < 1, "RT4 needs at least one quadrature node per hemisphere, got nmu = {}", nquad);
  ARTS_USER_ERROR_IF(stdr::any_of(p.extra_mu, [](Numeric mu) { return not(mu > 0.0 and mu <= 1.0); }),
                     "extra_mu values must be in (0, 1]");
  ARTS_USER_ERROR_IF(n > max_vector,
                     "RT4 requires nstokes * (nmu + extra_mu.size()) <= {}, got {} * ({} + {}) = {}",
                     max_vector,
                     ns,
                     nquad,
                     nextra,
                     n);
  ARTS_USER_ERROR_IF(nlay < 1, "height needs at least 2 interfaces (1 layer), got {}", p.height.size());
  ARTS_USER_ERROR_IF(nlay > max_layers, "RT4 supports at most {} layers, got {}", max_layers, nlay);
  ARTS_USER_ERROR_IF((nlay + 1) * n * n > max_layer_matrix,
                     "RT4 requires (nlay + 1) * (nstokes * (nmu + extra_mu.size()))^2 <= {}, got ({} + 1) * {}^2 = {}",
                     max_layer_matrix,
                     nlay,
                     n,
                     (nlay + 1) * n * n);
  ARTS_USER_ERROR_IF(not(p.max_delta_tau > 0.0), "max_delta_tau must be positive, got {}", p.max_delta_tau);
  ARTS_USER_ERROR_IF(
      not(p.normalisation_tolerance >= 0.0), "normalisation_tolerance must be >= 0, got {}", p.normalisation_tolerance);
  ARTS_USER_ERROR_IF(not(p.frequency > 0.0), "frequency must be positive, got {} Hz", p.frequency);
  ARTS_USER_ERROR_IF(static_cast<Index>(p.temperature.size()) != nlay + 1 or
                         static_cast<Index>(p.gas_extinction.size()) != nlay or
                         static_cast<Index>(p.layer_optics_index.size()) != nlay,
                     "With {} heights (nlay = {}), temperature needs nlay + 1 values and gas_extinction and "
                     "layer_optics_index nlay values; got {}, {} and {}",
                     p.height.size(),
                     nlay,
                     p.temperature.size(),
                     p.gas_extinction.size(),
                     p.layer_optics_index.size());
  ARTS_USER_ERROR_IF(stdr::any_of(p.temperature, [](Numeric t) { return not(t > 0.0); }),
                     "temperature values must be positive (RT4's linear-in-tau source of a scattering layer "
                     "degenerates to a constant when the Planck function at its top is 0)");
  ARTS_USER_ERROR_IF(stdr::any_of(p.gas_extinction, [](Numeric k) { return not(k >= 0.0); }),
                     "gas_extinction values must be non-negative (RT4 would silently clip them to 0)");
  ARTS_USER_ERROR_IF(stdr::any_of(p.layer_optics_index, [nsl](Index i) { return i >= nsl; }),
                     "layer_optics_index values must be < optics.size() = {} (negative means gas-only)",
                     nsl);

  const auto quad_nodes = get_quadrature(nquad, p.quad);
  for (Index iset = 0; iset < nsl; iset++) {
    const auto& o = p.optics[iset];
    ARTS_USER_ERROR_IF(o.extinction.shape() != (std::array<Index, 4>{2, nmu, ns, ns}) or
                           o.absorption.shape() != (std::array<Index, 3>{2, nmu, ns}) or
                           o.phase.shape() != (std::array<Index, 6>{2, 2, nmu, nmu, ns, ns}),
                       "Optics set {} must have extinction [2, {}, {}, {}], absorption [2, {}, {}] and phase "
                       "[2, 2, {}, {}, {}, {}] (nmu_total = {}, nstokes = {}); got {:B,}, {:B,} and {:B,}",
                       iset,
                       nmu,
                       ns,
                       ns,
                       nmu,
                       ns,
                       nmu,
                       nmu,
                       ns,
                       ns,
                       nmu,
                       ns,
                       o.extinction.shape(),
                       o.absorption.shape(),
                       o.phase.shape());
    check_mirror_symmetry(o, iset, nmu, ns);
    check_normalisation(o, iset, quad_nodes, nquad, p.normalisation_tolerance);
  }

  // RADTRANO arguments (radtran4.h).  A row-major [a, b, c] array is the
  // Fortran column-major (c, b, a) array.  RADTRANO clips the gas
  // extinction at zero in place, so it gets a copy.
  Vector gas_extinction = p.gas_extinction;

  // SCATLAYERS(layer): 1-based optics set, 0 for gas-only
  Vector scatlayers(nlay);
  for (Index l = 0; l < nlay; l++)
    scatlayers[l] = p.layer_optics_index[l] < 0 ? 0.0 : static_cast<Numeric>(p.layer_optics_index[l] + 1);

  // EXTINCT_MATRIX(row, col, mu, hem, set) is [set, hem, mu, col, row],
  // EMIS_VECTOR(s, mu, hem, set) is [set, hem, mu, s] and
  // SCATTER_MATRIX(out s, out mu, in s, in mu, q, set) is
  // [set, q, in mu, in s, out mu, out s] with q = 2 * out_hem + in_hem
  // (0-based; RT4 q = 1 +<-+, 2 +<--, 3 -<-+, 4 -<--)
  const Index nset = std::max<Index>(nsl, 1);
  Tensor5     extinct(nset, 2, nmu, ns, ns, 0.0);
  Tensor4     emis(nset, 2, nmu, ns, 0.0);
  Tensor6     scatter(nset, 4, nmu, ns, nmu, ns, 0.0);
  for (Index iset = 0; iset < nsl; iset++) {
    const auto& o = p.optics[iset];
    emis[iset]    = o.absorption;
    for (Index h = 0; h < 2; h++)
      for (Index i = 0; i < nmu; i++) extinct[iset, h, i] = transpose(o.extinction[h, i]);
    for (Index ho = 0; ho < 2; ho++) {
      for (Index hi = 0; hi < 2; hi++) {
        const Index q = 2 * ho + hi;
        for (Index io = 0; io < nmu; io++)
          for (Index ii = 0; ii < nmu; ii++)
            for (Index so = 0; so < ns; so++)
              for (Index si = 0; si < ns; si++) scatter[iset, q, ii, si, io, so] = o.phase[ho, hi, io, ii, so, si];
      }
    }
  }

  // The streams, as RADTRANO makes them (MU_VALUES and its weights): the
  // nquad nodes followed by the extra angles, which have weight 0
  Vector mu(nmu), weights(nmu, 0.0);
  mu[Range{0, nquad}]      = quad_nodes.mu;
  mu[Range{nquad, nextra}] = p.extra_mu;
  weights[Range{0, nquad}] = quad_nodes.weights;

  // The ground as RADTRANO's external surface: SURF_REFLECT(out s, out mu,
  // in s, in mu) is [in mu, in s, out mu, out s], GND_RADIANCE(s, mu) is
  // [mu, s]
  Tensor4 surf_reflect(nmu, ns, nmu, ns);
  Matrix  gnd_radiance(nmu, ns);

  // UP_RAD/DOWN_RAD(s, mu, level) is the row-major [level, mu, s] layout
  result r{.mu      = Vector(nmu, 0.0),
           .weights = Vector(nmu, 0.0),
           .up      = Tensor3(nlay + 1, nmu, ns, 0.0),
           .down    = Tensor3(nlay + 1, nmu, ns, 0.0)};

  ground_surface(p.ground, mu, weights, p.frequency, p.surface_temperature, surf_reflect, gnd_radiance);

  rt4_workdata work;
  radtrano(p.max_delta_tau,
           p.quad,
           surf_reflect,
           gnd_radiance,
           p.sky_temperature,
           p.frequency,
           p.height,
           p.temperature,
           gas_extinction,
           scatlayers,
           extinct,
           emis,
           scatter,
           p.extra_mu,
           mu,
           r.up,
           r.down,
           work);

  r.mu      = mu;
  r.weights = weights;
  return r;
#else
  (void)p;
  return {};
#endif
}
}  // namespace rt4
