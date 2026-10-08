#include "rt3.h"

#include <arts_constants.h>
#include <debug.h>
#include <integration.h>

#include <algorithm>
#include <cmath>
#include <mutex>

#ifdef ARTS_HAS_RT3
#include "radtran3.h"
#include "radutil3.h"
#endif

namespace rt3 {
namespace {
#ifdef ARTS_HAS_RT3
//! Serialises every call.  The Fortran RADTRAN needed it (COMMON blocks,
//! SAVEd FFT tables and large static local arrays); its C++ port,
//! rt3::radtran, calls no Fortran and keeps no state, but the lock is still
//! taken.  RT4 has its own.
std::mutex fortran_mutex;

//! Fixed sizes in radtran3.f (MAXV, MAXLAY, MAXLM, MAXLEG, MAXSBUF, MAXDBUF)
constexpr Index max_vector       = 64;
constexpr Index max_layers       = 200;
constexpr Index max_layer_matrix = 101 * 4096;
constexpr Index max_coefficients = 1024;
constexpr Index max_scat_buffer  = 16 * max_layers * 2 * 4096;
constexpr Index max_direct_buf   = 16 * max_layers * 2 * max_vector;
//! The phase matrices are sampled at NUMPTS azimuths; with aziorder > 0 that
//! is 2 * 2^int(log2(degree + 4) + 1), which must fit FFT1DR's MAXN = 512
//! (and DIRECT_SCATTERING's 2 * MAXLEG = 512 samples): degree + 4 < 256.
constexpr Index max_fft_degree = 251;
//! Azimuth basis buffers: 2 * aziorder + 1 entries in DIRECT_SCATTERING's
//! BASIS_MATRIX (2 * 256) and in FOURIER_MATRIX's BASIS_VECTOR (4 * 256)
constexpr Index max_basis_direct = 512;
constexpr Index max_basis        = 1024;
//! RT3's CHECK_NORM (rt3::check_norm) throws above 1e-7; the discrete
//! normalisation equals legendre[0, 0] - 1 up to round-off for a series
//! within NLEGLIM.
constexpr Numeric normalisation_tolerance = 1e-9;

const char* quad_name(quadrature_type type) {
  switch (type) {
    case quadrature_type::gauss:        return "gauss";
    case quadrature_type::double_gauss: return "double_gauss";
    case quadrature_type::lobatto:      return "lobatto";
  }
  return "unknown";
}

//! The highest row of a [nleg + 1, 6] series with a non-zero coefficient (0 if none)
Index stripped_degree(const Matrix& legendre) {
  for (Index l = legendre.nrows() - 1; l > 0; l--)
    for (Index c = 0; c < 6; c++)
      if (legendre[l, c] != 0.0) return l;
  return 0;
}
#endif
}  // namespace

bool available() {
#ifdef ARTS_HAS_RT3
  return true;
#else
  return false;
#endif
}

quadrature get_quadrature(Index nmu, quadrature_type type) {
  ARTS_USER_ERROR_IF(nmu < 1, "RT3 needs at least one quadrature node per hemisphere, got nmu = {}", nmu);

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
  ARTS_USER_ERROR("Unknown RT3 quadrature type {}", static_cast<int>(type));
}

Index max_legendre_degree(Index nmu, quadrature_type type) {
  switch (type) {
    case quadrature_type::double_gauss: return std::max<Index>(2 * nmu - 3, 1);
    case quadrature_type::lobatto:      return std::max<Index>(4 * nmu - 5, 1);
    case quadrature_type::gauss:        return std::max<Index>(4 * nmu - 3, 1);
  }
  ARTS_USER_ERROR("Unknown RT3 quadrature type {}", static_cast<int>(type));
}

result solve(const problem& p) {
  ARTS_USER_ERROR_IF(not available(), "RT3 requires ENABLE_RT3=ON");

#ifdef ARTS_HAS_RT3
  const Index ns     = p.nstokes;
  const Index nquad  = p.nmu;
  const Index nextra = static_cast<Index>(p.extra_mu.size());
  const Index nmu    = nquad + nextra;
  const Index n      = ns * nmu;
  const Index nlay   = static_cast<Index>(p.height.size()) - 1;
  const Index nsl    = static_cast<Index>(p.scattering_sets.size());
  const Index nazi   = p.aziorder + 1;
  const bool  beam   = p.direct_flux > 0.0;

  ARTS_USER_ERROR_IF(ns < 1 or ns > 4, "RT3 supports nstokes 1 to 4, got {}", ns);
  ARTS_USER_ERROR_IF(nquad < 1, "RT3 needs at least one quadrature node per hemisphere, got nmu = {}", nquad);
  ARTS_USER_ERROR_IF(nextra > 0 and p.quad != quadrature_type::gauss,
                     "RT3 adds extra_mu angles only to the gauss quadrature (its 'E' type); got {} extra angles "
                     "with another quadrature",
                     nextra);
  ARTS_USER_ERROR_IF(stdr::any_of(p.extra_mu, [](Numeric mu) { return not(mu > 0.0 and mu <= 1.0); }),
                     "extra_mu values must be in (0, 1]");
  ARTS_USER_ERROR_IF(n > max_vector,
                     "RT3 requires nstokes * (nmu + extra_mu.size()) <= {}, got {} * ({} + {}) = {}",
                     max_vector,
                     ns,
                     nquad,
                     nextra,
                     n);
  ARTS_USER_ERROR_IF(p.aziorder < 0, "aziorder must be >= 0, got {}", p.aziorder);
  ARTS_USER_ERROR_IF(2 * p.aziorder + 1 > (beam ? max_basis_direct : max_basis),
                     "RT3 requires 2 * aziorder + 1 <= {} {}, got aziorder = {}",
                     beam ? max_basis_direct : max_basis,
                     beam ? "with a direct beam" : "without a direct beam",
                     p.aziorder);
  ARTS_USER_ERROR_IF(nlay < 1, "height needs at least 2 interfaces (1 layer), got {}", p.height.size());
  ARTS_USER_ERROR_IF(nlay > max_layers, "RT3 supports at most {} layers, got {}", max_layers, nlay);
  ARTS_USER_ERROR_IF((nlay + 1) * n * n > max_layer_matrix,
                     "RT3 requires (nlay + 1) * (nstokes * (nmu + extra_mu.size()))^2 <= {}, got ({} + 1) * {}^2 = {}",
                     max_layer_matrix,
                     nlay,
                     n,
                     (nlay + 1) * n * n);
  ARTS_USER_ERROR_IF(nsl > max_layers, "RT3 supports at most {} scattering sets, got {}", max_layers, nsl);
  ARTS_USER_ERROR_IF(nsl * nazi * 2 * n * n > max_scat_buffer,
                     "RT3 requires scattering_sets.size() * (aziorder + 1) * 2 * (nstokes * nmu_total)^2 <= {}, "
                     "got {} * {} * 2 * {}^2 = {}",
                     max_scat_buffer,
                     nsl,
                     nazi,
                     n,
                     nsl * nazi * 2 * n * n);
  ARTS_USER_ERROR_IF(beam and nazi * 2 * n * std::max(nlay, nsl) > max_direct_buf,
                     "With a direct beam RT3 requires (aziorder + 1) * 2 * nstokes * nmu_total * "
                     "max(nlay, scattering_sets.size()) <= {}, got {} * 2 * {} * {} = {}",
                     max_direct_buf,
                     nazi,
                     n,
                     std::max(nlay, nsl),
                     nazi * 2 * n * std::max(nlay, nsl));
  ARTS_USER_ERROR_IF(not(p.max_delta_tau > 0.0), "max_delta_tau must be positive, got {}", p.max_delta_tau);
  ARTS_USER_ERROR_IF(not(p.frequency > 0.0), "frequency must be positive, got {} Hz", p.frequency);
  ARTS_USER_ERROR_IF(not(p.direct_flux >= 0.0), "direct_flux must be >= 0 (0 for no beam), got {}", p.direct_flux);
  ARTS_USER_ERROR_IF(beam and not(p.direct_mu > 0.0 and p.direct_mu <= 1.0),
                     "direct_mu must be in (0, 1] with a direct beam, got {}",
                     p.direct_mu);
  ARTS_USER_ERROR_IF(beam and not std::holds_alternative<lambertian_surface>(p.ground),
                     "RT3 supports a direct beam only over a Lambertian surface");
  ARTS_USER_ERROR_IF(static_cast<Index>(p.temperature.size()) != nlay + 1 or
                         static_cast<Index>(p.gas_extinction.size()) != nlay or
                         static_cast<Index>(p.layer_scattering_index.size()) != nlay,
                     "With {} heights (nlay = {}), temperature needs nlay + 1 values and gas_extinction and "
                     "layer_scattering_index nlay values; got {}, {} and {}",
                     p.height.size(),
                     nlay,
                     p.temperature.size(),
                     p.gas_extinction.size(),
                     p.layer_scattering_index.size());
  ARTS_USER_ERROR_IF(p.thermal and stdr::any_of(p.temperature, [](Numeric t) { return not(t > 0.0); }),
                     "temperature values must be positive with thermal emission (RT3's linear-in-tau source of a "
                     "scattering layer degenerates to a constant when the Planck function at its top is 0)");
  ARTS_USER_ERROR_IF(stdr::any_of(p.gas_extinction, [](Numeric k) { return not(k >= 0.0); }),
                     "gas_extinction values must be non-negative (RT3 would silently clip them to 0)");
  ARTS_USER_ERROR_IF(stdr::any_of(p.layer_scattering_index, [nsl](Index i) { return i >= nsl; }),
                     "layer_scattering_index values must be < scattering_sets.size() = {} (negative means gas-only)",
                     nsl);

  // The series as RT3 uses it: trailing zero rows dropped, delta-M scaled
  // exactly as GET_SCAT_SET (READ_SCAT_FILE) does, then truncated to NLEGLIM.
  const Index  nleglim = max_legendre_degree(nquad, p.quad);
  const Index  mdm     = 2 * nmu;  // delta-M order M = 2 NUMMU, NUMMU including the extra angles
  ArrayOfIndex degree(nsl);
  for (Index iset = 0; iset < nsl; iset++) {
    const auto& s = p.scattering_sets[iset];
    ARTS_USER_ERROR_IF(s.legendre.ncols() != 6 or s.legendre.nrows() < 1,
                       "Scattering set {} must have legendre [nleg + 1, 6], got {:B,}",
                       iset,
                       s.legendre.shape());
    const Index nleg = stripped_degree(s.legendre);
    degree[iset]     = nleg;
    ARTS_USER_ERROR_IF(nleg + 1 > max_coefficients,
                       "RT3 holds at most {} Legendre coefficients per series; scattering set {} has degree {} "
                       "(after dropping trailing zero rows)",
                       max_coefficients,
                       iset,
                       nleg);

    Numeric f = 0.0;
    if (p.delta_m and mdm <= nleg) f = s.legendre[mdm, 0] / static_cast<Numeric>(2 * mdm + 1);
    ARTS_USER_ERROR_IF(p.delta_m and 1.0 - f == 0.0,
                       "Delta-M scaling of scattering set {} divides by 1 - f, with f = legendre[{}, 0] / {} = 1",
                       iset,
                       mdm,
                       2 * mdm + 1);
    // The coefficient c of row l and column k as RT3 sums it
    const auto scaled = [&](Index l, Index k) {
      const Numeric c = l <= nleg ? s.legendre[l, k] : 0.0;
      if (not p.delta_m) return c;
      const auto m    = static_cast<Numeric>(2 * l + 1);
      const bool diag = k == 0 or k == 2 or k == 4 or k == 5;
      return diag ? m * (c / m - f) / (1.0 - f) : m * (c / m) / (1.0 - f);
    };
    const Index rt3_degree = p.delta_m ? mdm - 1 : nleg;

    const Numeric c0 = scaled(0, 0);
    ARTS_USER_ERROR_IF(not(std::abs(c0 - 1.0) <= normalisation_tolerance),
                       "The phase function of scattering set {} must be normalised: legendre[0, 0] (F11, l = 0){} "
                       "must be 1 to {}, got {} (RT3's CHECK_NORM would reject it)",
                       iset,
                       p.delta_m ? " after delta-M scaling" : "",
                       normalisation_tolerance,
                       c0);

    bool dropped = false;
    for (Index l = nleglim + 1; l <= rt3_degree; l++)
      for (Index k = 0; k < 6; k++) dropped = dropped or scaled(l, k) != 0.0;
    ARTS_USER_ERROR_IF(dropped and p.delta_m,
                       "RT3 would silently truncate the delta-M scaled Legendre series of scattering set {} from "
                       "degree 2 * nmu_total - 1 = {} to {} (NLEGLIM of the {} quadrature with nmu = {}), dropping "
                       "non-zero coefficients.  With delta_m this always happens for double_gauss (nmu >= 2) and for "
                       "gauss with nmu or more extra angles; use gauss with fewer extra angles, lobatto, or no "
                       "delta_m",
                       iset,
                       rt3_degree,
                       nleglim,
                       quad_name(p.quad),
                       nquad);
    ARTS_USER_ERROR_IF(dropped,
                       "RT3 would silently truncate the Legendre series of scattering set {} from degree {} to {} "
                       "(NLEGLIM of the {} quadrature with nmu = {}), dropping non-zero coefficients.  Use more "
                       "quadrature nodes, a series of degree <= {}, or delta_m",
                       iset,
                       rt3_degree,
                       nleglim,
                       quad_name(p.quad),
                       nquad,
                       nleglim);

    const Index summed = std::min(rt3_degree, nleglim);
    ARTS_USER_ERROR_IF(p.aziorder > 0 and summed > max_fft_degree,
                       "With aziorder > 0 RT3 can sum Legendre series of degree <= {} (its FFT holds 512 azimuth "
                       "samples); scattering set {} is summed to degree {}",
                       max_fft_degree,
                       iset,
                       summed);
  }

  // RADTRAN arguments (radtran3.h).  A row-major [a, b, c] array is the
  // Fortran column-major (c, b, a) array.

  // SCATLAYERS(layer): 1-based set, 0 for gas-only
  ArrayOfIndex scatlayers(nlay);
  for (Index l = 0; l < nlay; l++)
    scatlayers[l] = p.layer_scattering_index[l] < 0 ? 0 : p.layer_scattering_index[l] + 1;

  // OUTLEVELS: every level, 1-based
  ArrayOfIndex outlevels(nlay + 1);
  for (Index l = 0; l <= nlay; l++) outlevels[l] = l + 1;

  // SCAT_COEF(6, LDCOEF, set) is [set, LDCOEF, 6]: the [nleg + 1, 6] legendre
  // of each set is its leading block
  const Index  ldcoef = nsl > 0 ? *stdr::max_element(degree) + 1 : 1;
  Vector       scat_extinct(nsl, 0.0), scat_scatter(nsl, 0.0);
  ArrayOfIndex scat_nlegen(nsl, 0);
  Tensor3      scat_coef(nsl, ldcoef, 6, 0.0);
  for (Index iset = 0; iset < nsl; iset++) {
    const auto& s                               = p.scattering_sets[iset];
    scat_extinct[iset]                          = s.extinction;
    scat_scatter[iset]                          = s.scattering;
    scat_nlegen[iset]                           = degree[iset];
    scat_coef[iset, Range(0, degree[iset] + 1)] = s.legendre[Range(0, degree[iset] + 1)];
  }

  // MU_VALUES and their weights, as RADTRAN makes them: the nquad nodes
  // followed by the extra angles, which have weight 0
  const auto q = get_quadrature(nquad, p.quad);
  Vector     mu(nmu), weights(nmu, 0.0);
  mu[Range{0, nquad}]      = q.mu;
  mu[Range{nquad, nextra}] = p.extra_mu;
  weights[Range{0, nquad}] = q.weights;

  // The ground as RADTRAN's input: SURF_REFLECT(out s, out mu, in s, in mu)
  // of each mode is [mode, in mu, in s, out mu, out s], GND_RADIANCE(s, mu)
  // [mode, mu, s]
  Tensor5 surf_reflect(nazi, nmu, ns, nmu, ns);
  Tensor3 gnd_radiance(nazi, nmu, ns), direct_reflect(nazi, nmu, ns);

  // UP_RAD/DOWN_RAD(s, mu, m + 1, level) is the row-major [level, m, mu, s]
  // layout, UP_FLUX/DOWN_FLUX(s, level) the row-major [level, s]
  result r{.mu        = Vector(nmu, 0.0),
           .weights   = Vector(nmu, 0.0),
           .up        = Tensor4(nlay + 1, nazi, nmu, ns, 0.0),
           .down      = Tensor4(nlay + 1, nazi, nmu, ns, 0.0),
           .up_flux   = Matrix(nlay + 1, ns, 0.0),
           .down_flux = Matrix(nlay + 1, ns, 0.0)};

  const std::int64_t src_code = (beam ? 1 : 0) + (p.thermal ? 2 : 0);
  ground_surface(
      p.ground, src_code, mu, weights, p.frequency, p.surface_temperature, surf_reflect, gnd_radiance, direct_reflect);
  {
    std::lock_guard lock(fortran_mutex);
    radtran(p.max_delta_tau,
            src_code,
            p.quad,
            p.delta_m,
            p.direct_flux,
            beam ? p.direct_mu : 1.0,
            surf_reflect,
            gnd_radiance,
            direct_reflect,
            p.sky_temperature,
            p.frequency,
            p.height,
            p.temperature,
            p.gas_extinction,
            scat_extinct,
            scat_scatter,
            scat_nlegen,
            scat_coef,
            scatlayers,
            outlevels,
            p.extra_mu,
            mu,
            r.up_flux,
            r.down_flux,
            r.up,
            r.down);
  }

  r.mu      = mu;
  r.weights = weights;
  return r;
#else
  (void)p;
  return {};
#endif
}

Tensor4 azimuth_radiance(const Tensor4& coefficients, const Vector& phi) {
  const auto [nlev, nmode, nmu, ns] = coefficients.shape();
  ARTS_USER_ERROR_IF(ns > 4, "coefficients must have at most 4 Stokes components, got {}", ns);
  const Index nphi = static_cast<Index>(phi.size());

  Tensor4 out(nlev, nphi, nmu, ns, 0.0);
  for (Index k = 0; k < nphi; k++) {
    for (Index m = 0; m < nmode; m++) {
      const Numeric c = std::cos(static_cast<Numeric>(m) * phi[k]);
      const Numeric s = std::sin(static_cast<Numeric>(m) * phi[k]);
      for (Index l = 0; l < nlev; l++)
        for (Index i = 0; i < nmu; i++)
          for (Index st = 0; st < ns; st++) out[l, k, i, st] += coefficients[l, m, i, st] * (st < 2 ? c : s);
    }
  }
  return out;
}
}  // namespace rt3
