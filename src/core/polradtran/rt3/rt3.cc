#include "rt3.h"

#include <arts_constants.h>
#include <debug.h>

#include <algorithm>
#include <cmath>

#include "radtran3.h"
#include "radutil3.h"

namespace polradtran::rt3 {
namespace {
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
}  // namespace

Index max_legendre_degree(Index nmu, quadrature_type type) {
  switch (type) {
    case quadrature_type::double_gauss: return std::max<Index>(2 * nmu - 3, 1);
    case quadrature_type::lobatto:      return std::max<Index>(4 * nmu - 5, 1);
    case quadrature_type::gauss:        return std::max<Index>(4 * nmu - 3, 1);
  }
  ARTS_USER_ERROR("Unknown RT3 quadrature type {}", static_cast<int>(type));
}

result solve(const problem& p) {
  const Index ns     = p.nstokes;
  const Index nquad  = p.nmu;
  const Index nextra = static_cast<Index>(p.extra_mu.size());
  const Index nmu    = nquad + nextra;
  const Index nlay   = static_cast<Index>(p.height.size()) - 1;
  const Index nsl    = static_cast<Index>(p.scattering_sets.size());
  const Index nazi   = p.aziorder + 1;
  const bool  beam   = p.direct_flux > 0.0;

  ARTS_USER_ERROR_IF(ns < 1 or ns > 4, "RT3 supports nstokes 1 to 4, got {}", ns);
  ARTS_USER_ERROR_IF(nquad < 1, "RT3 needs at least one quadrature node per hemisphere, got nmu = {}", nquad);
  ARTS_USER_ERROR_IF(stdr::any_of(p.extra_mu, [](Numeric mu) { return not(mu > 0.0 and mu <= 1.0); }),
                     "extra_mu values must be in (0, 1]");
  ARTS_USER_ERROR_IF(p.aziorder < 0, "aziorder must be >= 0, got {}", p.aziorder);
  ARTS_USER_ERROR_IF(nlay < 1, "height needs at least 2 interfaces (1 layer), got {}", p.height.size());
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
    Numeric f        = 0.0;
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
  rt3_workdata work;
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
          r.down,
          work);

  r.mu      = mu;
  r.weights = weights;
  return r;
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
}  // namespace polradtran::rt3
