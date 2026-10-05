#include "rt4.h"

#include <arts_constants.h>
#include <debug.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <type_traits>

#ifdef ARTS_HAS_RT4
// ISO_C_BINDING entry points of 3rdparty/polradtran/rt4/rt4_c_interface.f90.
// The arrays are Fortran column-major; see the layouts in solve().
extern "C" {
void rt4_radtrano(std::int64_t nstokes,
                  std::int64_t nummu,
                  std::int64_t nuummu,
                  double       max_delta_tau,
                  char         quad_type,
                  double       ground_temp,
                  char         ground_type,
                  double       ground_albedo,
                  double       ground_index_re,
                  double       ground_index_im,
                  double*      ground_reflec,
                  double*      surf_reflect,
                  double*      gnd_radiance,
                  double       sky_temp,
                  double       wavelength,
                  std::int64_t num_layers,
                  double*      height,
                  double*      temperatures,
                  double*      gas_extinct,
                  std::int64_t nsl,
                  double*      scatlayers,
                  double*      extinct_matrix,
                  double*      emis_vector,
                  double*      scatter_matrix,
                  double*      mu_values,
                  double*      up_rad,
                  double*      down_rad);
void rt4_double_gauss_quadrature(std::int64_t num, double* abscissas, double* weights);
void rt4_gauss_legendre_quadrature(std::int64_t num, double* abscissas, double* weights);
void rt4_lobatto_quadrature(std::int64_t num, double* abscissas, double* weights);
}
#endif

namespace rt4 {
namespace {
#ifdef ARTS_HAS_RT4
//! RT4 uses COMMON blocks and large static local arrays: serialise every call.
std::mutex fortran_mutex;

//! Fixed sizes in radtran4.f (MAXV, MAXLAY, MAXLM).  MAXM = MAXV^2 and the
//! MINVERT limit of 256 are implied by MAXV.
constexpr Index max_vector       = 64;
constexpr Index max_layers       = 400;
constexpr Index max_layer_matrix = 301 * 4096;

char quad_char(quadrature_type type) {
  switch (type) {
    case quadrature_type::double_gauss: return 'D';
    case quadrature_type::gauss:        return 'G';
    case quadrature_type::lobatto:      return 'L';
  }
  ARTS_USER_ERROR("Unknown RT4 quadrature type {}", static_cast<int>(type));
}

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
  ARTS_USER_ERROR_IF(not available(), "RT4 requires ENABLE_RT4=ON");
  ARTS_USER_ERROR_IF(nmu < 1, "RT4 needs at least one quadrature node per hemisphere, got nmu = {}", nmu);

  quadrature q{.mu = Vector(nmu, 0.0), .weights = Vector(nmu, 0.0)};
#ifdef ARTS_HAS_RT4
  const char c = quad_char(type);

  std::lock_guard lock(fortran_mutex);
  if (c == 'D')
    rt4_double_gauss_quadrature(nmu, q.mu.data_handle(), q.weights.data_handle());
  else if (c == 'L')
    rt4_lobatto_quadrature(nmu, q.mu.data_handle(), q.weights.data_handle());
  else
    rt4_gauss_legendre_quadrature(nmu, q.mu.data_handle(), q.weights.data_handle());
#else
  (void)type;
#endif
  return q;
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
  }

  // RADTRANO arguments.  A row-major [a, b, c] array is the Fortran
  // column-major (c, b, a) array.  The legacy code declares no intent, so
  // the inputs are passed as copies.
  Vector height         = p.height;
  Vector temperature    = p.temperature;
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

  const Numeric wavelength_um = 1e6 * Constant::c / p.frequency;
  // RT4 radiance is per micrometre: B_nu = B_lambda[um^-1] * lambda[um] / f
  const Numeric per_um_to_per_hz = wavelength_um / p.frequency;

  // Surface.  GROUND_REFLEC(in, out) is read transposed by SPECULAR_SURFACE,
  // so the row-major R(out, in) is passed as is.  SURF_REFLECT(out s, out mu,
  // in s, in mu) is [in mu, in s, out mu, out s]; GND_RADIANCE(s, mu) is
  // [mu, s], input for 'A' and output otherwise.
  char    ground_type = 'L';
  Numeric albedo      = 0.0;
  Complex index{1.0, 0.0};
  Matrix  ground_reflec(ns, ns, 0.0);
  Tensor4 surf_reflect(nmu, ns, nmu, ns, 0.0);
  Matrix  gnd_radiance(nmu, ns, 0.0);
  std::visit(
      [&](const auto& g) {
        using T = std::remove_cvref_t<decltype(g)>;
        if constexpr (std::is_same_v<T, lambertian_surface>) {
          ground_type = 'L';
          albedo      = g.albedo;
        } else if constexpr (std::is_same_v<T, fresnel_surface>) {
          ground_type = 'F';
          index       = g.refractive_index;
        } else if constexpr (std::is_same_v<T, specular_surface>) {
          ARTS_USER_ERROR_IF(g.reflectivity.shape() != (std::array<Index, 2>{ns, ns}),
                             "specular_surface reflectivity must be [{}, {}] (nstokes = {}), got {:B,}",
                             ns,
                             ns,
                             ns,
                             g.reflectivity.shape());
          ground_type   = 'S';
          ground_reflec = g.reflectivity;
        } else {
          static_assert(std::is_same_v<T, discrete_surface>);
          ARTS_USER_ERROR_IF(g.reflection.shape() != (std::array<Index, 4>{nmu, nmu, ns, ns}) or
                                 g.emission.shape() != (std::array<Index, 2>{nmu, ns}),
                             "discrete_surface must have reflection [{}, {}, {}, {}] and emission [{}, {}] "
                             "(nmu_total = {}, nstokes = {}); got {:B,} and {:B,}",
                             nmu,
                             nmu,
                             ns,
                             ns,
                             nmu,
                             ns,
                             nmu,
                             ns,
                             g.reflection.shape(),
                             g.emission.shape());
          ground_type   = 'A';
          gnd_radiance  = g.emission;
          gnd_radiance /= per_um_to_per_hz;
          for (Index io = 0; io < nmu; io++)
            for (Index ii = 0; ii < nmu; ii++)
              for (Index so = 0; so < ns; so++)
                for (Index si = 0; si < ns; si++) surf_reflect[ii, si, io, so] = g.reflection[io, ii, so, si];
        }
      },
      p.ground);

  // MU_VALUES: RT4 writes the nquad nodes, the caller supplies the extra angles
  Vector mu(nmu, 0.0);
  mu[Range(nquad, nextra)] = p.extra_mu;

  // UP_RAD/DOWN_RAD(s, mu, level) is the row-major [level, mu, s] layout
  result r{.mu      = Vector(nmu, 0.0),
           .weights = Vector(nmu, 0.0),
           .up      = Tensor3(nlay + 1, nmu, ns, 0.0),
           .down    = Tensor3(nlay + 1, nmu, ns, 0.0)};

  {
    std::lock_guard lock(fortran_mutex);
    rt4_radtrano(ns,
                 nmu,
                 nextra,
                 p.max_delta_tau,
                 quad_char(p.quad),
                 p.surface_temperature,
                 ground_type,
                 albedo,
                 index.real(),
                 index.imag(),
                 ground_reflec.data_handle(),
                 surf_reflect.data_handle(),
                 gnd_radiance.data_handle(),
                 p.sky_temperature,
                 wavelength_um,
                 nlay,
                 height.data_handle(),
                 temperature.data_handle(),
                 gas_extinction.data_handle(),
                 nsl,
                 scatlayers.data_handle(),
                 extinct.data_handle(),
                 emis.data_handle(),
                 scatter.data_handle(),
                 mu.data_handle(),
                 r.up.data_handle(),
                 r.down.data_handle());
  }

  r.up   *= per_um_to_per_hz;
  r.down *= per_um_to_per_hz;
  r.mu    = mu;

  // RADTRANO does not return its weights; they are a function of the
  // quadrature alone, the extra angles having weight 0.
  const auto q               = get_quadrature(nquad, p.quad);
  r.weights[Range(0, nquad)] = q.weights;
  return r;
#else
  (void)p;
  return {};
#endif
}
}  // namespace rt4
