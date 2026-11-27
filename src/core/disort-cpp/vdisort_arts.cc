#include "vdisort_arts.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <legendre.h>
#include <phase_matrix.h>
#include <physics_funcs.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "common.h"
#include "vdisort-brdf.h"

namespace vdisort {
namespace {
/* The propagation zenith angles [deg] of the direction cosines mu (> 0
   upward), za = acos(mu), as a grid that is distinct and ascending, as
   ARTS's angular grids must be, and the grid index of each cosine. */
std::pair<Vector, ArrayOfIndex> zenith_angle_grid(const Vector& mu) {
  std::vector<Numeric> grid(mu.size());
  stdr::transform(mu, grid.begin(), [](Numeric x) { return Conversion::rad2deg(std::acos(x)); });
  stdr::sort(grid);
  grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
  const Vector za(std::move(grid));

  ArrayOfIndex index(mu.size());
  for (Size i = 0; i < mu.size(); i++)
    index[i] = static_cast<Index>(stdr::lower_bound(za, Conversion::rad2deg(std::acos(mu[i]))) - za.begin());
  return {za, index};
}

void check_cosines(const Vector& mu, std::string_view name) {
  if (stdr::any_of(mu, [](Numeric x) { return not(std::abs(x) <= 1.0 and x != 0.0); }))
    arts::user_error("The direction cosines {} must be in [-1, 0) or (0, 1]", name);
}

/* The expansion coefficients ("Greek constants") of the scattering matrix in generalized spherical functions,
   Lin et al. (2022, Eqs. 62-67), to the degree of its Legendre moments,
     F11 = sum_l a1_l P^{0,0}_l,  F22 + F33 = sum_l (a2 + a3)_l P^{2,2}_l,  F22 - F33 = sum_l (a2 - a3)_l P^{2,-2}_l,
     F44 = sum_l a4_l P^{0,0}_l,  F12 = sum_l b1_l P^{0,2}_l,      F34 = sum_l b2_l P^{0,2}_l,
   with a2, a3, b1 and b2 zero below degree 2.  The moments give F = sum_l (2 l + 1) beta_l P_l of degree
   L - 1 = moments.size() - 1, so a1_l = (2 l + 1) beta_l,11 and a4_l = (2 l + 1) beta_l,44, and the others are
   (2 l + 1) / 2 int F P^{m,n}_l dcos(Theta), whose integrands are polynomials of degree at most 2 L - 2 that the
   L-node Gauss-Legendre rule integrates exactly.  Every P^{m,n}_l of these is a polynomial of degree l, and the
   moments are those of the scattering matrix to degree L - 1, so the coefficients are those of the scattering
   matrix itself; a truncated Legendre series of each element is not regular at 0 and 180 deg, and this
   expansion is.  P^{m,n}_l are here Wigner's d^l_mn (Legendre::wigner_d); the convention of Hovenier et al. (2004)
   has the opposite sign of P^{0,2}_l, which cancels between b1, b2 and the functions they multiply. */
struct greek_constants {
  Vector a1, a2, a3, a4, b1, b2;

  explicit greek_constants(const CompactPlanarMuelmatConstVectorView& moments)
      : a1(moments.size(), 0.0),
        a2(moments.size(), 0.0),
        a3(moments.size(), 0.0),
        a4(moments.size(), 0.0),
        b1(moments.size(), 0.0),
        b2(moments.size(), 0.0) {
    const auto L = static_cast<Index>(moments.size());
    for (Index l = 0; l < L; l++) {
      a1[l] = static_cast<Numeric>(2 * l + 1) * moments[l].F11();
      a4[l] = static_cast<Numeric>(2 * l + 1) * moments[l].F44();
    }
    if (L < 3) return;

    Vector x(L), w(L), p(L), d22(L), d2m2(L), d02(L);
    Legendre::GaussLegendre(x, w);
    for (Index k = 0; k < L; k++) {
      Legendre::legendre_polynomials(p, x[k]);
      CompactPlanarMuelmat F{};
      for (Index l = 0; l < L; l++) F += (static_cast<Numeric>(2 * l + 1) * p[l]) * moments[l];
      const Numeric theta = std::acos(x[k]);
      Legendre::wigner_d(d22, 2, 2, theta);
      Legendre::wigner_d(d2m2, 2, -2, theta);
      Legendre::wigner_d(d02, 0, 2, theta);
      for (Index l = 2; l < L; l++) {
        const Numeric h  = 0.5 * static_cast<Numeric>(2 * l + 1) * w[k];
        a2[l]           += h * (0.5 * (F.F22() + F.F33()) * d22[l] + 0.5 * (F.F22() - F.F33()) * d2m2[l]);
        a3[l]           += h * (0.5 * (F.F22() + F.F33()) * d22[l] - 0.5 * (F.F22() - F.F33()) * d2m2[l]);
        b1[l]           += h * F.F12() * d02[l];
        b2[l]           += h * F.F34() * d02[l];
      }
    }
  }
};

/* The generalized spherical functions of Fourier order m in the matrix P^m_l(mu) of Lin et al. (2022, Eq. 59),
     P^m_l = [[P, 0, 0, 0], [0, R, -T, 0], [0, -T, R, 0], [0, 0, 0, P]],
   P = P^{m,0}_l, R = (P^{m,2}_l + P^{m,-2}_l) / 2 and T = (P^{m,2}_l - P^{m,-2}_l) / 2, [direction, l] */
struct order_functions {
  Matrix P, R, T;

  order_functions(Index m, const Vector& mu, Index L) : P(mu.size(), L), R(mu.size(), L), T(mu.size(), L) {
    Vector plus(L), minus(L);
    for (Size j = 0; j < mu.size(); j++) {
      const Numeric theta = std::acos(mu[j]);
      Legendre::wigner_d(P[j], m, 0, theta);
      Legendre::wigner_d(plus, m, 2, theta);
      Legendre::wigner_d(minus, m, -2, theta);
      for (Index l = 0; l < L; l++) {
        R[j, l] = 0.5 * (plus[l] + minus[l]);
        T[j, l] = 0.5 * (plus[l] - minus[l]);
      }
    }
  }
};
}  // namespace

fourier_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                 const AtmPoint&                 atm_point,
                                 const AtmSpectralPoint&         atm_spectral_point,
                                 Numeric                         frequency,
                                 const Vector&                   mu_out,
                                 const Vector&                   mu_in,
                                 Index                           nfourier,
                                 Numeric                         normalisation_tolerance) {
  const Index no = static_cast<Index>(mu_out.size()), ni = static_cast<Index>(mu_in.size());
  check_cosines(mu_out, "mu_out");
  check_cosines(mu_in, "mu_in");
  if (nfourier < 1) arts::user_error("nfourier must be >= 1, got {}", nfourier);
  if (not(normalisation_tolerance >= 0.0))
    arts::user_error("normalisation_tolerance must be >= 0, got {}", normalisation_tolerance);
  if (not(frequency > 0.0)) arts::user_error("frequency must be positive, got {} Hz", frequency);

  fourier_optics out{.extinction = 0.0,
                     .scattering = 0.0,
                     .cosine     = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0}),
                     .sine       = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0})};
  if (scattering_species.species.empty()) return out;

  if (no == 0 or ni == 0) arts::user_error("VDISORT's scattering optics need at least one mu_out and one mu_in");

  // The azimuthal Fourier modes Z = sum_m C_m cos(m phi) + S_m sin(m phi) of ARTS's laboratory-frame phase
  // matrix at the propagation zenith angles, with the extinction, the absorption and the phase integral
  // int Z11 dOmega at the incidence angles
  const auto [za_inc, inc]  = zenith_angle_grid(mu_in);
  const auto [za_scat, sca] = zenith_angle_grid(mu_out);
  const auto lab =
      scattering_species.aro_fourier(atm_point, atm_spectral_point, Vector{frequency}, za_inc, za_scat, nfourier - 1);
  if (not lab.phase_matrix.has_value())
    arts::user_error(
        "VDISORT needs the phase matrix of every scattering species; the bulk scattering properties "
        "have none");
  const auto& pha      = *lab.phase_matrix;
  const auto& ext      = lab.extinction_matrix;
  const auto& abs      = lab.absorption_vector;
  const auto& integral = pha.get_phase_integral();

  out.extinction      = ext[0, 0, 0].A();
  out.scattering      = out.extinction - abs[0, 0, 0].I();
  const Numeric sigma = integral[0, 0, 0];

  // VDISORT's optical depth and albedo are scalars: the particles' extinction must not depend on the direction
  // nor polarize, and their absorption and phase integral must not depend on the direction
  // An infinite normalisation_tolerance checks nothing; sigma must still be known to normalise P
  const bool    check = not std::isinf(normalisation_tolerance);
  const Numeric tol   = normalisation_tolerance * out.extinction;
  for (Index j = 0; j < static_cast<Index>(za_inc.size()); j++) {
    if (std::isnan(integral[0, 0, j]))
      arts::user_error(
          "The scattering species do not give the phase integral int Z11 dOmega at the incidence zenith "
          "angle {} deg (ARO data must cover all scattering zenith angles, [0, 180] deg)",
          za_inc[j]);
    if (not check) continue;
    const Numeric dk   = std::abs(ext[0, 0, j].A() - out.extinction),
                  dpol = std::max(std::abs(ext[0, 0, j].B()), std::abs(ext[0, 0, j].W()));
    const Numeric da = std::abs(abs[0, 0, j].I() - abs[0, 0, 0].I()), da2 = std::abs(abs[0, 0, j].Q());
    const Numeric ds = std::abs(integral[0, 0, j] - sigma);
    if (not(std::max({dk, dpol, da, da2, ds}) <= tol))
      arts::user_error(
          "VDISORT's optical depth and single-scattering albedo are scalars, but at the incidence zenith "
          "angle {} deg the particles' K11 differs by {} per m, K12 or K34 is {} per m, a1 differs by {} "
          "per m, a2 is {} per m, or the phase integral differs by {} per m from those at {} deg, more "
          "than normalisation_tolerance times K11, {} * {} per m (azimuthally randomly oriented "
          "particles generally give such optics, which VDISORT cannot transport)",
          za_inc[j],
          dk,
          dpol,
          da,
          da2,
          ds,
          za_inc[0],
          normalisation_tolerance,
          out.extinction);
  }

  if (check and not(std::abs(sigma - out.scattering) <= tol))
    arts::user_error(
        "The scattering coefficient from the phase matrix, int Z11 dOmega = {} per m, and the "
        "extinction minus the absorption, {} per m, must agree to normalisation_tolerance times the "
        "extinction, {} * {} per m (VDISORT normalises the phase matrix and takes the albedo from the "
        "latter)",
        sigma,
        out.scattering,
        normalisation_tolerance,
        out.extinction);
  if (sigma == 0.0) return out;

  // VDISORT's coefficients are (1 / 2 pi) int P cos(m phi) dphi with P = 4 pi Z / sigma, i.e. C_0 and half of
  // C_m and S_m for m > 0 (cpp.fast.vdisort-arts-test, V2)
  for (Index m = 0; m < nfourier; m++) {
    const Numeric scale = (m == 0 ? 1.0 : 0.5) * 4.0 * Constant::pi / sigma;
    for (Index o = 0; o < no; o++) {
      for (Index i = 0; i < ni; i++) {
        out.cosine[m, o, i] = scale * pha[0, 0, inc[i], sca[o], m, 0];
        out.sine[m, o, i]   = scale * pha[0, 0, inc[i], sca[o], m, 1];
      }
    }
  }
  return out;
}

fourier_optics legendre_fourier_optics(const CompactPlanarMuelmatConstVectorView& moments,
                                       const Vector&                              mu_out,
                                       const Vector&                              mu_in,
                                       Index                                      nfourier) {
  const Index no = static_cast<Index>(mu_out.size()), ni = static_cast<Index>(mu_in.size());
  const Index nleg = moments.size();
  check_cosines(mu_out, "mu_out");
  check_cosines(mu_in, "mu_in");
  if (nfourier < 1) arts::user_error("nfourier must be >= 1, got {}", nfourier);
  if (nleg < 1) arts::user_error("The Legendre moments need at least one degree");
  if (no == 0 or ni == 0) arts::user_error("The Fourier coefficients need at least one mu_out and one mu_in");

  fourier_optics out{.extinction = 0.0,
                     .scattering = 0.0,
                     .cosine     = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0}),
                     .sine       = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0})};

  /* The addition theorem of the generalized spherical functions (Lin et al. 2022, Eqs. 56-58),
       A^m = sum_{l=m}^{L-1} P^m_l(mu_out) Lambda_l P^m_l(mu_in),
       Lambda_l = [[a1, b1, 0, 0], [b1, a2, 0, 0], [0, 0, a3, b2], [0, 0, -b2, a4]]_l,
     gives VDISORT's coefficients (1 / 2 pi) int P cos(m phi) dphi and (1 / 2 pi) int P sin(m phi) dphi, phi the
     azimuth of incidence less that of scattering: C^m = (A^m + D A^m D) / 2, the elements of A^m in the I, Q and
     the U, V blocks, and S^m = (D A^m - A^m D) / 2, its other elements with the sign of their row in D,
     D = diag(1, 1, -1, -1). */
  const greek_constants g(moments);
  for (Index m = 0; m < std::min(nfourier, nleg); m++) {
    const order_functions fo(m, mu_out, nleg), fi(m, mu_in, nleg);
    for (Index o = 0; o < no; o++) {
      for (Index i = 0; i < ni; i++) {
        std::array<Numeric, 16> A{};
        for (Index l = m; l < nleg; l++) {
          const Numeric p_o = fo.P[o, l], R_o = fo.R[o, l], T_o = fo.T[o, l];
          const Numeric p_i = fi.P[i, l], R_i = fi.R[i, l], T_i = fi.T[i, l];
          const Numeric a1 = g.a1[l], a2 = g.a2[l], a3 = g.a3[l], a4 = g.a4[l], b1 = g.b1[l], b2 = g.b2[l];
          A[0]  += p_o * a1 * p_i;
          A[1]  += p_o * b1 * R_i;
          A[2]  -= p_o * b1 * T_i;
          A[4]  += R_o * b1 * p_i;
          A[5]  += R_o * a2 * R_i + T_o * a3 * T_i;
          A[6]  -= R_o * a2 * T_i + T_o * a3 * R_i;
          A[7]  -= T_o * b2 * p_i;
          A[8]  -= T_o * b1 * p_i;
          A[9]  -= T_o * a2 * R_i + R_o * a3 * T_i;
          A[10] += T_o * a2 * T_i + R_o * a3 * R_i;
          A[11] += R_o * b2 * p_i;
          A[13] += p_o * b2 * T_i;
          A[14] -= p_o * b2 * R_i;
          A[15] += p_o * a4 * p_i;
        }
        rtepack::muelmat C{0.0}, S{0.0};
        for (Index r = 0; r < 4; r++) {
          for (Index q = 0; q < 4; q++) {
            if ((r < 2) == (q < 2))
              C[r, q] = A[r * 4 + q];
            else
              S[r, q] = (r < 2 ? 1.0 : -1.0) * A[r * 4 + q];
          }
        }
        out.cosine[m, o, i] = C;
        out.sine[m, o, i]   = S;
      }
    }
  }
  return out;
}

main_data main_data_from_path(const ArrayOfPropagationPathPoint& ray_path,
                              const ArrayOfAtmPoint&             atm_path,
                              const AtmSpectralField&            atm_spectral_field,
                              const ArrayOfPropmatVector&        spectral_propmat_path,
                              const AscendingGrid&               freq_grid,
                              Index                              freq_index,
                              const ArrayOfScatteringSpecies&    scattering_species,
                              const path_settings&               settings,
                              const surface&                     ground,
                              Numeric                            surface_temperature,
                              Numeric                            sky_temperature) {
  const Index nlev = static_cast<Index>(ray_path.size());
  const Index nf   = static_cast<Index>(freq_grid.size());

  if (nlev < 2) arts::user_error("ray_path needs at least 2 points (1 layer), got {}", nlev);
  if (static_cast<Index>(atm_path.size()) != nlev or static_cast<Index>(spectral_propmat_path.size()) != nlev)
    arts::user_error(
        "ray_path, atm_path and spectral_propmat_path must have one entry per level; they have {}, {} "
        "and {}",
        nlev,
        atm_path.size(),
        spectral_propmat_path.size());
  if (freq_index < 0 or freq_index >= nf)
    arts::user_error("freq_index must be in [0, {}) for a freq_grid of {} frequencies, got {}", nf, nf, freq_index);
  if (stdr::any_of(spectral_propmat_path, [nf](const PropmatVector& v) { return static_cast<Index>(v.size()) != nf; }))
    arts::user_error("Every spectral_propmat_path level must have freq_grid.size() = {} propagation matrices", nf);
  for (Index l = 0; l < nlev - 1; l++) {
    if (not(ray_path[l].altitude() > ray_path[l + 1].altitude()))
      arts::user_error(
          "The ray_path altitudes must decrease strictly from the first point (top of the atmosphere) "
          "to the last (surface)");
  }
  if (stdr::any_of(spectral_propmat_path,
                   [freq_index](const PropmatVector& v) { return v[freq_index].is_polarized(); }))
    arts::user_error(
        "VDISORT's gas extinction is scalar: the gas propagation matrices in spectral_propmat_path must "
        "not be polarized (only A may be non-zero) at frequency index {}",
        freq_index);
  if (settings.nquad < 2 or settings.nquad % 2 != 0)
    arts::user_error("nquad must be a positive even number, got {}", settings.nquad);
  if (settings.nfourier < 1) arts::user_error("nfourier must be >= 1, got {}", settings.nfourier);
  if (not(settings.beam_flux >= 0.0))
    arts::user_error("beam_flux must be >= 0 (0 for no beam), got {}", settings.beam_flux);

  const bool    beam      = settings.beam_flux > 0.0;
  const Index   nlay      = nlev - 1;
  const Index   NQ        = settings.nquad;
  const Index   N         = NQ / 2;
  const Index   NF        = settings.nfourier;
  const Numeric frequency = freq_grid[freq_index];

  if (beam and not(settings.beam_mu > 0.0 and settings.beam_mu <= 1.0))
    arts::user_error("beam_mu must be in (0, 1] with a beam, got {}", settings.beam_mu);

  // VDISORT's own streams: upward (mu > 0) first, then mu = -mu[i]
  Vector mu(NQ), inv_mu(NQ), W(N);
  disort_common::initialize_streams(mu, inv_mu, W);

  // The incident directions: the streams, then the beam
  Vector mu_in(NQ + (beam ? 1 : 0));
  for (Index i = 0; i < NQ; i++) mu_in[i] = mu[i];
  if (beam) mu_in[NQ] = -settings.beam_mu;

  std::vector<fourier_optics> level;
  level.reserve(nlev);
  for (Index i = 0; i < nlev; i++) {
    const Vector3& pos = ray_path[i].pos;
    level.push_back(scattering_optics(scattering_species,
                                      atm_path[i],
                                      atm_spectral_field.at(AscendingGrid{frequency}, pos[0], pos[1], pos[2]),
                                      frequency,
                                      mu,
                                      mu_in,
                                      NF,
                                      settings.normalisation_tolerance));
  }

  Vector                   tau(nlay), omega(nlay);
  rtepack::muelmat_tensor4 C(NF, nlay, NQ, NQ, rtepack::muelmat{0.0}), S = C;
  rtepack::muelmat_tensor3 Cb(NF, nlay, NQ, rtepack::muelmat{0.0}), Sb   = Cb;
  Numeric                  t = 0.0;
  for (Index l = 0; l < nlay; l++) {
    const auto&   a = level[l];
    const auto&   b = level[l + 1];
    const Numeric gas =
        std::midpoint(spectral_propmat_path[l][freq_index].A(), spectral_propmat_path[l + 1][freq_index].A());
    const Numeric k  = gas + std::midpoint(a.extinction, b.extinction);
    const Numeric dz = ray_path[l].altitude() - ray_path[l + 1].altitude();
    if (not(k * dz > 0.0))
      arts::user_error(
          "VDISORT needs a positive optical thickness in every layer: the gas plus particle extinction "
          "must be positive along the whole path");
    t        += k * dz;
    tau[l]    = t;
    omega[l]  = std::midpoint(a.scattering, b.scattering) / k;

    const Numeric s = a.scattering + b.scattering;
    if (s == 0.0) continue;
    const Numeric wa = a.scattering / s, wb = b.scattering / s;
    for (Index m = 0; m < NF; m++) {
      for (Index o = 0; o < NQ; o++) {
        for (Index i = 0; i < NQ; i++) {
          C[m, l, o, i] = wa * a.cosine[m, o, i] + wb * b.cosine[m, o, i];
          S[m, l, o, i] = wa * a.sine[m, o, i] + wb * b.sine[m, o, i];
        }
        if (beam) {
          Cb[m, l, o] = wa * a.cosine[m, o, NQ] + wb * b.cosine[m, o, NQ];
          Sb[m, l, o] = wa * a.sine[m, o, NQ] + wb * b.sine[m, o, NQ];
        }
      }
    }
  }

  // B(tau) = c0 + c1 tau in the global optical depth, linear within each layer
  rtepack::stokvec_matrix source(nlay, 2);
  source = rtepack::stokvec{};
  if (settings.thermal) {
    Numeric tau_top = 0.0;
    for (Index l = 0; l < nlay; l++) {
      const Numeric B0 = planck(frequency, atm_path[l].temperature);
      const Numeric B1 = planck(frequency, atm_path[l + 1].temperature);
      const Numeric c1 = (B1 - B0) / (tau[l] - tau_top);
      source[l, 0]     = {B0 - c1 * tau_top, 0.0, 0.0, 0.0};
      source[l, 1]     = {c1, 0.0, 0.0, 0.0};
      tau_top          = tau[l];
    }
  }

  rtepack::stokvec_tensor3 bottom(2, NF, N), top(2, NF, N);
  bottom               = rtepack::stokvec{};
  top                  = rtepack::stokvec{};
  const Numeric B_sky  = planck(frequency, sky_temperature);
  const Numeric B_surf = settings.thermal ? planck(frequency, surface_temperature) : 0.0;
  // The surface's Planck radiance; VDISORT transmits it with 1 - R, Kirchhoff's emission
  for (Index i = 0; i < N; i++) {
    top[cosine_mode, 0, i]    = {B_sky, 0.0, 0.0, 0.0};
    bottom[cosine_mode, 0, i] = {B_surf, 0.0, 0.0, 0.0};
  }

  // The surface's reflection on the streams, and of the beam if there is one
  const ConstVectorView streams = mu[Range(0, N)];
  const Numeric         mu0     = beam ? settings.beam_mu : std::numeric_limits<Numeric>::quiet_NaN();
  surface_reflection    surface = std::visit(
      [&](const auto& g) {
        using T = std::remove_cvref_t<decltype(g)>;
        if constexpr (std::is_same_v<T, lambertian_surface>) {
          return brdf::lambertian_fourier_modes(g.albedo, NF, streams, streams, mu0);
        } else {
          static_assert(std::is_same_v<T, fresnel_surface>);
          return brdf::fresnel_fourier_modes(g.refractive_index, NF, streams, streams, mu0);
        }
      },
      ground);

  const rtepack::stokvec beam_stokes{beam ? settings.beam_flux / settings.beam_mu : 0.0, 0.0, 0.0, 0.0};
  return main_data(NQ,
                   NF,
                   AscendingGrid{std::move(tau)},
                   std::move(omega),
                   combine_phase_matrices(C, S),
                   std::move(bottom),
                   std::move(top),
                   std::move(source),
                   std::move(surface),
                   settings.beam_mu,
                   beam_stokes,
                   settings.beam_azimuth,
                   combine_beam_phase_matrices(Cb, Sb));
}
}  // namespace vdisort
