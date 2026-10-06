#include "vdisort_arts.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <physics_funcs.h>

#include <algorithm>
#include <cmath>
#include <memory>
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
  ARTS_USER_ERROR_IF(stdr::any_of(mu, [](Numeric x) { return not(std::abs(x) <= 1.0 and x != 0.0); }),
                     "The direction cosines {} must be in [-1, 0) or (0, 1]",
                     name);
}
}  // namespace

fourier_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                 const AtmPoint&                 atm_point,
                                 Numeric                         frequency,
                                 const Vector&                   mu_out,
                                 const Vector&                   mu_in,
                                 Index                           nfourier,
                                 Numeric                         normalisation_tolerance) {
  const Index no = static_cast<Index>(mu_out.size()), ni = static_cast<Index>(mu_in.size());
  check_cosines(mu_out, "mu_out");
  check_cosines(mu_in, "mu_in");
  ARTS_USER_ERROR_IF(nfourier < 1, "nfourier must be >= 1, got {}", nfourier);
  ARTS_USER_ERROR_IF(
      not(normalisation_tolerance >= 0.0), "normalisation_tolerance must be >= 0, got {}", normalisation_tolerance);
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "frequency must be positive, got {} Hz", frequency);

  fourier_optics out{.extinction = 0.0,
                     .scattering = 0.0,
                     .cosine     = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0}),
                     .sine       = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0})};
  if (scattering_species.species.empty()) return out;

  ARTS_USER_ERROR_IF(no == 0 or ni == 0, "VDISORT's scattering optics need at least one mu_out and one mu_in");

  // The azimuthal Fourier modes Z = sum_m C_m cos(m phi) + S_m sin(m phi) of ARTS's laboratory-frame phase
  // matrix at the propagation zenith angles, with the extinction, the absorption and the phase integral
  // int Z11 dOmega at the incidence angles
  const auto [za_inc, inc]  = zenith_angle_grid(mu_in);
  const auto [za_scat, sca] = zenith_angle_grid(mu_out);
  const auto lab            = scattering_species.get_bulk_scattering_properties_aro_fourier(
      atm_point, Vector{frequency}, za_inc, za_scat, nfourier - 1);
  ARTS_USER_ERROR_IF(not lab.phase_matrix.has_value(),
                     "VDISORT needs the phase matrix of every scattering species; the bulk scattering properties "
                     "have none");
  const auto& pha      = *lab.phase_matrix;
  const auto& ext      = lab.extinction_matrix;
  const auto& abs      = lab.absorption_vector;
  const auto& integral = pha.get_phase_integral();

  out.extinction      = ext[0, 0, 0, 0];
  out.scattering      = out.extinction - abs[0, 0, 0, 0];
  const Numeric sigma = integral[0, 0, 0];

  // VDISORT's optical depth and albedo are scalars: the particles' extinction must not depend on the direction
  // nor polarize, and their absorption and phase integral must not depend on the direction
  const Numeric tol = normalisation_tolerance * out.extinction;
  for (Index j = 0; j < static_cast<Index>(za_inc.size()); j++) {
    ARTS_USER_ERROR_IF(std::isnan(integral[0, 0, j]),
                       "The scattering species do not give the phase integral int Z11 dOmega at the incidence zenith "
                       "angle {} deg (ARO data must cover all scattering zenith angles, [0, 180] deg)",
                       za_inc[j])
    const Numeric dk = std::abs(ext[0, 0, j, 0] - out.extinction), dpol = std::max(std::abs(ext[0, 0, j, 1]),
                                                                                std::abs(ext[0, 0, j, 2]));
    const Numeric da = std::abs(abs[0, 0, j, 0] - abs[0, 0, 0, 0]), da2 = std::abs(abs[0, 0, j, 1]);
    const Numeric ds = std::abs(integral[0, 0, j] - sigma);
    ARTS_USER_ERROR_IF(not(std::max({dk, dpol, da, da2, ds}) <= tol),
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
                       out.extinction)
  }

  ARTS_USER_ERROR_IF(not(std::abs(sigma - out.scattering) <= tol),
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
        out.cosine[m, o, i] = scale * rtepack::muelmat{pha[0, 0, inc[i], sca[o], m, 0]};
        out.sine[m, o, i]   = scale * rtepack::muelmat{pha[0, 0, inc[i], sca[o], m, 1]};
      }
    }
  }
  return out;
}

main_data main_data_from_path(const ArrayOfPropagationPathPoint& ray_path,
                              const ArrayOfAtmPoint&             atm_path,
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

  ARTS_USER_ERROR_IF(nlev < 2, "ray_path needs at least 2 points (1 layer), got {}", nlev);
  ARTS_USER_ERROR_IF(
      static_cast<Index>(atm_path.size()) != nlev or static_cast<Index>(spectral_propmat_path.size()) != nlev,
      "ray_path, atm_path and spectral_propmat_path must have one entry per level; they have {}, {} "
      "and {}",
      nlev,
      atm_path.size(),
      spectral_propmat_path.size());
  ARTS_USER_ERROR_IF(freq_index < 0 or freq_index >= nf,
                     "freq_index must be in [0, {}) for a freq_grid of {} frequencies, got {}",
                     nf,
                     nf,
                     freq_index);
  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_propmat_path, [nf](const PropmatVector& v) { return static_cast<Index>(v.size()) != nf; }),
      "Every spectral_propmat_path level must have freq_grid.size() = {} propagation matrices",
      nf);
  for (Index l = 0; l < nlev - 1; l++)
    ARTS_USER_ERROR_IF(not(ray_path[l].altitude() > ray_path[l + 1].altitude()),
                       "The ray_path altitudes must decrease strictly from the first point (top of the atmosphere) "
                       "to the last (surface)");
  ARTS_USER_ERROR_IF(stdr::any_of(spectral_propmat_path,
                                  [freq_index](const PropmatVector& v) { return v[freq_index].is_polarized(); }),
                     "VDISORT's gas extinction is scalar: the gas propagation matrices in spectral_propmat_path must "
                     "not be polarized (only A may be non-zero) at frequency index {}",
                     freq_index);
  ARTS_USER_ERROR_IF(
      settings.nquad < 2 or settings.nquad % 2 != 0, "nquad must be a positive even number, got {}", settings.nquad);
  ARTS_USER_ERROR_IF(settings.nfourier < 1, "nfourier must be >= 1, got {}", settings.nfourier);
  ARTS_USER_ERROR_IF(
      not(settings.beam_flux >= 0.0), "beam_flux must be >= 0 (0 for no beam), got {}", settings.beam_flux);

  const bool    beam      = settings.beam_flux > 0.0;
  const Index   nlay      = nlev - 1;
  const Index   NQ        = settings.nquad;
  const Index   N         = NQ / 2;
  const Index   NF        = settings.nfourier;
  const Numeric frequency = freq_grid[freq_index];

  ARTS_USER_ERROR_IF(beam and not(settings.beam_mu > 0.0 and settings.beam_mu <= 1.0),
                     "beam_mu must be in (0, 1] with a beam, got {}",
                     settings.beam_mu);

  // VDISORT's own streams: upward (mu > 0) first, then mu = -mu[i]
  Vector mu(NQ), inv_mu(NQ), W(N);
  disort_common::initialize_streams(mu, inv_mu, W);

  // The incident directions: the streams, then the beam
  Vector mu_in(NQ + (beam ? 1 : 0));
  for (Index i = 0; i < NQ; i++) mu_in[i] = mu[i];
  if (beam) mu_in[NQ] = -settings.beam_mu;

  std::vector<fourier_optics> level;
  level.reserve(nlev);
  for (const auto& atm : atm_path)
    level.push_back(scattering_optics(scattering_species,
                                      atm,
                                      frequency,
                                      mu,
                                      mu_in,
                                      NF,
                                      settings.normalisation_tolerance));

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
    ARTS_USER_ERROR_IF(not(k * dz > 0.0),
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
  for (Index i = 0; i < N; i++) top[cosine_mode, 0, i] = {B_sky, 0.0, 0.0, 0.0};

  std::vector<BDRF> brdf = std::visit(
      [&](const auto& g) {
        using T = std::remove_cvref_t<decltype(g)>;
        if constexpr (std::is_same_v<T, lambertian_surface>) {
          for (Index i = 0; i < N; i++) bottom[cosine_mode, 0, i] = {(1.0 - g.albedo) * B_surf, 0.0, 0.0, 0.0};
          return brdf::lambertian_fourier_modes(g.albedo, NF);
        } else {
          static_assert(std::is_same_v<T, fresnel_surface>);
          const brdf::Fresnel fresnel{.refractive_index = g.refractive_index};
          for (Index i = 0; i < N; i++) {
            const auto R              = fresnel(mu[i]);
            bottom[cosine_mode, 0, i] = {
                (1.0 - R[0, 0]) * B_surf, -R[1, 0] * B_surf, -R[2, 0] * B_surf, -R[3, 0] * B_surf};
          }
          return brdf::fresnel_fourier_modes(g.refractive_index, NF);
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
                   std::move(brdf),
                   settings.beam_mu,
                   beam_stokes,
                   settings.beam_azimuth,
                   combine_beam_phase_matrices(Cb, Sb));
}
}  // namespace vdisort
