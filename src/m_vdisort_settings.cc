#include <array_algo.h>
#include <arts_omp.h>
#include <disort-brdf.h>
#include <legendre.h>
#include <matpack.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <sun_methods.h>
#include <vdisort-brdf.h>
#include <vdisort_settings.h>
#include <workspace.h>

#include <algorithm>
#include <concepts>
#include <functional>
#include <limits>
#include <numeric>
#include <string_view>
#include <type_traits>
#include <utility>

/* The VDISORT counterparts of the DISORT settings methods (m_disort_settings.cc), one for one.  They take the
   same inputs and do the same, with Stokes vectors where DISORT has scalars: an azimuth-independent Stokes
   vector [I, Q, U, V] is I and Q in the combined cosine system and U and V in the combined sine system, of mode
   0 (vdisort.h). */

namespace {
//! The azimuth-independent Stokes vector s as VDISORT's combined cosine and sine parts of mode 0
void set_mode_zero(StokvecTensor3View boundary, Index i, const Stokvec& s) {
  boundary[vdisort::cosine_mode, 0, i] = Stokvec{s.I(), s.Q(), 0.0, 0.0};
  boundary[vdisort::sine_mode, 0, i]   = Stokvec{0.0, 0.0, s.U(), s.V()};
}

void check_boundary(const StokvecTensor4& limit, Index nv, const char* which) {
  ARTS_USER_ERROR_IF(nv != limit.extent(0),
                     "Frequency grid size does not match the {} boundary condition size: {} vs {}",
                     which,
                     nv,
                     limit.extent(0))
  ARTS_USER_ERROR_IF(
      limit.extent(2) < 1, "Must have at least one fourier mode to use the {} boundary condition.", which)
}
}  // namespace

////////////////////////////////////////////////////////////////////////////////
// VDISORT settings initialization
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsInit(VDisortSettings&                   vdisort_settings,
                          const AscendingGrid&               freq_grid,
                          const ArrayOfPropagationPathPoint& ray_path,
                          const Index&                       quadrature_dimension,
                          const Index&                       legendre_polynomial_dimension,
                          const Index&                       fourier_mode_dimension) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(quadrature_dimension % 2, "Quadrature dimension ({}) must be even", quadrature_dimension);
  ARTS_USER_ERROR_IF(ray_path.size() < 2, "Must have at least one layer (two levels) to use the vdisort solver");

  vdisort_settings.resize(quadrature_dimension,
                          legendre_polynomial_dimension,
                          fourier_mode_dimension,
                          freq_grid,
                          DescendingGrid{ray_path.begin(), ray_path.end(), [](auto& v) { return v.altitude(); }});
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT solar source
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsNoSun(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  vdisort_settings.solar_source        = Stokvec{};
  vdisort_settings.solar_zenith_angle  = 0.0;
  vdisort_settings.solar_azimuth_angle = 0.0;
}

void vdisort_settingsSetSun(VDisortSettings&            vdisort_settings,
                            const AscendingGrid&        freq_grid,
                            const SurfaceField&         surf_field,
                            const Sun&                  sun,
                            const PropagationPathPoint& ray_point) {
  ARTS_TIME_REPORT

  const Numeric h = surf_field.single_value(SurfaceKey::h, sun.latitude, sun.longitude);

  const Vector3 sun_pos{sun.distance - h, sun.latitude, sun.longitude};
  const Vector2 los = geometric_los(ray_point.pos, sun_pos, surf_field.ellipsoid);

  const Numeric sin2_alpha = sun.sin_alpha_squared(ray_point.pos, surf_field.ellipsoid);

  const Size nv = freq_grid.size();

  ARTS_USER_ERROR_IF(vdisort_settings.solar_source.size() != nv or static_cast<Size>(sun.spectrum.nrows()) != nv or
                         sun.spectrum.ncols() != 4,
                     R"(Solar spectrum not agreeing with frequency grids:

freq_grid.size():                     {}
vdisort_settings.solar_source.size(): {}
sun.spectrum.shape():                 {:B,} (must be [freq_grid.size(), 4])
)",
                     nv,
                     vdisort_settings.solar_source.size(),
                     sun.spectrum.shape())

  for (Size iv = 0; iv < nv; iv++) {
    vdisort_settings.solar_source[iv] =
        sin2_alpha * Stokvec{sun.spectrum[iv, 0], sun.spectrum[iv, 1], sun.spectrum[iv, 2], sun.spectrum[iv, 3]};
  }

  vdisort_settings.solar_zenith_angle  = los[0];
  vdisort_settings.solar_azimuth_angle = los[1];
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT source polynomial
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsNoLayerThermalEmission(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  vdisort_settings.source_polynomial.resize(vdisort_settings.frequency_count(), vdisort_settings.layer_count(), 0);
}

void vdisort_settingsLayerThermalEmissionLinearInTau(VDisortSettings&       vdisort_settings,
                                                     const ArrayOfAtmPoint& atm_path,
                                                     const AscendingGrid&   freq_grid) {
  ARTS_TIME_REPORT

  const Size nv = freq_grid.size();
  const Size N  = atm_path.size();

  ARTS_USER_ERROR_IF(not same_shape({nv, N - 1}, vdisort_settings.optical_thicknesses),
                     "Incorrect shape: [{}, {}] vs {:B,}",
                     nv,
                     N - 1,
                     vdisort_settings.optical_thicknesses.shape());

  vdisort_settings.source_polynomial.resize(nv, N - 1, 2);

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const Numeric& f = freq_grid[iv];

    for (Size i = 0; i < N - 1; i++) {
      const Numeric y0 = planck(f, atm_path[i + 0].temperature);
      const Numeric y1 = planck(f, atm_path[i + 1].temperature);

      const Numeric x0 = i == 0 ? 0.0 : vdisort_settings.optical_thicknesses[iv, i - 1];
      const Numeric x1 = vdisort_settings.optical_thicknesses[iv, i];

      const Numeric b                              = (y1 - y0) / (x1 - x0);
      vdisort_settings.source_polynomial[iv, i, 0] = Stokvec{y0 - b * x0, 0.0, 0.0, 0.0};
      vdisort_settings.source_polynomial[iv, i, 1] = Stokvec{b, 0.0, 0.0, 0.0};
    }
  }
}

void vdisort_settingsLayerNonThermalEmissionLinearInTau(VDisortSettings&            vdisort_settings,
                                                        const ArrayOfAtmPoint&      atm_path,
                                                        const ArrayOfPropmatVector& spectral_propmat_path,
                                                        const ArrayOfStokvecVector& spectral_nlte_srcvec_path,
                                                        const AscendingGrid&        freq_grid) {
  ARTS_TIME_REPORT

  const Size nv = freq_grid.size();
  const Size N  = atm_path.size();

  ARTS_USER_ERROR_IF(not same_shape({nv, N - 1}, vdisort_settings.optical_thicknesses),
                     "Incorrect shape: [{}, {}] vs {:B,}",
                     nv,
                     N - 1,
                     vdisort_settings.optical_thicknesses.shape());

  vdisort_settings.source_polynomial.resize(nv, N - 1, 2);

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const Numeric& f = freq_grid[iv];

    for (Size i = 0; i < N - 1; i++) {
      Stokvec y0 = Stokvec{planck(f, atm_path[i + 0].temperature), 0.0, 0.0, 0.0} +
                   inv(spectral_propmat_path[i + 0][iv]) * spectral_nlte_srcvec_path[i + 0][iv];
      Stokvec y1 = Stokvec{planck(f, atm_path[i + 1].temperature), 0.0, 0.0, 0.0} +
                   inv(spectral_propmat_path[i + 1][iv]) * spectral_nlte_srcvec_path[i + 1][iv];

      // Numerically, these may be slightly downward even though they can mathematically only be exactly 0.
      if (y0[0] < 0) y0[0] = 0.0;
      if (y1[0] < 0) y1[0] = 0.0;

      const Numeric x0 = i == 0 ? 0.0 : vdisort_settings.optical_thicknesses[iv, i - 1];
      const Numeric x1 = vdisort_settings.optical_thicknesses[iv, i];

      const Stokvec b                              = (1.0 / (x1 - x0)) * (y1 - y0);
      vdisort_settings.source_polynomial[iv, i, 0] = y0 - x0 * b;
      vdisort_settings.source_polynomial[iv, i, 1] = b;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT upward boundary condition - from "surface" or "below"
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsNoSurfaceEmission(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  vdisort_settings.upward_boundary_condition = Stokvec{};
}

void vdisort_settingsSurfaceEmissionByTemperature(VDisortSettings&            vdisort_settings,
                                                  const AscendingGrid&        freq_grid,
                                                  const PropagationPathPoint& ray_point,
                                                  const SurfaceField&         surf_field) {
  ARTS_TIME_REPORT

  const Index   nv = freq_grid.size();
  const Numeric T  = surf_field.single_value(SurfaceKey::t, ray_point.latitude(), ray_point.longitude());

  auto& limit = vdisort_settings.upward_boundary_condition = Stokvec{};
  check_boundary(limit, nv, "upward");

  for (Index iv = 0; iv < nv; iv++) {
    const Stokvec B{planck(freq_grid[iv], T), 0.0, 0.0, 0.0};
    for (Index i = 0; i < limit.extent(3); i++) set_mode_zero(limit[iv], i, B);
  }
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT downward boundary condition - from "space" or "above"
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsNoSpaceEmission(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  vdisort_settings.downward_boundary_condition = Stokvec{};
}

void vdisort_settingsCosmicMicrowaveBackgroundRadiation(VDisortSettings&     vdisort_settings,
                                                        const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  const Index nv = freq_grid.size();

  auto& limit = vdisort_settings.downward_boundary_condition = Stokvec{};
  check_boundary(limit, nv, "downward");

  for (Index iv = 0; iv < nv; iv++) {
    const Stokvec B{planck(freq_grid[iv], Constant::cosmic_microwave_background_temperature), 0.0, 0.0, 0.0};
    for (Index i = 0; i < limit.extent(3); i++) set_mode_zero(limit[iv], i, B);
  }
}

void vdisort_settingsDownwellingObserver(const Workspace&                   ws,
                                         VDisortSettings&                   vdisort_settings,
                                         const AscendingGrid&               freq_grid,
                                         const ArrayOfPropagationPathPoint& ray_path,
                                         const AtmField&                    atm_field,
                                         const SurfaceField&                surf_field,
                                         const SubsurfaceField&             subsurf_field,
                                         const Agenda&                      spectral_rad_observer_agenda) {
  ARTS_TIME_REPORT

  const auto& ray_point = ray_path.front();

  auto& limit = vdisort_settings.downward_boundary_condition = Stokvec{};

  const Index nv = freq_grid.size();
  const Index N  = vdisort_settings.quadrature_dimension / 2;

  check_boundary(limit, nv, "downward");
  ARTS_USER_ERROR_IF(limit.extent(3) != N, "Must have at least one quadrature dimension.")

  Vector mu(N);
  Vector W(N);
  Legendre::PositiveDoubleGaussLegendre(mu, W);

  StokvecVector               spectral_rad;
  StokvecMatrix               spectral_rad_jac;
  ArrayOfPropagationPathPoint ray_path_up;
  const JacobianTargets       jac_targets{};

  String error{};

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(spectral_rad, spectral_rad_jac, ray_path_up)
  for (Index i = 0; i < N; i++) {
    try {
      spectral_rad_observer_agendaExecute(ws,
                                          spectral_rad,
                                          spectral_rad_jac,
                                          ray_path_up,
                                          freq_grid,
                                          jac_targets,
                                          ray_point.pos,
                                          {Conversion::acosd(mu[i]), ray_point.azimuth()},
                                          atm_field,
                                          surf_field,
                                          subsurf_field,
                                          spectral_rad_observer_agenda);

      for (Index iv = 0; iv < nv; iv++) set_mode_zero(limit[iv], i, spectral_rad[iv]);
    } catch (std::exception& e) {
#pragma omp critical
      error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), error)
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT Legendre coefficients
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsNoLegendre(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(vdisort_settings.legendre_coefficients.nrows() < 1,
                     "Must have at least one Legendre mode to use the Legendre coefficients.")

  vdisort_settings.legendre_coefficients                     = 0.0;
  vdisort_settings.legendre_coefficients[joker, joker, 0, 0] = 1.0;
}

void vdisort_settingsLegendreCoefficientsFromPath(
    VDisortSettings& vdisort_settings, const ArrayOfCompactPlanarMuelmatMatrix& spectral_phamat_spectral_path) try {
  ARTS_TIME_REPORT

  const Size  N = vdisort_settings.layer_count();
  const Index F = vdisort_settings.frequency_count();
  const Index L = spectral_phamat_spectral_path.front().ncols();

  ARTS_USER_ERROR_IF(
      (N + 1) != spectral_phamat_spectral_path.size(),
      R"(The number of levels in spectral_phamat_spectral_path must be one more than the number of layers in vdisort_settings.

  vdisort_settings.layer_count() + 1:   {}
  spectral_phamat_spectral_path.size(): {}
)",
      N + 1,
      spectral_phamat_spectral_path.size());

  ARTS_USER_ERROR_IF(
      L < vdisort_settings.legendre_polynomial_dimension or not all_same_shape({F, L}, spectral_phamat_spectral_path),
      "The shape of spectral_phamat_spectral_path must be {:B,}, at least one is not",
      std::array{F, L});

  vdisort_settings.legendre_coefficients.resize(F, N, L, 6);

  Vector invfac(L);
  for (Index j = 0; j < L; j++) { invfac[j] = 1.0 / std::sqrt(2 * j + 1); }

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      auto c = vdisort_settings.legendre_coefficients[iv, i];
      for (Index j = 0; j < L; j++) {
        // The six elements [F11, F12, F22, F33, F34, F44] of each coefficient on Y_l0
        for (Index e = 0; e < 6; e++) {
          c[j, e] = invfac[j] * std::midpoint(spectral_phamat_spectral_path[i][iv, j].data[e],
                                              spectral_phamat_spectral_path[i + 1][iv, j].data[e]);
        }
      }

      // VDISORT wants the first F11 value to be 1.0, so we normalize, as for DISORT
      if (std::isnormal(c[0, 0])) {
        c /= Numeric{c[0, 0]};
      } else {
        c       = 0.0;
        c[0, 0] = 1.0;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

////////////////////////////////////////////////////////////////////////////////
// VDISORT optical thicknesses
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsOpticalThicknessFromPath(VDisortSettings&                   vdisort_settings,
                                              const ArrayOfPropagationPathPoint& ray_path,
                                              const ArrayOfPropmatVector&        spectral_propmat_path,
                                              const Numeric&                     min_optical_depth) {
  ARTS_TIME_REPORT

  const Index N  = vdisort_settings.layer_count();
  const Index nv = vdisort_settings.frequency_count();

  ARTS_USER_ERROR_IF(ray_path.size() != static_cast<Size>(N + 1),
                     "*ray_path* must have one more point than *vdisort_settings* has layers: {} vs {}",
                     ray_path.size(),
                     N)

  if (N == 0) return;

  ARTS_USER_ERROR_IF(not all_same_shape({nv}, spectral_propmat_path),
                     "Propagation matrices and frequency grids must have the same shape.")

  // VDISORT's extinction is scalar, as DISORT's
  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_propmat_path,
                   [](const PropmatVector& pms) { return stdr::any_of(pms, Cmp::eq<true>(), &Propmat::is_polarized); }),
      "No implementation for polarized propagation matrices: VDISORT's optical thickness is scalar.");

  Vector r(N);
  for (Index i = 0; i < N; i++) { r[i] = ray_path[i].altitude() - ray_path[i + 1].altitude(); }

  ARTS_USER_ERROR_IF(stdr::any_of(r, Cmp::le<0.0>()),
                     R"(Atmospheric layer thickness must be positive.

Values:   {:B,}

Ray path points must be sorted by decreasing altitude.

ray_path: {:B,}
)",
                     r,
                     ray_path);

  for (Index iv = 0; iv < nv; iv++) {
    for (Index i = 0; i < N; i++) {
      vdisort_settings.optical_thicknesses[iv, i] =
          std::max(r[i] * std::midpoint(spectral_propmat_path[i + 1][iv].A(), spectral_propmat_path[i + 0][iv].A()),
                   min_optical_depth);
      if (i > 0) {
        vdisort_settings.optical_thicknesses[iv, i] += vdisort_settings.optical_thicknesses[iv, i - 1];

        ARTS_USER_ERROR_IF(
            (vdisort_settings.optical_thicknesses[iv, i] <= vdisort_settings.optical_thicknesses[iv, i - 1]),
            R"(
Not strictly increasing optical thicknesses between layers.

Check *spectral_propmat_path* contain zeroes or negative values for A().

Value:                {}
Frequency grid index: {}
)",
            vdisort_settings.optical_thicknesses[iv, i - 1],
            iv);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT BRDF / BDRF
////////////////////////////////////////////////////////////////////////////////

/* The surface is held as its reflection between the streams for each frequency, and its reflection of the sun's
   beam at the sun's zenith angle, so the sun must be set first.  A surface with the same parameters at every
   frequency is evaluated once, and again only where the sun's zenith angle changes. */

namespace {
//! The cosines of the upward streams
Vector upward_streams(const VDisortSettings& vdisort_settings) {
  const Index N = vdisort_settings.quadrature_dimension / 2;
  Vector      mu(2 * N), inv_mu(2 * N), W(N);
  disort_common::initialize_streams(mu, inv_mu, W);
  return Vector{mu[Range(0, N)]};
}

//! The cosine of the sun's zenith angle of frequency iv, or NaN for a sun below the horizon
Numeric sun_cosine(const VDisortSettings& vdisort_settings, const Index iv) {
  const Numeric mu0 = Conversion::cosd(vdisort_settings.solar_zenith_angle[iv]);
  return mu0 > 0.0 ? mu0 : std::numeric_limits<Numeric>::quiet_NaN();
}

//! Size the surface for NBRDF Fourier modes, with a specular part if specular
void resize_surface(VDisortSettings& vdisort_settings, const Index NBRDF, const bool specular) {
  const Index nv = vdisort_settings.frequency_count();
  const Index N  = vdisort_settings.quadrature_dimension / 2;
  ARTS_USER_ERROR_IF(vdisort_settings.solar_zenith_angle.size() != static_cast<Size>(nv),
                     "The surface reflects the sun, so the solar zenith angle must be set for all {} frequencies "
                     "first; it has {} values",
                     nv,
                     vdisort_settings.solar_zenith_angle.size());
  vdisort_settings.surface_brdf.resize(nv, 2, NBRDF, N, N);
  vdisort_settings.surface_beam_brdf.resize(nv, 2, NBRDF, N);
  vdisort_settings.surface_beam_mu.resize(nv);
  vdisort_settings.surface_specular.resize(nv, specular ? N : 0);
}

//! Set the surface of frequency iv
void set_surface(VDisortSettings& vdisort_settings, const Index iv, const vdisort::surface_reflection& surface) {
  vdisort_settings.surface_brdf[iv]      = surface.brdf();
  vdisort_settings.surface_beam_brdf[iv] = surface.beam_brdf();
  vdisort_settings.surface_beam_mu[iv]   = surface.beam_mu();
  vdisort_settings.surface_specular[iv]  = surface.specular();
}

//! Throw unless a spectral parameter has one value per frequency
void check_spectral(const VDisortSettings& vdisort_settings, const Vector& value, const std::string_view name) {
  ARTS_USER_ERROR_IF(static_cast<Index>(value.size()) != vdisort_settings.frequency_count(),
                     "Spectral BRDF parameter {} must have one value per frequency: expected {}, got {}",
                     name,
                     vdisort_settings.frequency_count(),
                     value.size());
}

//! A scalar BRDF as a fully depolarizing BPrDF
vdisort::brdf::RawFunction depolarizing(const disort::brdf::RawFunction& scalar) {
  return [scalar](const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) {
    rtepack::muelmat result{0.0};
    result[0, 0] = scalar(outgoing_mu, incoming_mu, relative_azimuth);
    return result;
  };
}

//! Set the surface of every frequency to the BPrDF raw
void set_brdf_surface(VDisortSettings&                  vdisort_settings,
                      const vdisort::brdf::RawFunction& raw,
                      const Index                       azimuth_quadrature_points) {
  const Vector mu = upward_streams(vdisort_settings);
  resize_surface(vdisort_settings, vdisort_settings.fourier_mode_dimension, false);
  vdisort::surface_reflection surface;
  for (Index iv = 0; iv < vdisort_settings.frequency_count(); iv++) {
    if (iv == 0 or vdisort_settings.solar_zenith_angle[iv] != vdisort_settings.solar_zenith_angle[iv - 1])
      surface = vdisort::brdf::fourier_modes(raw,
                                             vdisort_settings.fourier_mode_dimension,
                                             mu,
                                             mu,
                                             sun_cosine(vdisort_settings, iv),
                                             azimuth_quadrature_points);
    set_surface(vdisort_settings, iv, surface);
  }
}

//! Set the surface of each frequency iv to the BPrDF raw(iv)
void set_spectral_brdf_surface(VDisortSettings&                                        vdisort_settings,
                               const std::function<vdisort::brdf::RawFunction(Index)>& raw,
                               const Index                                             azimuth_quadrature_points) {
  const Vector mu = upward_streams(vdisort_settings);
  resize_surface(vdisort_settings, vdisort_settings.fourier_mode_dimension, false);
  for (Index iv = 0; iv < vdisort_settings.frequency_count(); iv++)
    set_surface(vdisort_settings,
                iv,
                vdisort::brdf::fourier_modes(raw(iv),
                                             vdisort_settings.fourier_mode_dimension,
                                             mu,
                                             mu,
                                             sun_cosine(vdisort_settings, iv),
                                             azimuth_quadrature_points));
}
}  // namespace

void vdisort_settingsNoSurfaceScattering(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  resize_surface(vdisort_settings, 0, false);
  vdisort_settings.surface_beam_mu = std::numeric_limits<Numeric>::quiet_NaN();
}

void vdisort_settingsSurfaceLambertian(VDisortSettings&                           vdisort_settings,
                                       const Generic<const Numeric, const Vector> value) {
  ARTS_TIME_REPORT

  const Index nv     = vdisort_settings.frequency_count();
  const auto  albedo = std::visit(
      [nv](const auto& selected) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(*selected)>, Numeric>)
          return Vector(nv, *selected);
        else
          return Vector{*selected};
      },
      value);
  check_spectral(vdisort_settings, albedo, "value");

  const Vector mu = upward_streams(vdisort_settings);
  const Index  NF = vdisort_settings.fourier_mode_dimension;
  resize_surface(vdisort_settings, NF, false);
  for (Index iv = 0; iv < nv; iv++)
    set_surface(vdisort_settings,
                iv,
                vdisort::brdf::lambertian_fourier_modes(albedo[iv], NF, mu, mu, sun_cosine(vdisort_settings, iv)));
}

void vdisort_settingsSurfaceFresnel(VDisortSettings& vdisort_settings,
                                    const Numeric&   refractive_index_real,
                                    const Numeric&   refractive_index_imag) {
  ARTS_TIME_REPORT

  // A mirror cannot reflect the beam (see vdisort::surface_reflection), so it has no beam reflection
  const Vector mu = upward_streams(vdisort_settings);
  const Index  NF = vdisort_settings.fourier_mode_dimension;
  resize_surface(vdisort_settings, NF, true);
  const auto surface = vdisort::brdf::fresnel_fourier_modes(
      Complex{refractive_index_real, refractive_index_imag}, NF, mu, mu, std::numeric_limits<Numeric>::quiet_NaN());
  for (Index iv = 0; iv < vdisort_settings.frequency_count(); iv++) set_surface(vdisort_settings, iv, surface);
}

void vdisort_settingsSurfaceFresnelSpectral(VDisortSettings& vdisort_settings,
                                            const Vector&    refractive_index_real,
                                            const Vector&    refractive_index_imag) {
  ARTS_TIME_REPORT

  check_spectral(vdisort_settings, refractive_index_real, "refractive_index_real");
  check_spectral(vdisort_settings, refractive_index_imag, "refractive_index_imag");

  // A mirror cannot reflect the beam (see vdisort::surface_reflection), so it has no beam reflection
  const Vector mu = upward_streams(vdisort_settings);
  const Index  NF = vdisort_settings.fourier_mode_dimension;
  resize_surface(vdisort_settings, NF, true);
  for (Index iv = 0; iv < vdisort_settings.frequency_count(); iv++)
    set_surface(vdisort_settings,
                iv,
                vdisort::brdf::fresnel_fourier_modes(Complex{refractive_index_real[iv], refractive_index_imag[iv]},
                                                     NF,
                                                     mu,
                                                     mu,
                                                     std::numeric_limits<Numeric>::quiet_NaN()));
}

void vdisort_settingsSurfaceHapke(VDisortSettings& vdisort_settings,
                                  const Numeric&   opposition_amplitude,
                                  const Numeric&   opposition_width,
                                  const Numeric&   single_scattering_albedo,
                                  const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  set_brdf_surface(vdisort_settings,
                   depolarizing(disort::brdf::Hapke{.opposition_amplitude     = opposition_amplitude,
                                                    .opposition_width         = opposition_width,
                                                    .single_scattering_albedo = single_scattering_albedo}),
                   azimuth_quadrature_points);
}

void vdisort_settingsSurfaceHapkeSpectral(VDisortSettings& vdisort_settings,
                                          const Vector&    opposition_amplitude,
                                          const Vector&    opposition_width,
                                          const Vector&    single_scattering_albedo,
                                          const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  check_spectral(vdisort_settings, opposition_amplitude, "opposition_amplitude");
  check_spectral(vdisort_settings, opposition_width, "opposition_width");
  check_spectral(vdisort_settings, single_scattering_albedo, "single_scattering_albedo");
  set_spectral_brdf_surface(
      vdisort_settings,
      [&](const Index iv) {
        return depolarizing(disort::brdf::Hapke{.opposition_amplitude     = opposition_amplitude[iv],
                                                .opposition_width         = opposition_width[iv],
                                                .single_scattering_albedo = single_scattering_albedo[iv]});
      },
      azimuth_quadrature_points);
}

void vdisort_settingsSurfaceCoxMunk(VDisortSettings& vdisort_settings,
                                    const Numeric&   wind_speed,
                                    const Numeric&   refractive_index,
                                    const Index&     shadowing,
                                    const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(shadowing != 0 and shadowing != 1, "shadowing must be 0 or 1, got {}", shadowing);
  set_brdf_surface(vdisort_settings,
                   vdisort::brdf::CoxMunk{.wind_speed       = wind_speed,
                                          .refractive_index = Complex{refractive_index, 0.0},
                                          .shadowing        = static_cast<bool>(shadowing)},
                   azimuth_quadrature_points);
}

void vdisort_settingsSurfaceCoxMunkSpectral(VDisortSettings& vdisort_settings,
                                            const Numeric&   wind_speed,
                                            const Vector&    refractive_index,
                                            const Index&     shadowing,
                                            const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(shadowing != 0 and shadowing != 1, "shadowing must be 0 or 1, got {}", shadowing);
  check_spectral(vdisort_settings, refractive_index, "refractive_index");
  set_spectral_brdf_surface(
      vdisort_settings,
      [&](const Index iv) -> vdisort::brdf::RawFunction {
        return vdisort::brdf::CoxMunk{.wind_speed       = wind_speed,
                                      .refractive_index = Complex{refractive_index[iv], 0.0},
                                      .shadowing        = static_cast<bool>(shadowing)};
      },
      azimuth_quadrature_points);
}

void vdisort_settingsSurfaceRPV(VDisortSettings& vdisort_settings,
                                const Numeric&   rho0,
                                const Numeric&   kappa,
                                const Numeric&   asymmetry,
                                const Numeric&   hotspot,
                                const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  set_brdf_surface(
      vdisort_settings,
      depolarizing(disort::brdf::RPV{.rho0 = rho0, .kappa = kappa, .asymmetry = asymmetry, .hotspot = hotspot}),
      azimuth_quadrature_points);
}

void vdisort_settingsSurfaceRPVSpectral(VDisortSettings& vdisort_settings,
                                        const Vector&    rho0,
                                        const Vector&    kappa,
                                        const Vector&    asymmetry,
                                        const Vector&    hotspot,
                                        const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  check_spectral(vdisort_settings, rho0, "rho0");
  check_spectral(vdisort_settings, kappa, "kappa");
  check_spectral(vdisort_settings, asymmetry, "asymmetry");
  check_spectral(vdisort_settings, hotspot, "hotspot");
  set_spectral_brdf_surface(
      vdisort_settings,
      [&](const Index iv) {
        return depolarizing(disort::brdf::RPV{
            .rho0 = rho0[iv], .kappa = kappa[iv], .asymmetry = asymmetry[iv], .hotspot = hotspot[iv]});
      },
      azimuth_quadrature_points);
}

void vdisort_settingsSurfaceRossLi(VDisortSettings& vdisort_settings,
                                   const Numeric&   isotropic,
                                   const Numeric&   volumetric,
                                   const Numeric&   geometric,
                                   const Numeric&   hotspot_angle,
                                   const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  set_brdf_surface(vdisort_settings,
                   depolarizing(disort::brdf::RossLi{.isotropic     = isotropic,
                                                     .volumetric    = volumetric,
                                                     .geometric     = geometric,
                                                     .hotspot_angle = Conversion::deg2rad(hotspot_angle)}),
                   azimuth_quadrature_points);
}

void vdisort_settingsSurfaceRossLiSpectral(VDisortSettings& vdisort_settings,
                                           const Vector&    isotropic,
                                           const Vector&    volumetric,
                                           const Vector&    geometric,
                                           const Numeric&   hotspot_angle,
                                           const Index&     azimuth_quadrature_points) {
  ARTS_TIME_REPORT

  check_spectral(vdisort_settings, isotropic, "isotropic");
  check_spectral(vdisort_settings, volumetric, "volumetric");
  check_spectral(vdisort_settings, geometric, "geometric");
  set_spectral_brdf_surface(
      vdisort_settings,
      [&](const Index iv) {
        return depolarizing(disort::brdf::RossLi{.isotropic     = isotropic[iv],
                                                 .volumetric    = volumetric[iv],
                                                 .geometric     = geometric[iv],
                                                 .hotspot_angle = Conversion::deg2rad(hotspot_angle)});
      },
      azimuth_quadrature_points);
}

////////////////////////////////////////////////////////////////////////////////
// VDISORT single scattering albedo
////////////////////////////////////////////////////////////////////////////////

void vdisort_settingsNoSingleScatteringAlbedo(VDisortSettings& vdisort_settings) {
  ARTS_TIME_REPORT

  vdisort_settings.single_scattering_albedo = 0.0;
}

void vdisort_settingsSingleScatteringAlbedoFromPath(VDisortSettings&            vdisort_settings,
                                                    const ArrayOfPropmatVector& spectral_propmat_path,
                                                    const ArrayOfPropmatVector& spectral_propmat_scat_path,
                                                    const ArrayOfStokvecVector& spectral_absvec_scat_path) try {
  ARTS_TIME_REPORT

  const Size  N = vdisort_settings.layer_count();
  const Index F = vdisort_settings.frequency_count();

  ARTS_USER_ERROR_IF(
      (N + 1) != spectral_propmat_path.size() or (N + 1) != spectral_propmat_scat_path.size() or
          (N + 1) != spectral_absvec_scat_path.size(),
      R"(The number of levels of the paths must be one more than the number of layers in vdisort_settings.

  vdisort_settings.layer_count():    {}
  spectral_propmat_path.size():      {}
  spectral_propmat_scat_path.size(): {}
  spectral_absvec_scat_path.size():  {}
)",
      N,
      spectral_propmat_path.size(),
      spectral_propmat_scat_path.size(),
      spectral_absvec_scat_path.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape({F}, spectral_absvec_scat_path, spectral_propmat_scat_path),
      "The shape of spectral_propmat_scat_path and spectral_absvec_scat_path must be {:B,}, at least one is not",
      std::array{F});

  ARTS_USER_ERROR_IF(not same_shape({F, N}, vdisort_settings.single_scattering_albedo),
                     R"(The shape of vdisort_settings.single_scattering_albedo must be {:B,}, but is {:B,})",
                     std::array{F, static_cast<Index>(N)},
                     vdisort_settings.single_scattering_albedo.shape());

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      const Numeric ext = std::midpoint(spectral_propmat_path[i][iv][0], spectral_propmat_path[i + 1][iv][0]);
      const Numeric abs_scat =
          std::midpoint(spectral_absvec_scat_path[i][iv][0], spectral_absvec_scat_path[i + 1][iv][0]);
      const Numeric ext_scat =
          std::midpoint(spectral_propmat_scat_path[i][iv][0], spectral_propmat_scat_path[i + 1][iv][0]);

      const Numeric x                                  = (ext_scat - abs_scat) / ext;
      vdisort_settings.single_scattering_albedo[iv, i] = std::isnormal(x) ? x : 0.0;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

////////////////////////////////////////////////////////////////////////////////
// Convenience
////////////////////////////////////////////////////////////////////////////////

void legendre_degreeFromVDisortSettings(Index& legendre_degree, const VDisortSettings& vdisort_settings) {
  legendre_degree = vdisort_settings.legendre_polynomial_dimension - 1;
  ARTS_USER_ERROR_IF(legendre_degree < 0, "The legendre_degree must be non-negative, is {}", legendre_degree)
}

////////////////////////////////////////////////////////////////////////////////
// Agenda setup
////////////////////////////////////////////////////////////////////////////////

void vdisort_settings_agendaSetup(Agenda&        vdisort_settings_agenda,
                                  const String&  layer_emission_setting,
                                  const String&  scattering_setting,
                                  const String&  space_setting,
                                  const String&  sun_setting,
                                  const String&  surf_setting,
                                  const Vector&  surf_lambertian_value,
                                  const Numeric& min_optical_depth) {
  ARTS_TIME_REPORT

  AgendaCreator agenda("vdisort_settings_agenda");
  agenda.add("jac_targetsOff");
  // Clearsky absorption
  agenda.add("atm_pathFromPath");
  agenda.add("freq_grid_pathFromPath");
  agenda.add("spectral_propmat_pathFromPath");
  agenda.add("vdisort_settingsInit");
  switch (to<settings_agenda_setup_scattering_type>(scattering_setting)) {
    using enum settings_agenda_setup_scattering_type;
    case ScatteringSpecies:
      agenda.add("legendre_degreeFromVDisortSettings");
      agenda.add("spectral_propmat_scat_pathFromSpectralAgenda");
      agenda.add("spectral_propmat_pathAddScattering");
      agenda.add("vdisort_settingsLegendreCoefficientsFromPath");
      agenda.add("vdisort_settingsSingleScatteringAlbedoFromPath");
      break;
    case None:
      agenda.add("vdisort_settingsNoLegendre");
      agenda.add("vdisort_settingsNoSingleScatteringAlbedo");
      break;
    case ScatteringSpeciesDeltaMPlus:
      ARTS_USER_ERROR(R"(VDISORT has no delta-M-plus scaling.

scattering_setting: "{}"

Use "ScatteringSpecies" or "None".
)",
                      scattering_setting);
  }
  // We have both scattering and clearsky absorption, so we can set the optical thickness
  agenda.add("vdisort_settingsOpticalThicknessFromPath", SetWsv("min_optical_depth", min_optical_depth));
  // Since we have the optical thickness, we can set the thermal emission
  switch (to<settings_agenda_setup_layer_emission_type>(layer_emission_setting)) {
    using enum settings_agenda_setup_layer_emission_type;
    case LinearInTau:       agenda.add("vdisort_settingsLayerThermalEmissionLinearInTau"); break;
    case LinearInTauNonLTE: agenda.add("vdisort_settingsLayerNonThermalEmissionLinearInTau"); break;
    case None:              agenda.add("vdisort_settingsNoLayerThermalEmission"); break;
  }
  switch (to<settings_agenda_setup_space_type>(space_setting)) {
    using enum settings_agenda_setup_space_type;
    case None:                               agenda.add("vdisort_settingsNoSpaceEmission"); break;
    case CosmicMicrowaveBackgroundRadiation: agenda.add("vdisort_settingsCosmicMicrowaveBackgroundRadiation"); break;
  }
  // The sun before the surface: the surface's reflection of the direct beam is computed for the sun's direction
  switch (to<settings_agenda_setup_sun_type>(sun_setting)) {
    using enum settings_agenda_setup_sun_type;
    case None: agenda.add("vdisort_settingsNoSun"); break;
    case Sun:
      agenda.add("ray_pointHighestFromPath");
      agenda.add("vdisort_settingsSetSun");
      break;
  }
  switch (to<settings_agenda_setup_surface_type>(surf_setting)) {
    using enum settings_agenda_setup_surface_type;
    case None:
      agenda.add("vdisort_settingsNoSurfaceEmission");
      agenda.add("vdisort_settingsNoSurfaceScattering");
      break;
    case Thermal:
      agenda.add("ray_pointLowestFromPath");
      agenda.add("vdisort_settingsSurfaceEmissionByTemperature");
      agenda.add("vdisort_settingsNoSurfaceScattering");
      break;
    case ThermalLambertian:
      agenda.add("ray_pointLowestFromPath");
      agenda.add("vdisort_settingsSurfaceEmissionByTemperature");
      agenda.add("vdisort_settingsSurfaceLambertian", SetWsv{"value", surf_lambertian_value});
      break;
    case Lambertian:
      agenda.add("vdisort_settingsNoSurfaceEmission");
      agenda.add("vdisort_settingsSurfaceLambertian", SetWsv{"value", surf_lambertian_value});
      break;
  }
  vdisort_settings_agenda = std::move(agenda).finalize(false);
}
