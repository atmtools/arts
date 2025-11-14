#include <array_algo.h>
#include <arts_omp.h>
#include <workspace.h>

namespace {
void ray_path_spectral_radStepByStepEmissionForwardOnly(
    ArrayOfStokvecVector& ray_path_spectral_rad,
    const ArrayOfMuelmatVector& spectral_tramat_path,
    const ArrayOfStokvecVector& spectral_rad_srcvec_path,
    const StokvecVector& spectral_rad_bkg) try {
  ARTS_TIME_REPORT

  ray_path_spectral_rad.resize(spectral_tramat_path.size());
  arr::elemwise_resize(spectral_rad_bkg.size(), ray_path_spectral_rad);

  ray_path_spectral_rad.back() = spectral_rad_bkg;

  two_level_linear_emission_step_by_step_full(
      ray_path_spectral_rad, spectral_tramat_path, spectral_rad_srcvec_path);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_spectral_radClearskyEmission(
    const Workspace& ws,
    ArrayOfStokvecVector& ray_path_spectral_rad,
    const AtmField& atm_field,
    const AscendingGrid& freq_grid,
    const Agenda& spectral_propmat_agenda,
    const ArrayOfPropagationPathPoint& ray_path,
    const Agenda& spectral_rad_space_agenda,
    const Agenda& spectral_rad_surface_agenda,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field) try {
  ARTS_TIME_REPORT

  PropagationPathPoint ray_point;
  ray_pointBackground(ray_point, ray_path);
  StokvecVector spectral_rad_bkg;
  StokvecMatrix spectral_rad_bkg_jac;
  spectral_rad_bkgAgendasAtEndOfPath(ws,
                                     spectral_rad_bkg,
                                     spectral_rad_bkg_jac,
                                     freq_grid,
                                     {},
                                     ray_point,
                                     surf_field,
                                     subsurf_field,
                                     spectral_rad_space_agenda,
                                     spectral_rad_surface_agenda);
  ArrayOfAtmPoint atm_path;
  atm_pathFromPath(atm_path, ray_path, atm_field);
  ArrayOfAscendingGrid freq_grid_path;
  ArrayOfVector3 freq_wind_shift_jac_path;
  freq_grid_pathFromPath(
      freq_grid_path, freq_wind_shift_jac_path, freq_grid, ray_path, atm_path);
  ArrayOfPropmatVector spectral_propmat_path;
  ArrayOfStokvecVector spectral_srcvec_nlte_path;
  ArrayOfPropmatMatrix spectral_propmat_jac_path;
  ArrayOfStokvecMatrix spectral_srcvec_nlte_jac_path;
  spectral_propmat_pathFromPath(ws,
                                spectral_propmat_path,
                                spectral_srcvec_nlte_path,
                                spectral_propmat_jac_path,
                                spectral_srcvec_nlte_jac_path,
                                spectral_propmat_agenda,
                                freq_grid_path,
                                freq_wind_shift_jac_path,
                                {},
                                ray_path,
                                atm_path);
  ArrayOfMuelmatVector spectral_tramat_path;
  ArrayOfMuelmatTensor3 spectral_tramat_jac_path;
  spectral_tramat_pathFromPath(spectral_tramat_path,
                               spectral_tramat_jac_path,
                               spectral_propmat_path,
                               spectral_propmat_jac_path,
                               ray_path,
                               atm_path,
                               surf_field,
                               {},
                               0);
  ArrayOfStokvecVector spectral_rad_srcvec_path;
  ArrayOfStokvecMatrix spectral_rad_srcvec_jac_path;
  spectral_rad_srcvec_pathFromPropmat(spectral_rad_srcvec_path,
                                      spectral_rad_srcvec_jac_path,
                                      spectral_propmat_path,
                                      spectral_srcvec_nlte_path,
                                      spectral_propmat_jac_path,
                                      spectral_srcvec_nlte_jac_path,
                                      freq_grid_path,
                                      atm_path,
                                      {});
  ray_path_spectral_radStepByStepEmissionForwardOnly(ray_path_spectral_rad,
                                                     spectral_tramat_path,
                                                     spectral_rad_srcvec_path,
                                                     spectral_rad_bkg);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace

void spectral_flux_profileFromPathField(
    const Workspace& ws,
    Matrix& spectral_flux_profile,
    const ArrayOfArrayOfPropagationPathPoint& ray_path_field,
    const AtmField& atm_field,
    const Agenda& spectral_propmat_agenda,
    const Agenda& spectral_rad_space_agenda,
    const Agenda& spectral_rad_surface_agenda,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const AscendingGrid& freq_grid,
    const AscendingGrid& alt_grid) try {
  ARTS_TIME_REPORT

  const Size N = ray_path_field.size();
  const Size M = alt_grid.size();
  const Size K = freq_grid.size();

  spectral_flux_profile.resize(M, K);
  spectral_flux_profile = 0.0;

  ArrayOfArrayOfStokvecVector ray_path_spectral_rad_field(N);

  String error{};
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size n = 0; n < N; n++) {
    try {
      ray_path_spectral_radClearskyEmission(ws,
                                            ray_path_spectral_rad_field[n],
                                            atm_field,
                                            freq_grid,
                                            spectral_propmat_agenda,
                                            ray_path_field[n],
                                            spectral_rad_space_agenda,
                                            spectral_rad_surface_agenda,
                                            surf_field,
                                            subsurf_field);
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);

  ARTS_USER_ERROR_IF(
      not arr::each_same_size(ray_path_spectral_rad_field, ray_path_field),
      "Not all ray paths have the same altitude count")

  for (auto& sradf : ray_path_spectral_rad_field) {
    ARTS_USER_ERROR_IF(
        stdr::any_of(sradf,
                     [K](const StokvecVector& v) { return v.size() != K; }),
        "Not all ray path spectral radiances in the field have the same frequency count")
  }

  struct Zenith {
    Size outer;
    Size inner;
    Numeric angle;
  };

  std::vector<Zenith> zenith_angles(2 * M);  // Up and down
  for (Size m = 0; m < M; m++) {
    zenith_angles.resize(0);
    const Numeric alt = alt_grid[m];
    VectorView t      = spectral_flux_profile[m];

    for (Size i = 0; i < ray_path_field.size(); i++) {
      for (Size j = 0; j < ray_path_field[i].size(); j++) {
        if (alt == ray_path_field[i][j].altitude()) {
          zenith_angles.emplace_back(i, j, ray_path_field[i][j].zenith());
        }
      }
    }

    ARTS_USER_ERROR_IF(
        zenith_angles.size() == 0, "No ray paths intersects altitude {} m", alt)

    stdr::sort(zenith_angles, {}, &Zenith::angle);

    // Integrate
    using Constant::pi;
    using Conversion::cosd;
    for (Size i = 0; i < zenith_angles.size() - 1; i++) {
      const auto& z0 = zenith_angles[i];
      const auto& z1 = zenith_angles[i + 1];

      const auto& y0 = ray_path_spectral_rad_field[z0.outer][z0.inner];
      const auto& y1 = ray_path_spectral_rad_field[z1.outer][z1.inner];

      const Numeric w = 0.5 * pi * (cosd(z0.angle) - cosd(z1.angle));
      for (Size k = 0; k < K; k++) t[k] += w * (y0[k][0] + y1[k][0]);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void flux_profileIntegrate(Vector& flux_profile,
                           const Matrix& spectral_flux_profile,
                           const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      static_cast<Index>(freq_grid.size()) != spectral_flux_profile.extent(1),
      "Frequency grid and spectral flux profile size mismatch")

  const Size K = spectral_flux_profile.extent(0);

  flux_profile.resize(K);
  flux_profile = 0.0;

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size k = 0; k < K; k++) {
    auto s  = spectral_flux_profile[k];
    auto& y = flux_profile[k];
    for (Size i = 0; i < freq_grid.size() - 1; i++) {
      const auto w  = 0.5 * (freq_grid[i + 1] - freq_grid[i]);
      y            += w * (s[i] + s[i + 1]);
    }
  }
}

void nlte_line_flux_profileIntegrate(
    QuantumIdentifierVectorMap& nlte_line_flux_profile,
    const Matrix& spectral_flux_profile,
    const AbsorptionBands& abs_bands,
    const ArrayOfAtmPoint& atm_path,
    const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  const Size K = spectral_flux_profile.extent(0);
  const Size M = spectral_flux_profile.extent(1);

  ARTS_USER_ERROR_IF(freq_grid.size() != M,
                     "Frequency grid and spectral flux profile size mismatch")

  ARTS_USER_ERROR_IF(
      atm_path.size() != K,
      "Atmospheric point and spectral flux profile size mismatch");

  nlte_line_flux_profile.clear();
  for (const auto& [key, band] : abs_bands) {
    ARTS_USER_ERROR_IF(band.size() != 1, "Only one line per band is supported");

    ARTS_USER_ERROR_IF(not is_voigt(band.lineshape),
                       "Only one line per band is supported");

    Matrix weighted_spectral_flux_profile(spectral_flux_profile.shape(), 0.0);

    for (Size k = 0; k < K; k++) {
      lbl::compute_voigt(weighted_spectral_flux_profile[k],
                         band.lines.front(),
                         freq_grid,
                         atm_path[k],
                         key.isot.mass);
    }

    weighted_spectral_flux_profile *= spectral_flux_profile;

    flux_profileIntegrate(
        nlte_line_flux_profile[key], weighted_spectral_flux_profile, freq_grid);
  }
}

void spectral_flux_profileFromSpectralRadianceField(
    Matrix& spectral_flux_profile,
    const GriddedSpectralField6& spectral_rad_field,
    const Stokvec& pol) {
  ARTS_TIME_REPORT

  using Constant::pi;
  using Conversion::cosd;
  using Conversion::deg2rad;

  const Size NALT = spectral_rad_field.grid<0>().size();
  const Size NLAT = spectral_rad_field.grid<1>().size();
  const Size NLON = spectral_rad_field.grid<2>().size();
  const Size NZA  = spectral_rad_field.grid<3>().size();
  const Size NAA  = spectral_rad_field.grid<4>().size();
  const Size NFRE = spectral_rad_field.grid<5>().size();

  const auto& za = spectral_rad_field.grid<3>();
  const auto& aa = spectral_rad_field.grid<4>();

  ARTS_USER_ERROR_IF(not spectral_rad_field.ok(),
                     R"(Spectral radiance field is not OK (shape is wrong):

spectral_rad_field.data.shape(): {:B,}
Should be:                            [NALT, NLAT, NLON, NZA, NAA, NFRE]

NALT: {}
NLAT: {}
NLON: {}
NZA:  {}
NAA:  {}
NFRE: {}
)",
                     spectral_rad_field.data.shape(),
                     NALT,
                     NLAT,
                     NLON,
                     NZA,
                     NAA,
                     NFRE);
  ARTS_USER_ERROR_IF(NLAT != 1, "Only for one latitude point");
  ARTS_USER_ERROR_IF(NLON != 1, "Only for one longitude point");
  ARTS_USER_ERROR_IF(NZA < 2, "Must have more than one zenith angle.");
  ARTS_USER_ERROR_IF(NAA != 1, "Only for one azimuth angle.")

  spectral_flux_profile.resize(NALT, NFRE);
  spectral_flux_profile = 0.0;

  if (NAA == 1) {
    const auto&& s = spectral_rad_field.data.view_as(NALT, NZA, NFRE);

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
    for (Size i = 0; i < NALT; i++) {
      for (Size j = 0; j < NFRE; j++) {
        for (Size iza0 = 0; iza0 < NZA - 1; iza0++) {
          const Size iza1   = iza0 + 1;
          const Numeric wza = 0.5 * pi * (cosd(za[iza0]) - cosd(za[iza1]));

          spectral_flux_profile[i, j] +=
              wza * (dot(s[i, iza0, j], pol) + dot(s[i, iza1, j], pol));
        }
      }
    }
  } else {
    const auto&& s = spectral_rad_field.data.view_as(NALT, NZA, NAA, NFRE);

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
    for (Size i = 0; i < NALT; i++) {
      for (Size j = 0; j < NFRE; j++) {
        for (Size iza0 = 0; iza0 < NZA - 1; iza0++) {
          const Size iza1   = iza0 + 1;
          const Numeric wza = 0.25 * (cosd(za[iza0]) - cosd(za[iza1]));

          for (Size iaa0 = 0; iaa0 < NAA - 1; iaa0++) {
            const Size iaa1 = iaa0 + 1;
            const Numeric w = wza * deg2rad(aa[iaa0] - aa[iaa1]);

            spectral_flux_profile[i, j] +=
                w *
                (dot(s[i, iza0, iaa0, j], pol) + dot(s[i, iza1, iaa0, j], pol) +
                 dot(s[i, iza0, iaa1, j], pol) + dot(s[i, iza1, iaa1, j], pol));
          }
        }
      }
    }
  }
}
