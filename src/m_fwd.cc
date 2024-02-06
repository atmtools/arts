#include <fwd.h>
#include <workspace.h>

#include <memory>

#include "arts_omp.h"
#include "debug.h"
#include "fwd_spectral_radiance.h"
#include "matpack_view.h"
#include "rtepack.h"
#include "sorted_grid.h"

void spectral_radiance_operator1D(
    const Workspace& ws,
    SpectralRadianceOperator& spectral_radiance_operator,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const AscendingGrid& altitude_grid,
    const Numeric& latitude,
    const Numeric& longitude,
    const Numeric& cia_extrapolation,
    const Index& cia_robust) {
  ARTS_USER_ERROR_IF(altitude_grid.size() < 2, "Must have some type of path")

  using lines_t = AbsorptionBands;
  using cia_t = ArrayOfCIARecord;
  using xsec_t = ArrayOfXsecRecord;
  using predef_t = PredefinedModelData;

  const String lines_str = "absorption_bands";
  const String cia_str = "absorption_cia_data";
  const String xsec_str = "absorption_xsec_fit_data";
  const String predef_str = "absorption_predefined_model_data";

  auto lines = ws.wsv_and_contains(lines_str)
                   ? ws.share(lines_str)->share<lines_t>()
                   : std::shared_ptr<lines_t>{};
  auto cia = ws.wsv_and_contains(cia_str) ? ws.share(cia_str)->share<cia_t>()
                                          : std::shared_ptr<cia_t>{};
  auto xsec = ws.wsv_and_contains(xsec_str)
                  ? ws.share(xsec_str)->share<xsec_t>()
                  : std::shared_ptr<xsec_t>{};
  auto predef = ws.wsv_and_contains(predef_str)
                    ? ws.share(predef_str)->share<predef_t>()
                    : std::shared_ptr<predef_t>{};

  spectral_radiance_operator = SpectralRadianceOperator(altitude_grid,
                                                        latitude,
                                                        longitude,
                                                        atmospheric_field,
                                                        surface_field,
                                                        std::move(lines),
                                                        std::move(cia),
                                                        std::move(xsec),
                                                        std::move(predef),
                                                        cia_extrapolation,
                                                        cia_robust);
}

void spectral_radiance_field1DOperator(
    StokvecTensor4& spectral_radiance_field,
    const SpectralRadianceOperator& spectral_radiance_operator,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& zenith_grid,
    const AscendingGrid& azimuth_grid) {
  const ExhaustiveConstVectorView altitude_grid =
      spectral_radiance_operator.altitude_grid();
  const ExhaustiveConstVectorView latitude =
      spectral_radiance_operator.latitude_grid();
  const ExhaustiveConstVectorView longitude =
      spectral_radiance_operator.longitude_grid();

  const Numeric alt_low = altitude_grid.front();
  const Numeric alt_high = altitude_grid.back();

  ARTS_USER_ERROR_IF(altitude_grid.size() < 2, "Must have some type of path")
  ARTS_USER_ERROR_IF(latitude.size() != 1 and longitude.size() != 1,
                     "Must be 1D")

  const auto nz = zenith_grid.size();
  const auto na = azimuth_grid.size();
  const auto nfreq = frequency_grid.size();

  spectral_radiance_field.resize(nz, na, nfreq, altitude_grid.size());

  if (arts_omp_in_parallel()) {
    for (Index i = 0; i < nz; i++) {
      auto path = spectral_radiance_operator.geometric_planar(
          zenith_grid[i] < 90.0 ? alt_high : alt_low, zenith_grid[i]);
      for (Index j = 0; j < na; j++) {
        for (Index k = 0; k < nfreq; k++) {
          spectral_radiance_operator(spectral_radiance_field(i, j, k, joker),
                                     frequency_grid[j],
                                     {zenith_grid[i], azimuth_grid[j]},
                                     path);
        }
      }
    }
  } else {
    String error{};

#pragma omp parallel for collapse(2)
    for (Index i = 0; i < nz; i++) {
      for (Index j = 0; j < na; j++) {
        try {
          const auto path = spectral_radiance_operator.geometric_planar(
              zenith_grid[i] <= 90.0 ? alt_high : alt_low, zenith_grid[i]);
          for (Index k = 0; k < nfreq; k++) {
            spectral_radiance_operator(spectral_radiance_field(i, j, k, joker),
                                       frequency_grid[j],
                                       {zenith_grid[i], azimuth_grid[j]},
                                       path);
          }
        } catch (std::exception& e) {
#pragma omp critical
          error += e.what() + String("\n");
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)
  }
}
