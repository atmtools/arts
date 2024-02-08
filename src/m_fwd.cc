#include <fwd.h>
#include <workspace.h>

#include <algorithm>
#include <memory>

#include "arts_omp.h"
#include "debug.h"
#include "fwd_spectral_radiance.h"
#include "matpack_data.h"
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
                                                        {latitude},
                                                        {longitude},
                                                        atmospheric_field,
                                                        surface_field,
                                                        lines,
                                                        cia,
                                                        xsec,
                                                        predef,
                                                        cia_extrapolation,
                                                        cia_robust);
}

void spectral_radiance_fieldPlanarGeometricFromOperator(
    StokvecGriddedField6& spectral_radiance_field,
    const SpectralRadianceOperator& spectral_radiance_operator,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& zenith_grid,
    const AscendingGrid& azimuth_grid) {
  spectral_radiance_field.grid<0>() = zenith_grid;
  spectral_radiance_field.grid<1>() = azimuth_grid;
  const auto& altitude_grid = spectral_radiance_field.grid<2>() =
      spectral_radiance_operator.altitude();
  const auto& latitude_grid = spectral_radiance_field.grid<3>() =
      spectral_radiance_operator.latitude();
  const auto& longitude_grid = spectral_radiance_field.grid<4>() =
      spectral_radiance_operator.longitude();
  spectral_radiance_field.grid<5>() = frequency_grid;

  spectral_radiance_field.grid_names = {"Zenith angle",
                                        "Azimuth angle",
                                        "Altitude",
                                        "Latitude",
                                        "Longitude",
                                        "Frequency"};

  spectral_radiance_field.data_name = "Spectral Radiance Field";

  const Index nza = spectral_radiance_field.grid<0>().size();
  const Index naa = spectral_radiance_field.grid<1>().size();
  const Index nalt = spectral_radiance_field.grid<2>().size();
  const Index nlat = spectral_radiance_field.grid<3>().size();
  const Index nlon = spectral_radiance_field.grid<4>().size();
  const Index nfreq = spectral_radiance_field.grid<5>().size();

  spectral_radiance_field.data.resize(nza, naa, nalt, nlat, nlon, nfreq);

  ARTS_USER_ERROR_IF(altitude_grid.size() < 2, "Must have some type of path")
  ARTS_USER_ERROR_IF(latitude_grid.size() != 1, "Latitude must be scalar")
  ARTS_USER_ERROR_IF(longitude_grid.size() != 1, "Longitude must be scalar")
  ARTS_USER_ERROR_IF(std::ranges::binary_search(zenith_grid, 90.0),
                     "Zenith angle must not be 90 degrees")
  const Numeric alt_low = altitude_grid.front();
  const Numeric alt_high = altitude_grid.back();
  const Numeric lat = latitude_grid[0];
  const Numeric lon = longitude_grid[0];

  const std::vector<fwd::path> upwards =
      spectral_radiance_operator.geometric_planar({alt_low, lat, lon},
                                                  {180, 0});
  const std::vector<fwd::path> downwards =
      spectral_radiance_operator.geometric_planar({alt_high, lat, lon}, {0, 0});

  matpack::matpack_data<std::vector<fwd::path>, 2> paths(nza, naa);

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        if (zenith_grid[i] > 90.0) {
          paths(i, j) = upwards;
        } else {
          paths(i, j) = downwards;
        }
        const Numeric scl = std::abs(1.0 / Conversion::cosd(zenith_grid[i]));
        for (auto& path : paths(i, j)) {
          path.point.azimuth() = azimuth_grid[j];
          path.point.zenith() = zenith_grid[i];
          path.distance *= scl;
        }
      }
    }

    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        for (Index n = 0; n < nfreq; ++n) {
          auto srad =
              spectral_radiance_operator(frequency_grid[n],
                                         paths(i, j),
                                         SpectralRadianceOperator::as_vector{});
          if (zenith_grid[i] < 90.0) {
            spectral_radiance_field.data(i, j, joker, 0, 0, n) = srad;
          } else {
            std::ranges::reverse(srad);
            spectral_radiance_field.data(i, j, joker, 0, 0, n) = srad;
          }
        }
      }
    }
  } else {
    String error{};

#pragma omp parallel for collapse(2)
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        try {
          if (zenith_grid[i] > 90.0) {
            paths(i, j) = upwards;
          } else {
            paths(i, j) = downwards;
          }
          const Numeric scl = std::abs(1.0 / Conversion::cosd(zenith_grid[i]));
          for (auto& path : paths(i, j)) {
            path.point.azimuth() = azimuth_grid[j];
            path.point.zenith() = zenith_grid[i];
            path.distance *= scl;
          }
        } catch (std::exception& e) {
#pragma omp critical
          error += e.what() + String{"\n"};
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)

#pragma omp parallel for collapse(3)
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        for (Index n = 0; n < nfreq; ++n) {
          try {
            auto srad = spectral_radiance_operator(
                frequency_grid[n],
                paths(i, j),
                SpectralRadianceOperator::as_vector{});
            if (zenith_grid[i] < 90.0) {
              spectral_radiance_field.data(i, j, joker, 0, 0, n) = srad;
            } else {
              std::ranges::reverse(srad);
              spectral_radiance_field.data(i, j, joker, 0, 0, n) = srad;
            }
          } catch (std::exception& e) {
#pragma omp critical
            error += e.what() + String{"\n"};
          }
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)
  }
}
