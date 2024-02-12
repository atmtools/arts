#include <fwd.h>
#include <workspace.h>

#include <algorithm>
#include <exception>
#include <memory>

#include "arts_omp.h"
#include "debug.h"
#include "fwd_path.h"
#include "fwd_spectral_radiance.h"
#include "path_point.h"
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

void spectral_radiance_fieldFromOperatorPlanarGeometric(
    StokvecGriddedField6& spectral_radiance_field,
    const SpectralRadianceOperator& spectral_radiance_operator,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& zenith_grid,
    const AscendingGrid& azimuth_grid) {
  const AscendingGrid& altitude_grid = spectral_radiance_operator.altitude();
  const AscendingGrid& latitude_grid = spectral_radiance_operator.latitude();
  const AscendingGrid& longitude_grid = spectral_radiance_operator.longitude();

  const Index nza = zenith_grid.size();
  const Index naa = azimuth_grid.size();
  const Index nalt = altitude_grid.size();
  const Index nlat = latitude_grid.size();
  const Index nlon = longitude_grid.size();
  const Index nfreq = frequency_grid.size();

  spectral_radiance_field = StokvecGriddedField6{
      .data_name = "Spectral Radiance Field",
      .data = StokvecTensor6(
          nza, naa, nalt, nlat, nlon, nfreq, Stokvec{0.0, 0.0, 0.0, 0.0}),
      .grid_names = {"Zenith angle",
                     "Azimuth angle",
                     "Altitude",
                     "Latitude",
                     "Longitude",
                     "Frequency"},
      .grids = {zenith_grid,
                azimuth_grid,
                altitude_grid,
                latitude_grid,
                longitude_grid,
                frequency_grid}};

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

  const auto pathstep = [&upwards, &downwards](const Numeric za,
                                               const Numeric az) {
    auto path = (za > 90.0) ? upwards : downwards;
    const Numeric scl = std::abs(1.0 / Conversion::cosd(za));
    for (auto& pp : path) {
      pp.point.azimuth() = az;
      pp.point.zenith() = za;
      pp.distance *= scl;
    }
    return path;
  };

  const auto freqstep = [&spectral_radiance_operator](
                            const Numeric freq,
                            const Numeric za,
                            const std::vector<fwd::path>& path) {
    if (za < 90.0) {
      return spectral_radiance_operator(
          freq, path, SpectralRadianceOperator::as_vector{});
    }
    auto srad = spectral_radiance_operator(
        freq, path, SpectralRadianceOperator::as_vector{});
    std::ranges::reverse(srad);
    return srad;
  };

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        const auto path = pathstep(zenith_grid[i], azimuth_grid[j]);
        for (Index n = 0; n < nfreq; ++n) {
          spectral_radiance_field.data(i, j, joker, 0, 0, n) =
              freqstep(frequency_grid[n], zenith_grid[i], path);
        }
      }
    }
  } else {
    String error{};

#pragma omp parallel for collapse(2)
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        try {
          const auto path = pathstep(zenith_grid[i], azimuth_grid[j]);
          for (Index n = 0; n < nfreq; ++n) {
            spectral_radiance_field.data(i, j, joker, 0, 0, n) =
                freqstep(frequency_grid[n], zenith_grid[i], path);
          }
        } catch (std::exception& e) {
#pragma omp critical
          error += e.what() + String{"\n"};
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)
  }
}

void spectral_radiance_fieldFromOperatorPath(
    const Workspace& ws,
    StokvecGriddedField6& spectral_radiance_field,
    const SpectralRadianceOperator& spectral_radiance_operator,
    const Agenda& propagation_path_observer_agenda,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& zenith_grid,
    const AscendingGrid& azimuth_grid) {
  const AscendingGrid& altitude_grid = spectral_radiance_operator.altitude();
  const AscendingGrid& latitude_grid = spectral_radiance_operator.latitude();
  const AscendingGrid& longitude_grid = spectral_radiance_operator.longitude();

  const Index nza = zenith_grid.size();
  const Index naa = azimuth_grid.size();
  const Index nalt = altitude_grid.size();
  const Index nlat = latitude_grid.size();
  const Index nlon = longitude_grid.size();
  const Index nfreq = frequency_grid.size();

  spectral_radiance_field = StokvecGriddedField6{
      .data_name = "Spectral Radiance Field",
      .data = StokvecTensor6(
          nza, naa, nalt, nlat, nlon, nfreq, Stokvec{0.0, 0.0, 0.0, 0.0}),
      .grid_names = {"Zenith angle",
                     "Azimuth angle",
                     "Altitude",
                     "Latitude",
                     "Longitude",
                     "Frequency"},
      .grids = {zenith_grid,
                azimuth_grid,
                altitude_grid,
                latitude_grid,
                longitude_grid,
                frequency_grid}};

  if (arts_omp_in_parallel()) {
    for (Index iza = 0; iza < nza; ++iza) {
      for (Index iaa = 0; iaa < naa; ++iaa) {
        for (Index ialt = 0; ialt < nalt; ++ialt) {
          for (Index ilat = 0; ilat < nlat; ++ilat) {
            for (Index ilon = 0; ilon < nlon; ++ilon) {
              ArrayOfPropagationPathPoint propagation_path;
              propagation_path_observer_agendaExecute(
                  ws,
                  propagation_path,
                  {altitude_grid[ialt],
                   latitude_grid[ilat],
                   longitude_grid[ilon]},
                  {zenith_grid[iza], azimuth_grid[iaa]},
                  propagation_path_observer_agenda);
              std::transform(
                  frequency_grid.begin(),
                  frequency_grid.end(),
                  spectral_radiance_field(iza, iaa, ialt, ilat, ilon, joker)
                      .begin(),
                  [path = fwd::path_from_propagation_path(propagation_path,
                                                          altitude_grid,
                                                          latitude_grid,
                                                          longitude_grid),
                   &spectral_radiance_operator](Numeric f) {
                    return spectral_radiance_operator(f, path);
                  });
            }
          }
        }
      }
    }
  } else {
    String errors{};

#pragma omp parallel for collapse(5)
    for (Index iza = 0; iza < nza; ++iza) {
      for (Index iaa = 0; iaa < naa; ++iaa) {
        for (Index ialt = 0; ialt < nalt; ++ialt) {
          for (Index ilat = 0; ilat < nlat; ++ilat) {
            for (Index ilon = 0; ilon < nlon; ++ilon) {
              try {
                ArrayOfPropagationPathPoint propagation_path;
                propagation_path_observer_agendaExecute(
                    ws,
                    propagation_path,
                    {altitude_grid[ialt],
                     latitude_grid[ilat],
                     longitude_grid[ilon]},
                    {zenith_grid[iza], azimuth_grid[iaa]},
                    propagation_path_observer_agenda);
                std::transform(
                    frequency_grid.begin(),
                    frequency_grid.end(),
                    spectral_radiance_field(iza, iaa, ialt, ilat, ilon, joker)
                        .begin(),
                    [path = fwd::path_from_propagation_path(propagation_path,
                                                            altitude_grid,
                                                            latitude_grid,
                                                            longitude_grid),
                     &spectral_radiance_operator](Numeric f) {
                      return spectral_radiance_operator(f, path);
                    });
              } catch (std::exception& e) {
#pragma omp critical
                errors += e.what() + String("\n");
              }
            }
          }
        }
      }
    }

    ARTS_USER_ERROR_IF(not errors.empty(), errors)
  }
}
