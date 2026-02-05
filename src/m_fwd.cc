#include <arts_omp.h>
#include <debug.h>
#include <fwd.h>
#include <fwd_path.h>
#include <fwd_spectral_radiance.h>
#include <obsel.h>
#include <path_point.h>
#include <rtepack.h>
#include <workspace.h>

#include <algorithm>
#include <exception>
#include <memory>
#include <ranges>

#include "workspace_class.h"

void spectral_rad_operatorClearsky1D(
    const Workspace& ws,
    SpectralRadianceOperator& spectral_rad_operator,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const AscendingGrid& alt_grid,
    const Numeric& latitude,
    const Numeric& longitude,
    const Numeric& cia_extrapolation,
    const Index& cia_robust) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ARTS_USER_ERROR_IF(alt_grid.size() < 2, "Must have some type of path")

  using lines_t  = AbsorptionBands;
  using cia_t    = CIARecords;
  using xsec_t   = XsecRecords;
  using predef_t = PredefinedModelData;

  const String lines_str  = "abs_bands";
  const String cia_str    = "abs_cia_data";
  const String xsec_str   = "abs_xfit_data";
  const String predef_str = "abs_predef_data";

  auto lines = ws.wsv_and_contains(lines_str)
                   ? ws.share(lines_str).share<lines_t>()
                   : std::shared_ptr<lines_t>{};
  auto cia   = ws.wsv_and_contains(cia_str) ? ws.share(cia_str).share<cia_t>()
                                            : std::shared_ptr<cia_t>{};
  auto xsec = ws.wsv_and_contains(xsec_str) ? ws.share(xsec_str).share<xsec_t>()
                                            : std::shared_ptr<xsec_t>{};
  auto predef = ws.wsv_and_contains(predef_str)
                    ? ws.share(predef_str).share<predef_t>()
                    : std::shared_ptr<predef_t>{};

  spectral_rad_operator = SpectralRadianceOperator(alt_grid,
                                                   Vector{latitude},
                                                   Vector{longitude},
                                                   atm_field,
                                                   surf_field,
                                                   lines,
                                                   cia,
                                                   xsec,
                                                   predef,
                                                   cia_extrapolation,
                                                   cia_robust);
}

void spectral_rad_fieldFromOperatorPlanarGeometric(
    GriddedSpectralField6& spectral_rad_field,
    const SpectralRadianceOperator& spectral_rad_operator,
    const AscendingGrid& freq_grid,
    const ZenGrid& zen_grid,
    const AziGrid& azi_grid) {
  ARTS_TIME_REPORT

  const AscendingGrid& alt_grid = spectral_rad_operator.altitude();
  const LatGrid& lat_grid       = spectral_rad_operator.latitude();
  const LonGrid& lon_grid       = spectral_rad_operator.longitude();

  const Index nza   = zen_grid.size();
  const Index naa   = azi_grid.size();
  const Index nalt  = alt_grid.size();
  const Index nlat  = lat_grid.size();
  const Index nlon  = lon_grid.size();
  const Index nfreq = freq_grid.size();

  spectral_rad_field = GriddedSpectralField6{
      .data_name = "Spectral Radiance Field",
      .data      = StokvecTensor6(
          nalt, nlat, nlon, nza, naa, nfreq, Stokvec{0.0, 0.0, 0.0, 0.0}),
      .grid_names = {"Altitude",
                     "Latitude",
                     "Longitude",
                     "Zenith angle",
                     "Azimuth angle",
                     "Frequency"},
      .grids = {alt_grid, lat_grid, lon_grid, zen_grid, azi_grid, freq_grid}};

  ARTS_USER_ERROR_IF(alt_grid.size() < 2, "Must have some type of path")
  ARTS_USER_ERROR_IF(lat_grid.size() != 1, "Latitude must be scalar")
  ARTS_USER_ERROR_IF(lon_grid.size() != 1, "Longitude must be scalar")
  ARTS_USER_ERROR_IF(stdr::binary_search(zen_grid, 90.0),
                     "Zenith angle must not be 90 degrees")
  const Numeric alt_low  = alt_grid.front();
  const Numeric alt_high = alt_grid.back();
  const Numeric lat      = lat_grid[0];
  const Numeric lon      = lon_grid[0];

  const std::vector<fwd::path> upwards =
      spectral_rad_operator.geometric_planar({alt_low, lat, lon}, {180, 0});
  const std::vector<fwd::path> downwards =
      spectral_rad_operator.geometric_planar({alt_high, lat, lon}, {0, 0});

  const auto pathstep = [&upwards, &downwards](const Numeric za,
                                               const Numeric az) {
    auto path         = (za > 90.0) ? upwards : downwards;
    const Numeric scl = std::abs(1.0 / Conversion::cosd(za));
    for (auto& pp : path) {
      pp.point.azimuth()  = az;
      pp.point.zenith()   = za;
      pp.distance        *= scl;
    }
    return path;
  };

  const auto freqstep = [&spectral_rad_operator](
                            const Numeric freq,
                            const Numeric za,
                            const std::vector<fwd::path>& path) {
    StokvecVector srad;
    if (za < 90.0) {
      srad = spectral_rad_operator(
          freq, path, SpectralRadianceOperator::as_vector{});
    } else {
      srad = spectral_rad_operator(
          freq, path, SpectralRadianceOperator::as_vector{});
      stdr::reverse(srad);
    }
    return srad;
  };

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        const auto path = pathstep(zen_grid[i], azi_grid[j]);
        for (Index n = 0; n < nfreq; ++n) {
          spectral_rad_field.data[joker, 0, 0, i, j, n] =
              freqstep(freq_grid[n], zen_grid[i], path);
        }
      }
    }
  } else {
    String error{};

#pragma omp parallel for collapse(2)
    for (Index i = 0; i < nza; ++i) {
      for (Index j = 0; j < naa; ++j) {
        try {
          const auto path = pathstep(zen_grid[i], azi_grid[j]);
          for (Index n = 0; n < nfreq; ++n) {
            spectral_rad_field.data[joker, 0, 0, i, j, n] =
                freqstep(freq_grid[n], zen_grid[i], path);
          }
        } catch (std::exception& e) {
#pragma omp critical
          error += e.what() + String{"\n"};
        }
      }
    }

    if (not error.empty()) throw std::runtime_error(error);
  }
}

void spectral_rad_fieldFromOperatorPath(
    const Workspace& ws,
    GriddedSpectralField6& spectral_rad_field,
    const SpectralRadianceOperator& spectral_rad_operator,
    const Agenda& ray_path_observer_agenda,
    const AscendingGrid& freq_grid,
    const ZenGrid& zen_grid,
    const AziGrid& azi_grid) {
  ARTS_TIME_REPORT

  const AscendingGrid& alt_grid = spectral_rad_operator.altitude();
  const LatGrid& lat_grid       = spectral_rad_operator.latitude();
  const LonGrid& lon_grid       = spectral_rad_operator.longitude();

  const Index nza   = zen_grid.size();
  const Index naa   = azi_grid.size();
  const Index nalt  = alt_grid.size();
  const Index nlat  = lat_grid.size();
  const Index nlon  = lon_grid.size();
  const Index nfreq = freq_grid.size();

  spectral_rad_field = GriddedSpectralField6{
      .data_name = "spectral_rad_fieldFromOperatorPath",
      .data      = StokvecTensor6(
          nalt, nlat, nlon, nza, naa, nfreq, Stokvec{0.0, 0.0, 0.0, 0.0}),
      .grid_names = {"Altitude",
                     "Latitude",
                     "Longitude",
                     "Zenith angle",
                     "Azimuth angle",
                     "Frequency"},
      .grids = {alt_grid, lat_grid, lon_grid, zen_grid, azi_grid, freq_grid}};

  if (arts_omp_in_parallel()) {
    for (Index ialt = 0; ialt < nalt; ++ialt) {
      for (Index ilat = 0; ilat < nlat; ++ilat) {
        for (Index ilon = 0; ilon < nlon; ++ilon) {
          for (Index iza = 0; iza < nza; ++iza) {
            for (Index iaa = 0; iaa < naa; ++iaa) {
              ArrayOfPropagationPathPoint ray_path;
              ray_path_observer_agendaExecute(
                  ws,
                  ray_path,
                  {alt_grid[ialt], lat_grid[ilat], lon_grid[ilon]},
                  {zen_grid[iza], azi_grid[iaa]},
                  ray_path_observer_agenda);
              stdr::transform(
                  freq_grid,
                  spectral_rad_field[ialt, ilat, ilon, iza, iaa, joker].begin(),
                  [path = spectral_rad_operator.from_path(ray_path),
                   &spectral_rad_operator](Numeric f) {
                    return spectral_rad_operator(f, path);
                  });
            }
          }
        }
      }
    }
  } else {
    String errors{};

#pragma omp parallel for collapse(5)
    for (Index ialt = 0; ialt < nalt; ++ialt) {
      for (Index ilat = 0; ilat < nlat; ++ilat) {
        for (Index ilon = 0; ilon < nlon; ++ilon) {
          for (Index iza = 0; iza < nza; ++iza) {
            for (Index iaa = 0; iaa < naa; ++iaa) {
              try {
                ArrayOfPropagationPathPoint ray_path;
                ray_path_observer_agendaExecute(
                    ws,
                    ray_path,
                    {alt_grid[ialt], lat_grid[ilat], lon_grid[ilon]},
                    {zen_grid[iza], azi_grid[iaa]},
                    ray_path_observer_agenda);
                stdr::transform(
                    freq_grid,
                    spectral_rad_field[ialt, ilat, ilon, iza, iaa, joker]
                        .begin(),
                    [path = spectral_rad_operator.from_path(ray_path),
                     &spectral_rad_operator](Numeric f) {
                      return spectral_rad_operator(f, path);
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

    ARTS_USER_ERROR_IF(not errors.empty(), "{}", errors)
  }
}

void measurement_vecFromOperatorPath(
    const Workspace& ws,
    Vector& measurement_vec,
    const ArrayOfSensorObsel& measurement_vec_sensor,
    const SpectralRadianceOperator& spectral_rad_operator,
    const Agenda& ray_path_observer_agenda) try {
  ARTS_TIME_REPORT

  measurement_vec.resize(measurement_vec_sensor.size());
  measurement_vec = 0.0;
  if (measurement_vec_sensor.empty()) return;

  //! Check the observational elements that their dimensions are correct
  for (auto& obsel : measurement_vec_sensor) obsel.check();

  const SensorSimulations simulations =
      collect_simulations(measurement_vec_sensor);

  for (auto& [freq_grid, poslos_grid, iposlos] : simulations) {
    const auto& poslos = poslos_grid[iposlos];

    ArrayOfPropagationPathPoint ray_path;

    ray_path_observer_agendaExecute(
        ws, ray_path, poslos.pos, poslos.los, ray_path_observer_agenda);

    const StokvecVector spectral_rad{
        std::from_range,
        freq_grid |
            stdv::transform([path = spectral_rad_operator.from_path(ray_path),
                             &spectral_rad_operator](Numeric f) {
              return spectral_rad_operator(f, path);
            })};

    for (Size iv = 0; iv < measurement_vec_sensor.size(); ++iv) {
      const SensorObsel& obsel = measurement_vec_sensor[iv];
      if (obsel.same_freqs(freq_grid)) {
        measurement_vec[iv] += obsel.sumup(spectral_rad, iposlos);
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
