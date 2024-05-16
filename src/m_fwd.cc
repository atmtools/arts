#include <fwd.h>
#include <workspace.h>

#include <algorithm>
#include <exception>
#include <memory>

#include "arts_omp.h"
#include "debug.h"
#include "enums.h"
#include "fwd_path.h"
#include "fwd_spectral_radiance.h"
#include "matpack_constexpr.h"
#include "obsel.h"
#include "path_point.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "surf.h"
#include "workspace_class.h"

void spectral_radiance_operatorClearsky1D(
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

  using lines_t = ArrayOfAbsorptionBand;
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
                  [path =
                       spectral_radiance_operator.from_path(propagation_path),
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
                    [path =
                         spectral_radiance_operator.from_path(propagation_path),
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

void measurement_vectorFromOperatorPath(
    const Workspace& ws,
    Vector& measurement_vector,
    const ArrayOfSensorObsel& measurement_vector_sensor,
    const SpectralRadianceOperator& spectral_radiance_operator,
    const Agenda& propagation_path_observer_agenda,
    const Index& exhaustive_) try {
  //! Flag whether or not all frequency and pos-los grids are to be assumed the same
  const bool exhaustive = static_cast<bool>(exhaustive_);

  ARTS_USER_ERROR_IF(not all_ok(measurement_vector_sensor),
                     "Measurement vector sensor infromation is not OK")
  ARTS_USER_ERROR_IF(
      exhaustive and not sensor::is_exhaustive_like(measurement_vector_sensor),
      "Measurement vector sensor infromation is not exhaustive-like despite exhaustive flag")

  measurement_vector.resize(measurement_vector_sensor.size());
  measurement_vector = 0.0;
  if (measurement_vector_sensor.empty()) return;

  ArrayOfPropagationPathPoint propagation_path;
  AscendingGrid frequency_grid;
  std::vector<fwd::path> path;
  StokvecVector spectral_radiance;

  if (exhaustive) {
    frequency_grid = measurement_vector_sensor.front().f_grid;
    spectral_radiance.resize(frequency_grid.size());
  }

  for (const auto poslos :
       (exhaustive ? measurement_vector_sensor.front().poslos_grid
                   : sensor::collect_poslos(measurement_vector_sensor))) {
    propagation_path_observer_agendaExecute(ws,
                                            propagation_path,
                                            poslos.pos,
                                            poslos.los,
                                            propagation_path_observer_agenda);
    spectral_radiance_operator.from_path(path, propagation_path);

    if (not exhaustive) {
      sensor::collect_f_grid(frequency_grid, measurement_vector_sensor, poslos);
      spectral_radiance.resize(frequency_grid.size());
    }

    if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
      std::transform(frequency_grid.begin(),
                     frequency_grid.end(),
                     spectral_radiance.begin(),
                     [&path, &spectral_radiance_operator](Numeric f) {
                       return spectral_radiance_operator(f, path);
                     });
    } else {
      String error{};

#pragma omp parallel for
      for (Index i = 0; i < frequency_grid.size(); ++i) {
        try {
          spectral_radiance[i] =
              spectral_radiance_operator(frequency_grid[i], path);
        } catch (std::exception& e) {
#pragma omp critical
          error += e.what() + String{"\n"};
        }
      }

      ARTS_USER_ERROR_IF(not error.empty(), error)
    }

    if (exhaustive) {
      sensor::exhaustive_sumup(measurement_vector,
                               spectral_radiance,
                               measurement_vector_sensor,
                               poslos);
    } else {
      sensor::sumup(measurement_vector,
                    spectral_radiance,
                    frequency_grid,
                    measurement_vector_sensor,
                    poslos);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceClearsky(
    const Workspace& ws,
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AtmField& atmospheric_field,
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets,
    const Agenda& propagation_matrix_agenda,
    const ArrayOfPropagationPathPoint& propagation_path,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda,
    const String& spectral_radiance_unit,
    const SurfaceField& surface_field,
    const Numeric& rte_alonglos_v,
    const Index& hse_derivative) try {
  PropagationPathPoint propagation_path_point;
  propagation_path_pointBackground(propagation_path_point, propagation_path);
  StokvecVector spectral_radiance_background;
  StokvecMatrix spectral_radiance_background_jacobian;
  spectral_radiance_backgroundAgendasAtEndOfPath(
      ws,
      spectral_radiance_background,
      spectral_radiance_background_jacobian,
      frequency_grid,
      jacobian_targets,
      propagation_path_point,
      surface_field,
      spectral_radiance_space_agenda,
      spectral_radiance_surface_agenda);
  ArrayOfAtmPoint propagation_path_atmospheric_point;
  propagation_path_atmospheric_pointFromPath(
      propagation_path_atmospheric_point, propagation_path, atmospheric_field);
  ArrayOfAscendingGrid propagation_path_frequency_grid;
  propagation_path_frequency_gridFromPath(propagation_path_frequency_grid,
                                          frequency_grid,
                                          propagation_path,
                                          propagation_path_atmospheric_point,
                                          rte_alonglos_v);
  ArrayOfPropmatVector propagation_path_propagation_matrix;
  ArrayOfStokvecVector propagation_path_propagation_matrix_source_vector_nonlte;
  ArrayOfPropmatMatrix propagation_path_propagation_matrix_jacobian;
  ArrayOfStokvecMatrix
      propagation_path_propagation_matrix_source_vector_nonlte_jacobian;
  propagation_path_propagation_matrixFromPath(
      ws,
      propagation_path_propagation_matrix,
      propagation_path_propagation_matrix_source_vector_nonlte,
      propagation_path_propagation_matrix_jacobian,
      propagation_path_propagation_matrix_source_vector_nonlte_jacobian,
      propagation_matrix_agenda,
      jacobian_targets,
      propagation_path_frequency_grid,
      propagation_path,
      propagation_path_atmospheric_point);
  ArrayOfMuelmatVector propagation_path_transmission_matrix;
  ArrayOfArrayOfMuelmatMatrix propagation_path_transmission_matrix_jacobian;
  propagation_path_transmission_matrixFromPath(
      propagation_path_transmission_matrix,
      propagation_path_transmission_matrix_jacobian,
      propagation_path_propagation_matrix,
      propagation_path_propagation_matrix_jacobian,
      propagation_path,
      propagation_path_atmospheric_point,
      surface_field,
      jacobian_targets,
      hse_derivative);
  ArrayOfMuelmatVector propagation_path_transmission_matrix_cumulative;
  propagation_path_transmission_matrix_cumulativeForward(
      propagation_path_transmission_matrix_cumulative,
      propagation_path_transmission_matrix);
  ArrayOfStokvecVector propagation_path_spectral_radiance_source;
  ArrayOfStokvecMatrix propagation_path_spectral_radiance_source_jacobian;
  propagation_path_spectral_radiance_sourceFromPropmat(
      propagation_path_spectral_radiance_source,
      propagation_path_spectral_radiance_source_jacobian,
      propagation_path_propagation_matrix,
      propagation_path_propagation_matrix_source_vector_nonlte,
      propagation_path_propagation_matrix_jacobian,
      propagation_path_propagation_matrix_source_vector_nonlte_jacobian,
      propagation_path_frequency_grid,
      propagation_path_atmospheric_point,
      jacobian_targets);
  // propagation_path_spectral_radiance_sourceAddSuns(
  //     propagation_path_spectral_radiance_source,
  //     propagation_path_spectral_radiance_source_jacobian,
  //     propagation_path_propagation_matrix,
  //     propagation_path,
  //     jacobian_targets);
  ArrayOfStokvecVector propagation_path_spectral_radiance;
  ArrayOfStokvecMatrix propagation_path_spectral_radiance_jacobian;
  propagation_path_spectral_radianceCalcEmission(
      propagation_path_spectral_radiance,
      propagation_path_spectral_radiance_jacobian,
      spectral_radiance_background,
      propagation_path_spectral_radiance_source,
      propagation_path_spectral_radiance_source_jacobian,
      propagation_path_transmission_matrix,
      propagation_path_transmission_matrix_cumulative,
      propagation_path_transmission_matrix_jacobian);
  MuelmatVector background_transmittance;
  background_transmittanceFromPathPropagationBack(
      background_transmittance,
      propagation_path_transmission_matrix_cumulative);
  spectral_radianceFromPathPropagation(spectral_radiance,
                                       propagation_path_spectral_radiance);
  spectral_radiance_jacobianFromBackground(
      spectral_radiance_jacobian,
      spectral_radiance_background_jacobian,
      background_transmittance);
  spectral_radiance_jacobianAddPathPropagation(
      spectral_radiance_jacobian,
      propagation_path_spectral_radiance_jacobian,
      jacobian_targets,
      atmospheric_field,
      propagation_path);
  propagation_path_pointForeground(propagation_path_point, propagation_path);
  spectral_radiance_jacobianApplyUnit(spectral_radiance_jacobian,
                                      spectral_radiance,
                                      frequency_grid,
                                      propagation_path_point,
                                      spectral_radiance_unit);
  spectral_radianceApplyUnit(spectral_radiance,
                             frequency_grid,
                             propagation_path_point,
                             spectral_radiance_unit);
}
ARTS_METHOD_ERROR_CATCH
