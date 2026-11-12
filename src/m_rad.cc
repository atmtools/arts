#include <arts_omp.h>
#include <atm.h>
#include <enums.h>
#include <jac_polyfit.h>
#include <jacobian.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <sensor.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>
#include <exception>

void spectral_rad_jacEmpty(StokvecMatrix         &spectral_rad_jac,
                           const AscendingGrid   &freq_grid,
                           const JacobianTargets &jac_targets) try {
  ARTS_TIME_REPORT

  spectral_rad_jac.resize(jac_targets.x_size(), freq_grid.size());
  spectral_rad_jac = Stokvec{0.0, 0.0, 0.0, 0.0};
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_jacFromBackground(StokvecMatrix             &spectral_rad_jac,
                                    const StokvecMatrix       &spectral_rad_bkg_jac,
                                    const TransmittanceMatrix &spectral_tramat) try {
  ARTS_TIME_REPORT

  const auto [nf, np, nq] = spectral_tramat.shape();

  //! The radiance derivative shape is the background shape
  spectral_rad_jac.resize(spectral_rad_bkg_jac.shape());

  if (nq == 0 or np == 0) return;

  const auto &&background_transmittance = spectral_tramat.P[joker, np - 1];

  //! Set the background radiance derivative as that which is seen after "this" swath
  for (Index i = 0; i < spectral_rad_jac.nrows(); i++) {
    const auto b = spectral_rad_bkg_jac[i];
    auto       s = spectral_rad_jac[i];
    std::transform(
        background_transmittance.begin(), background_transmittance.end(), b.begin(), s.begin(), std::multiplies<>());
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_jacAddPathPropagation(StokvecMatrix                     &spectral_rad_jac,
                                        const StokvecTensor3              &spectral_rad_jac_path,
                                        const JacobianTargets             &jac_targets,
                                        const AtmField                    &atm_field,
                                        const ArrayOfPropagationPathPoint &ray_path) try {
  ARTS_TIME_REPORT

  const Size np = spectral_rad_jac_path.nrows();
  const Size nt = jac_targets.target_count();

  if (nt == 0) return;

  //! The derivative part from the atmosphere
  for (auto &atm_block : jac_targets.atm) {
    if (not atm_field.contains(atm_block.type))
      arts::user_error("No {} in atm_field but in jac_targets", atm_block.type);
    const auto &data = atm_field[atm_block.type];
    for (Size ip = 0; ip < np; ip++) {
      const auto weights = data.flat_weight(ray_path[ip].pos);
      const auto local   = spectral_rad_jac_path[joker, ip, atm_block.target_pos];

      for (auto &w : weights) {
        if (w.second != 0.0) {
          const auto i  = w.first + atm_block.x_start;
          auto       sr = spectral_rad_jac[i];
          std::transform(local.begin(),
                         local.end(),
                         sr.begin(),
                         sr.begin(),
                         [x = w.second](const Stokvec &a, const Stokvec &b) -> Stokvec { return fma(x, a, b); });
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_transform_operatorSet(SpectralRadianceTransformOperator &spectral_rad_transform_operator,
                                        const SpectralRadianceUnitType    &x) try {
  ARTS_TIME_REPORT

  spectral_rad_transform_operator = SpectralRadianceTransformOperator(x);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radApplyUnit(StokvecVector                           &spectral_rad,
                           StokvecMatrix                           &spectral_rad_jac,
                           const AscendingGrid                     &freq_grid,
                           const PropagationPathPoint              &ray_point,
                           const SpectralRadianceTransformOperator &spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  spectral_rad_transform_operator(spectral_rad, spectral_rad_jac, freq_grid, ray_point);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radApplyForwardUnit(StokvecVector                           &spectral_rad,
                                  const AscendingGrid                     &freq_grid,
                                  const PropagationPathPoint              &ray_point,
                                  const SpectralRadianceTransformOperator &spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecMatrix spectral_rad_jac(0, freq_grid.size());

  spectral_radApplyUnit(spectral_rad, spectral_rad_jac, freq_grid, ray_point, spectral_rad_transform_operator);
}
ARTS_METHOD_ERROR_CATCH

namespace {
/** Whether the local workspace of an agenda can be kept from one call to the next
 *
 * An agenda copies the variables that it both reads from the calling workspace
 * and modifies.  A kept local workspace would start the next call from the
 * values that the previous call left in them, so only an agenda without such
 * variables keeps it.
 */
bool keep_local_workspace(const Agenda &agenda) { return agenda.get_copy().empty(); }

/** The radiance and its Jacobian of geometry ip of a sensor, in the unit of the transform
 *
 * lws is the local workspace of spectral_rad_observer_agenda, or nullptr for a
 * new one.  A thread that keeps it from one geometry to the next lets the
 * variables that the agenda computes keep their memory, instead of allocating
 * them anew for every geometry, see keep_local_workspace.
 */
void simulate(const Workspace                         &ws,
              Workspace                               *lws,
              StokvecVector                           &spectral_rad,
              StokvecMatrix                           &spectral_rad_jac,
              const Sensor                            &sensor,
              const Size                               ip,
              const JacobianTargets                   &jac_targets,
              const AtmField                          &atm_field,
              const SurfaceField                      &surf_field,
              const SubsurfaceField                   &subsurf_field,
              const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
              const Agenda                            &spectral_rad_observer_agenda) {
  ArrayOfPropagationPathPoint ray_path;

  Workspace used = spectral_rad_observer_agendaExecute(ws,
                                                       spectral_rad,
                                                       spectral_rad_jac,
                                                       ray_path,
                                                       sensor.freq_grid,
                                                       jac_targets,
                                                       sensor.geometry[ip].pos,
                                                       sensor.geometry[ip].los,
                                                       atm_field,
                                                       surf_field,
                                                       subsurf_field,
                                                       spectral_rad_observer_agenda,
                                                       lws);
  if (lws) *lws = std::move(used);

  if (ray_path.empty()) arts::user_error("No ray path found");
  spectral_rad_transform_operator(spectral_rad, spectral_rad_jac, sensor.freq_grid, ray_path.front());
}

void low_memory(const Workspace                         &ws,
                Vector                                  &measurement_vec,
                Matrix                                  &measurement_jac,
                const ArrayOfSensor                     &measurement_sensors,
                const JacobianTargets                   &jac_targets,
                const AtmField                          &atm_field,
                const SurfaceField                      &surf_field,
                const SubsurfaceField                   &subsurf_field,
                const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
                const Agenda                            &spectral_rad_observer_agenda) {
  ARTS_TIME_REPORT

  Size N = 0;
  for (const auto &sensor : measurement_sensors) N += sensor.geometry.size();

  std::string error{};

#pragma omp parallel if (arts_omp_parallel(-1, N > 1))
  {
    Workspace  lws{WorkspaceInitialization::Empty};
    Workspace *local = keep_local_workspace(spectral_rad_observer_agenda) ? &lws : nullptr;
    Size       start = 0;
    for (Size isensor = 0; isensor < measurement_sensors.size(); isensor++) {
      const auto &sensor = measurement_sensors[isensor];

      //! Dynamic scheduling as some simulations may take much longer time than others
#pragma omp for schedule(dynamic) nowait
      for (Size ip = 0; ip < sensor.geometry.size(); ip++) {
        try {
          StokvecVector spectral_rad;
          StokvecMatrix spectral_rad_jac;

          simulate(ws,
                   local,
                   spectral_rad,
                   spectral_rad_jac,
                   sensor,
                   ip,
                   jac_targets,
                   atm_field,
                   surf_field,
                   subsurf_field,
                   spectral_rad_transform_operator,
                   spectral_rad_observer_agenda);

#pragma omp critical
          for (Size ie = 0; ie < sensor.elements.size(); ++ie) {
            measurement_vec[start + ie] += sensor::sumup(sensor.elements[ie], spectral_rad, ip);

            sensor::sumup(measurement_jac[start + ie], sensor.elements[ie], spectral_rad_jac, ip);
          }
        } catch (const std::exception &e) {
#pragma omp critical
          if (error.empty())
            error = std::format(R"(Sensor "{}", geometry {}:
{})",
                                sensor.name(),
                                ip,
                                e.what());
        }
      }

      start += sensor.elements.size();
    }
  }

  if (not error.empty()) arts::user_error("Errors occurred:\n{:}", error);
}

void high_performance(const Workspace                         &ws,
                      Vector                                  &measurement_vec,
                      Matrix                                  &measurement_jac,
                      const ArrayOfSensor                     &measurement_sensors,
                      const JacobianTargets                   &jac_targets,
                      const AtmField                          &atm_field,
                      const SurfaceField                      &surf_field,
                      const SubsurfaceField                   &subsurf_field,
                      const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
                      const Agenda                            &spectral_rad_observer_agenda) {
  ARTS_TIME_REPORT

  const Size M = measurement_vec.size();
  const Size J = jac_targets.x_size();

  //! The radiance of every geometry of every sensor
  std::vector<std::vector<std::pair<StokvecVector, StokvecMatrix>>> cache(measurement_sensors.size());
  Size                                                              N = 0;
  for (Size isensor = 0; isensor < measurement_sensors.size(); isensor++) {
    cache[isensor].resize(measurement_sensors[isensor].geometry.size());
    N += measurement_sensors[isensor].geometry.size();
  }

  std::string error{};

#pragma omp parallel if (arts_omp_parallel(-1, N > 1))
  {
    Workspace  lws{WorkspaceInitialization::Empty};
    Workspace *local = keep_local_workspace(spectral_rad_observer_agenda) ? &lws : nullptr;
    for (Size isensor = 0; isensor < measurement_sensors.size(); isensor++) {
      const auto &sensor = measurement_sensors[isensor];

      //! Dynamic scheduling as some simulations may take much longer time than others
#pragma omp for schedule(dynamic) nowait
      for (Size ip = 0; ip < sensor.geometry.size(); ip++) {
        try {
          auto &[spectral_rad, spectral_rad_jac] = cache[isensor][ip];

          simulate(ws,
                   local,
                   spectral_rad,
                   spectral_rad_jac,
                   sensor,
                   ip,
                   jac_targets,
                   atm_field,
                   surf_field,
                   subsurf_field,
                   spectral_rad_transform_operator,
                   spectral_rad_observer_agenda);
        } catch (const std::exception &e) {
#pragma omp critical
          if (error.empty())
            error = std::format(R"(Sensor "{}", geometry {}:
{})",
                                sensor.name(),
                                ip,
                                e.what());
        }
      }
    }
  }

  if (not error.empty()) arts::user_error("Errors occurred:\n{:}", error);

#pragma omp parallel if (arts_omp_parallel(-1, M > 1))
  {
    Size start = 0;
    for (Size isensor = 0; isensor < measurement_sensors.size(); isensor++) {
      const auto &sensor = measurement_sensors[isensor];
      const auto &data   = cache[isensor];

#pragma omp for schedule(dynamic) nowait
      for (Size ie = 0; ie < sensor.elements.size(); ++ie) {
        const Size iv = start + ie;

        for (auto &sparse_weights : sensor.elements[ie]) {
          const Size ip    = sparse_weights.irow;
          const Size ifreq = sparse_weights.icol;

          const auto &[spectral_rad, spectral_rad_jac] = data[ip];

          measurement_vec[iv] += dot(sparse_weights.data, spectral_rad[ifreq]);

          auto jac   = measurement_jac[iv];
          auto srjac = spectral_rad_jac[joker, ifreq];
          for (Size ij = 0; ij < J; ij++) { jac[ij] += dot(sparse_weights.data, srjac[ij]); }
        }
      }

      start += sensor.elements.size();
    }
  }
}

void indexed(const Workspace                         &ws,
             Vector                                  &measurement_vec,
             Matrix                                  &measurement_jac,
             const ArrayOfSensor                     &measurement_sensors,
             const JacobianTargets                   &jac_targets,
             const AtmField                          &atm_field,
             const SurfaceField                      &surf_field,
             const SubsurfaceField                   &subsurf_field,
             const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
             const Agenda                            &spectral_rad_observer_agenda) {
  ARTS_TIME_REPORT

  const Size J = jac_targets.x_size();

  //! The elements of each sensor that have weights at each of its geometries
  std::vector<sensor::GeometryIndex> indices;
  indices.reserve(measurement_sensors.size());
  Size N = 0;
  for (const auto &sensor : measurement_sensors) {
    indices.emplace_back(sensor);
    N += sensor.geometry.size();
  }

  std::string error{};

#pragma omp parallel if (arts_omp_parallel(-1, N > 1))
  {
    Workspace  lws{WorkspaceInitialization::Empty};
    Workspace *local = keep_local_workspace(spectral_rad_observer_agenda) ? &lws : nullptr;
    Size       start = 0;
    for (Size isensor = 0; isensor < measurement_sensors.size(); isensor++) {
      const auto &sensor = measurement_sensors[isensor];
      const auto &index  = indices[isensor];

      //! Dynamic scheduling as some simulations may take much longer time than others
#pragma omp for schedule(dynamic) nowait
      for (Size ip = 0; ip < sensor.geometry.size(); ip++) {
        try {
          StokvecVector spectral_rad;
          StokvecMatrix spectral_rad_jac;

          simulate(ws,
                   local,
                   spectral_rad,
                   spectral_rad_jac,
                   sensor,
                   ip,
                   jac_targets,
                   atm_field,
                   surf_field,
                   subsurf_field,
                   spectral_rad_transform_operator,
                   spectral_rad_observer_agenda);

          //! An element with weights at this geometry only is written by this thread only
          Size nshared = 0;
          for (const auto &[ie, weights] : index[ip]) {
            if (index.geometry_count(ie) > 1) {
              nshared++;
              continue;
            }

            measurement_vec[start + ie] += sensor::sumup(weights, spectral_rad);
            sensor::sumup(measurement_jac[start + ie], weights, spectral_rad_jac);
          }

          if (nshared == 0) continue;

          //! The other elements are summed here, and only the adding to them is locked
          Vector shared_vec(nshared, 0.0);
          Matrix shared_jac(nshared, J, 0.0);
          Size   ishared = 0;
          for (const auto &[ie, weights] : index[ip]) {
            if (index.geometry_count(ie) == 1) continue;

            shared_vec[ishared] = sensor::sumup(weights, spectral_rad);
            sensor::sumup(shared_jac[ishared], weights, spectral_rad_jac);
            ishared++;
          }

#pragma omp critical
          {
            ishared = 0;
            for (const auto &entry : index[ip]) {
              if (index.geometry_count(entry.element) == 1) continue;

              measurement_vec[start + entry.element] += shared_vec[ishared];
              measurement_jac[start + entry.element] += shared_jac[ishared];
              ishared++;
            }
          }
        } catch (const std::exception &e) {
#pragma omp critical
          if (error.empty())
            error = std::format(R"(Sensor "{}", geometry {}:
{})",
                                sensor.name(),
                                ip,
                                e.what());
        }
      }

      start += sensor.elements.size();
    }
  }

  if (not error.empty()) arts::user_error("Errors occurred:\n{:}", error);
}
}  // namespace

void measurement_vecFromSensor(const Workspace                         &ws,
                               Vector                                  &measurement_vec,
                               Matrix                                  &measurement_jac,
                               const ArrayOfSensor                     &measurement_sensors,
                               const JacobianTargets                   &jac_targets,
                               const AtmField                          &atm_field,
                               const SurfaceField                      &surf_field,
                               const SubsurfaceField                   &subsurf_field,
                               const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
                               const Agenda                            &spectral_rad_observer_agenda,
                               const String                            &kernel) try {
  ARTS_TIME_REPORT

  const Size ny = measurement_size(measurement_sensors);

  measurement_vec.resize(ny);
  measurement_vec = 0.0;

  measurement_jac.resize(ny, jac_targets.x_size());
  measurement_jac = 0.0;

  if (ny == 0) return;

  for (auto &sensor : measurement_sensors) sensor.check();

  switch (to<MeasurementVectorSumupKernel>(kernel)) {
    using enum MeasurementVectorSumupKernel;
    case LowMem:
      low_memory(ws,
                 measurement_vec,
                 measurement_jac,
                 measurement_sensors,
                 jac_targets,
                 atm_field,
                 surf_field,
                 subsurf_field,
                 spectral_rad_transform_operator,
                 spectral_rad_observer_agenda);
      break;
    case HighPerf:
      high_performance(ws,
                       measurement_vec,
                       measurement_jac,
                       measurement_sensors,
                       jac_targets,
                       atm_field,
                       surf_field,
                       subsurf_field,
                       spectral_rad_transform_operator,
                       spectral_rad_observer_agenda);
      break;
    case Indexed:
      indexed(ws,
              measurement_vec,
              measurement_jac,
              measurement_sensors,
              jac_targets,
              atm_field,
              surf_field,
              subsurf_field,
              spectral_rad_transform_operator,
              spectral_rad_observer_agenda);
      break;
  }
}
ARTS_METHOD_ERROR_CATCH

void measurement_jacAddSensorJacobianPerturbations(const Workspace                         &ws,
                                                   Matrix                                  &measurement_jac,
                                                   const Vector                            &measurement_vec,
                                                   const ArrayOfSensor                     &measurement_sensors,
                                                   const JacobianTargets                   &jac_targets,
                                                   const AtmField                          &atm_field,
                                                   const SurfaceField                      &surf_field,
                                                   const SubsurfaceField                   &subsurf_field,
                                                   const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
                                                   const Agenda                            &spectral_rad_observer_agenda,
                                                   const String                            &kernel) try {
  ARTS_TIME_REPORT

  const Size ny = measurement_size(measurement_sensors);
  if (static_cast<Size>(measurement_vec.size()) != ny)
    arts::user_error("measurement_vec has {} elements but measurement_sensors has {}", measurement_vec.size(), ny);
  if (measurement_jac.shape() != std::array{static_cast<Index>(ny), static_cast<Index>(jac_targets.x_size())})
    arts::user_error("measurement_jac must have shape [{}, {}], got {:B,}", ny, jac_targets.x_size(), measurement_jac.shape());

  const JacobianTargets no_targets{};

  for (const auto &target : jac_targets.sensor) {
    if (target.type.sensor_elem < 0 or static_cast<Size>(target.type.sensor_elem) >= measurement_sensors.size())
      arts::user_error("Bad sensor {}, out-of-bounds", target);

    const Size    isensor = static_cast<Size>(target.type.sensor_elem);
    const Sensor &sensor  = measurement_sensors[isensor];

    //! The derivative of the sensor grid with regards to the coefficients of the polynomial offset
    const auto *poly = target.inverse_state.target<polyinv_sensor_offset_t>();
    if (poly == nullptr)
      arts::user_error(R"(Sensor "{}": target {} has no polynomial offset model)", sensor.name(), target.type);
    const Matrix dgrid = poly->normalized()();
    const Vector grid  = sensor.flat(target.type.type);
    if (static_cast<Size>(dgrid.nrows()) != grid.size() or static_cast<Size>(dgrid.ncols()) != target.x_size)
      arts::user_error(
          R"(Sensor "{}": target {} has a polynomial offset model of shape {:B,} for {} grid points and {} coefficients)",
          sensor.name(),
          target.type,
          dgrid.shape(),
          grid.size(),
          target.x_size);

    const Size    start = measurement_start(measurement_sensors, isensor);
    const Numeric d     = target.d;

    for (Size k = 0; k < target.x_size; k++) {
      ArrayOfSensor perturbed{sensor};

      Vector x = grid;
      for (Size i = 0; i < x.size(); i++) x[i] += d * dgrid[i, k];
      perturbed.front().unflatten(x, target.type.type);

      Vector y;
      Matrix jac;
      measurement_vecFromSensor(ws,
                                y,
                                jac,
                                perturbed,
                                no_targets,
                                atm_field,
                                surf_field,
                                subsurf_field,
                                spectral_rad_transform_operator,
                                spectral_rad_observer_agenda,
                                kernel);

      for (Size ie = 0; ie < y.size(); ie++) {
        measurement_jac[start + ie, target.x_start + k] += (y[ie] - measurement_vec[start + ie]) / d;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
