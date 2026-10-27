#include <arts_constants.h>
#include <jacobian.h>
#include <path_point.h>
#include <radar.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <sensor_meta_info.h>
#include <time_report.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

namespace {
struct Settings {
  Stokvec             transmitted;
  RadarRangeMode      range_mode;
  RadarUnit           unit;
  Numeric             ze_tref;
  Numeric             k2;
  Numeric             dbze_min;
  Numeric             pext_scaling;
  TransmittanceOption rte_option;
};

void check_settings(const ArrayOfSensorObsel&       sensor,
                    const Matrix&                   limits,
                    const JacobianTargets&          jac_targets,
                    const ArrayOfScatteringSpecies& scattering_species,
                    const Settings&                 settings) {
  const Size ny = sensor.size();
  if (ny == 0) arts::user_error("measurement_sensor cannot be empty for a radar calculation");
  if (limits.nrows() != static_cast<Index>(ny) or limits.ncols() != 2)
    arts::user_error("radar_range_limits must have shape [{}, 2], got {:B,}", ny, limits.shape());
  for (Size iy = 0; iy < ny; ++iy) {
    sensor[iy].check();
    if (not(limits[iy, 0] < limits[iy, 1]))
      arts::user_error(
          "Radar range limits in row {} are not strictly increasing: [{}, {}]", iy, limits[iy, 0], limits[iy, 1]);
  }
  if (settings.transmitted.I() != 1.0)
    arts::user_error("The transmitted Stokes I component must equal one, got {}", settings.transmitted.I());
  if (settings.pext_scaling < 0.0 or settings.pext_scaling > 2.0)
    arts::user_error("pext_scaling must be between 0 and 2; got {}", settings.pext_scaling);
  if (scattering_species.species.empty()) arts::user_error("Radar single scattering requires scattering species");
  if (not jac_targets.surf.empty() or not jac_targets.subsurf.empty() or not jac_targets.line.empty() or
      not jac_targets.sensor.empty() or not jac_targets.error.empty())
    arts::user_error("The radar forward model currently supports atmospheric Jacobian targets only");
}

RadarRangeMode resolved_range_mode(const RadarRangeMode mode, const Matrix& limits) {
  Numeric largest = -std::numeric_limits<Numeric>::infinity();
  for (Index i = 0; i < limits.nrows(); ++i)
    for (Index j = 0; j < limits.ncols(); ++j) largest = std::max(largest, limits[i, j]);
  return radar::resolve(mode, largest);
}

//! A sensor geometry, its path, and the samples of the gates of all sensor elements that use it
struct PreparedSimulation {
  const AscendingGrid*               frequency{};
  Size                               iposlos{};
  ArrayOfPropagationPathPoint        ray_path{};
  radar::PathRange                   range{};
  std::vector<Size>                  gates{};
  std::vector<rtepack::radar_sample> samples{};
  Size                               np{};
};

std::vector<PreparedSimulation> prepare_simulations(const Workspace&          ws,
                                                    const ArrayOfSensorObsel& sensor,
                                                    const Matrix&             limits,
                                                    const RadarRangeMode      range_mode,
                                                    const SurfaceField&       surf_field,
                                                    const Agenda&             ray_path_observer_agenda) {
  std::vector<PreparedSimulation> out;
  for (const auto& [frequency, poslos, iposlos] : collect_simulations(sensor)) {
    PreparedSimulation simulation{.frequency = &frequency, .iposlos = iposlos};
    ray_path_observer_agendaExecute(
        ws, simulation.ray_path, poslos[iposlos].pos, poslos[iposlos].los, ray_path_observer_agenda);
    simulation.range = radar::path_range(simulation.ray_path, poslos[iposlos].pos, surf_field.ellipsoid);

    const Vector coordinate = radar::range_coordinate(simulation.range, simulation.ray_path, range_mode);
    for (Size iy = 0; iy < sensor.size(); ++iy) {
      if (not sensor[iy].same_freqs(frequency) or not sensor[iy].same_poslos(poslos)) continue;
      radar::add_gate_samples(simulation.samples, coordinate, limits[iy, 0], limits[iy, 1], simulation.gates.size());
      simulation.gates.push_back(iy);
    }
    simulation.np = radar::sort_samples(simulation.samples);

    out.push_back(std::move(simulation));
  }
  return out;
}

struct RadarResult {
  Vector measurement;
  Matrix jacobian;
  Matrix auxiliary;
  Vector time;
};

RadarResult radar_forward(const Workspace&                       ws,
                          const ArrayOfSensorObsel&              sensor,
                          const std::vector<PreparedSimulation>& simulations,
                          const AtmField&                        atm_field,
                          const ArrayOfScatteringSpecies&        scattering_species,
                          const Agenda&                          spectral_propmat_agenda,
                          const Settings&                        settings,
                          const ArrayOfRadarAuxiliary&           aux_vars,
                          const JacobianTargets&                 jac_targets) {
  const Size        ny = sensor.size();
  RadarResult       out{.measurement = Vector(ny, 0.0),
                        .jacobian    = Matrix(ny, jac_targets.x_size(), 0.0),
                        .auxiliary   = Matrix(aux_vars.size(), ny, 0.0),
                        .time        = Vector(ny, 0.0)};
  std::vector<bool> covered(ny, false);
  Vector            scalar_weight_sum(ny, 0.0);

  for (const auto& simulation : simulations) {
    if (simulation.samples.empty()) continue;

    const auto& frequency = *simulation.frequency;
    const auto& ray_path  = simulation.ray_path;
    const auto& samples   = simulation.samples;
    const Size  nf        = frequency.size();
    const Size  ngate     = simulation.gates.size();
    const Size  iposlos   = simulation.iposlos;

    const radar::GasModel gas = [&](PropmatVector&              propmat,
                                    PropmatMatrix&              propmat_jac,
                                    const PropagationPathPoint& point,
                                    const AtmPoint&             atm_point) {
      StokvecVector srcvec;
      StokvecMatrix srcvec_jac;
      spectral_propmat_agendaExecute(ws,
                                     propmat,
                                     srcvec,
                                     propmat_jac,
                                     srcvec_jac,
                                     frequency,
                                     Vector3{0.0, 0.0, 0.0},
                                     jac_targets,
                                     {},
                                     point,
                                     atm_point,
                                     spectral_propmat_agenda);
    };
    const radar::Optics optics = radar::optics(
        ray_path, simulation.np, atm_field, frequency, scattering_species, jac_targets, settings.pext_scaling, gas);

    StokvecMatrix  signal;
    StokvecTensor4 signal_jac;
    radar::single_scattering(
        signal, signal_jac, settings.transmitted, optics, simulation.range.layer, samples, ngate, settings.rte_option);

    Vector factor(nf);
    for (Size iv = 0; iv < nf; ++iv)
      factor[iv] = radar::unit_factor(settings.unit, frequency[iv], settings.ze_tref, settings.k2);

    const Vector  coverage   = radar::gate_coverage(samples, ngate);
    const Vector  time       = radar::gate_sum(samples, simulation.range.time, ngate);
    const Numeric background = ray_path.back().los_type == PathPositionType::space
                                   ? 0.0
                                   : (ray_path.back().los_type == PathPositionType::surface ? 1.0 : 2.0);

    const auto contains = [&](const RadarAuxiliary key) { return stdr::contains(aux_vars, key); };

    StokvecMatrix backscatter(ngate, nf, Stokvec{});
    if (contains(RadarAuxiliary::Backscattering)) {
      for (const auto& [layer, x, weight, gate] : samples)
        for (Size iv = 0; iv < nf; ++iv)
          backscatter[gate, iv] +=
              weight * (((1.0 - x) * optics.Z[layer - 1][iv] + x * optics.Z[layer][iv]) * settings.transmitted);
    }

    Matrix gas_extinction(ngate, nf, 0.0);
    if (contains(RadarAuxiliary::AbsSpeciesExtinction))
      for (Size iv = 0; iv < nf; ++iv)
        gas_extinction[joker, iv] = radar::gate_sum(samples, optics.gas_extinction[iv], ngate);

    Matrix particle_extinction(ngate, nf, 0.0);
    if (contains(RadarAuxiliary::ParticleExtinction))
      for (Size iv = 0; iv < nf; ++iv)
        particle_extinction[joker, iv] = radar::gate_sum(samples, optics.particle_extinction[iv], ngate);

    StokvecVector bin(nf);
    for (Size ig = 0; ig < ngate; ++ig) {
      if (coverage[ig] == 0.0) continue;

      const Size  iy    = simulation.gates[ig];
      const auto& obsel = sensor[iy];
      covered[iy]       = true;

      for (Size iv = 0; iv < nf; ++iv) bin[iv] = factor[iv] * signal[ig, iv];
      out.measurement[iy] += obsel.sumup(bin, iposlos);

      // Extinction, background and time are means over the covered part of the gate
      const Numeric scalar_weight  = obsel.sumup(Stokvec{1.0, 0.0, 0.0, 0.0}, iposlos);
      scalar_weight_sum[iy]       += scalar_weight;
      out.time[iy]                += scalar_weight * time[ig] / coverage[ig];
      for (Size ia = 0; ia < aux_vars.size(); ++ia) {
        switch (aux_vars[ia]) {
          case RadarAuxiliary::RadiativeBackground: out.auxiliary[ia, iy] += scalar_weight * background; break;
          case RadarAuxiliary::Backscattering:
            for (Size iv = 0; iv < nf; ++iv) bin[iv] = factor[iv] * backscatter[ig, iv];
            out.auxiliary[ia, iy] += obsel.sumup(bin, iposlos);
            break;
          case RadarAuxiliary::AbsSpeciesExtinction:
            for (Size iv = 0; iv < nf; ++iv) bin[iv] = Stokvec{gas_extinction[ig, iv] / coverage[ig]};
            out.auxiliary[ia, iy] += obsel.sumup(bin, iposlos);
            break;
          case RadarAuxiliary::ParticleExtinction:
            for (Size iv = 0; iv < nf; ++iv) bin[iv] = Stokvec{particle_extinction[ig, iv] / coverage[ig]};
            out.auxiliary[ia, iy] += obsel.sumup(bin, iposlos);
            break;
        }
      }
    }

    for (const auto& target : jac_targets.atm) {
      const auto& field = atm_field[target.type];
      for (Size jp = 0; jp < simulation.np; ++jp) {
        const auto weights = field.flat_weight(ray_path[jp].pos);
        for (Size ig = 0; ig < ngate; ++ig) {
          if (coverage[ig] == 0.0) continue;

          const Size iy = simulation.gates[ig];
          for (Size iv = 0; iv < nf; ++iv) bin[iv] = factor[iv] * signal_jac[ig, jp, target.target_pos, iv];
          const Numeric local = sensor[iy].sumup(bin, iposlos);
          for (const auto& [index, weight] : weights) out.jacobian[iy, target.x_start + index] += weight * local;
        }
      }
    }
  }

  constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();
  for (Size iy = 0; iy < ny; ++iy) {
    if (not covered[iy]) {
      out.measurement[iy]      = nan;
      out.time[iy]             = nan;
      out.auxiliary[joker, iy] = nan;
      continue;
    }

    if (scalar_weight_sum[iy] == 0.0) continue;
    out.time[iy] /= scalar_weight_sum[iy];
    for (Size ia = 0; ia < aux_vars.size(); ++ia)
      if (aux_vars[ia] != RadarAuxiliary::Backscattering) out.auxiliary[ia, iy] /= scalar_weight_sum[iy];
  }

  if (settings.unit == RadarUnit::dBZe) {
    const Numeric ze_min = std::pow(10.0, settings.dbze_min / 10.0);
    for (Size iy = 0; iy < ny; ++iy) {
      if (std::isnan(out.measurement[iy])) continue;
      if (out.measurement[iy] <= ze_min) {
        out.measurement[iy]     = settings.dbze_min;
        out.jacobian[iy, joker] = 0.0;
      } else {
        out.jacobian[iy, joker] *= 10.0 / (std::log(10.0) * out.measurement[iy]);
        out.measurement[iy]      = 10.0 * std::log10(out.measurement[iy]);
      }
    }
    for (Size ia = 0; ia < aux_vars.size(); ++ia) {
      if (aux_vars[ia] != RadarAuxiliary::Backscattering) continue;
      for (Size iy = 0; iy < ny; ++iy) {
        auto& value = out.auxiliary[ia, iy];
        if (not std::isnan(value)) value = value <= ze_min ? settings.dbze_min : 10.0 * std::log10(value);
      }
    }
  }

  return out;
}

//! Apply the target transformations to a Jacobian with regards to the untransformed model state
void transform_jacobian(Matrix& jacobian, const JacobianTargets& jac_targets, const AtmField& atm_field) {
  const Size nx = jac_targets.x_size();
  if (nx == 0) return;

  Vector state(nx, 0.0);
  for (const auto& target : jac_targets.atm) target.update_state(state, atm_field);
  for (const auto& target : jac_targets.atm) target.update_jac(jacobian, state, atm_field);
}
}  // namespace

void measurement_sensorAddSimpleRadar(ArrayOfSensorObsel&    measurement_sensor,
                                      ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                      Matrix&                radar_range_limits,
                                      const AscendingGrid&   freq_grid,
                                      const Vector3&         pos,
                                      const Vector2&         los,
                                      const Stokvec&         pol,
                                      const AscendingGrid&   range_bins) try {
  ARTS_TIME_REPORT
  if (freq_grid.empty()) arts::user_error("Radar frequency grid cannot be empty");
  if (range_bins.size() < 2) arts::user_error("Radar range_bins needs at least two edges");

  const Size old_n = measurement_sensor.size();
  if (radar_range_limits.nrows() != 0 and
      (radar_range_limits.ncols() != 2 or static_cast<Size>(radar_range_limits.nrows()) != old_n))
    arts::user_error(
        "radar_range_limits must be empty or have shape [{}, 2], got {:B,}", old_n, radar_range_limits.shape());

  const Size nbins = range_bins.size() - 1;
  const Size nnew  = freq_grid.size() * nbins;
  Matrix     new_limits(old_n + nnew, 2);
  if (old_n) new_limits[Range(0, old_n), joker] = radar_range_limits;

  auto poslos = std::make_shared<const SensorPosLosVector>(SensorPosLosVector{SensorPosLos{.pos = pos, .los = los}});
  measurement_sensor.reserve(old_n + nnew);
  measurement_sensor_meta.reserve(measurement_sensor_meta.size() + freq_grid.size());

  Size row = old_n;
  for (Size iv = 0; iv < freq_grid.size(); ++iv) {
    auto                        one_frequency = std::make_shared<const AscendingGrid>(AscendingGrid{freq_grid[iv]});
    sensor::SparseStokvecMatrix weights(1, 1);
    weights[0, 0] = pol;

    for (Size ib = 0; ib < nbins; ++ib) {
      measurement_sensor.emplace_back(one_frequency, poslos, weights);
      new_limits[row, 0] = range_bins[ib];
      new_limits[row, 1] = range_bins[ib + 1];
      ++row;
    }

    SortedGriddedField1 meta;
    meta.data_name     = std::format("radar-{}-Hz", freq_grid[iv]);
    meta.gridname<0>() = "range";
    Vector centers(nbins);
    for (Size ib = 0; ib < nbins; ++ib) centers[ib] = std::midpoint(range_bins[ib], range_bins[ib + 1]);
    meta.grid<0>() = AscendingGrid{std::move(centers)};
    meta.data.resize(nbins);
    meta.data = 0.0;
    measurement_sensor_meta.push_back(SensorMetaInfo{.data = std::move(meta)});
  }
  radar_range_limits = std::move(new_limits);
}
ARTS_METHOD_ERROR_CATCH

void measurement_vecFromRadarSingleScattering(const Workspace&                ws,
                                              Vector&                         measurement_vec,
                                              Matrix&                         measurement_jac,
                                              Matrix&                         radar_aux,
                                              Vector&                         radar_time,
                                              const ArrayOfSensorObsel&       measurement_sensor,
                                              const Matrix&                   radar_range_limits,
                                              const JacobianTargets&          jac_targets,
                                              const AtmField&                 atm_field,
                                              const SurfaceField&             surf_field,
                                              const ArrayOfScatteringSpecies& scattering_species,
                                              const Agenda&                   ray_path_observer_agenda,
                                              const Agenda&                   spectral_propmat_agenda,
                                              const TransmittanceOption&      rte_option,
                                              const Stokvec&                  transmitted_stokes,
                                              const RadarRangeMode&           range_mode,
                                              const RadarUnit&                unit,
                                              const Numeric&                  ze_tref,
                                              const Numeric&                  k2,
                                              const Numeric&                  dbze_min,
                                              const Numeric&                  pext_scaling,
                                              const ArrayOfRadarAuxiliary&    aux_vars) try {
  ARTS_TIME_REPORT

  const Settings settings{.transmitted  = transmitted_stokes,
                          .range_mode   = resolved_range_mode(range_mode, radar_range_limits),
                          .unit         = unit,
                          .ze_tref      = ze_tref,
                          .k2           = k2,
                          .dbze_min     = dbze_min,
                          .pext_scaling = pext_scaling,
                          .rte_option   = rte_option};
  check_settings(measurement_sensor, radar_range_limits, jac_targets, scattering_species, settings);

  const auto simulations = prepare_simulations(
      ws, measurement_sensor, radar_range_limits, settings.range_mode, surf_field, ray_path_observer_agenda);
  RadarResult result = radar_forward(ws,
                                     measurement_sensor,
                                     simulations,
                                     atm_field,
                                     scattering_species,
                                     spectral_propmat_agenda,
                                     settings,
                                     aux_vars,
                                     jac_targets);

  measurement_vec = std::move(result.measurement);
  measurement_jac = std::move(result.jacobian);
  radar_aux       = std::move(result.auxiliary);
  radar_time      = std::move(result.time);
  transform_jacobian(measurement_jac, jac_targets, atm_field);
}
ARTS_METHOD_ERROR_CATCH

void model_state_vecFromRadarOnionPeeling(const Workspace&                ws,
                                          Vector&                         model_state_vec,
                                          Vector&                         measurement_vec_fit,
                                          Matrix&                         measurement_jac,
                                          AtmField&                       atm_field,
                                          const Vector&                   measurement_vec,
                                          const ArrayOfSensorObsel&       measurement_sensor,
                                          const Matrix&                   radar_range_limits,
                                          const JacobianTargets&          jac_targets,
                                          const SurfaceField&             surf_field,
                                          const ArrayOfScatteringSpecies& scattering_species,
                                          const Agenda&                   ray_path_observer_agenda,
                                          const Agenda&                   spectral_propmat_agenda,
                                          const TransmittanceOption&      rte_option,
                                          const Stokvec&                  transmitted_stokes,
                                          const RadarRangeMode&           range_mode,
                                          const RadarUnit&                unit,
                                          const Numeric&                  ze_tref,
                                          const Numeric&                  k2,
                                          const Numeric&                  dbze_min,
                                          const Numeric&                  pext_scaling,
                                          const Numeric&                  measurement_noise_floor,
                                          const Numeric&                  state_min,
                                          const Numeric&                  state_max,
                                          const Numeric&                  max_step,
                                          const Numeric&                  tolerance,
                                          const Index&                    max_iterations,
                                          const Index&                    max_sweeps) try {
  ARTS_TIME_REPORT
  const Size ny = measurement_sensor.size();
  const Size nx = jac_targets.x_size();

  const Settings settings{.transmitted  = transmitted_stokes,
                          .range_mode   = resolved_range_mode(range_mode, radar_range_limits),
                          .unit         = unit,
                          .ze_tref      = ze_tref,
                          .k2           = k2,
                          .dbze_min     = dbze_min,
                          .pext_scaling = pext_scaling,
                          .rte_option   = rte_option};
  check_settings(measurement_sensor, radar_range_limits, jac_targets, scattering_species, settings);
  if (measurement_vec.size() != ny)
    arts::user_error("measurement_vec has {} elements but there are {} radar observations", measurement_vec.size(), ny);
  if (nx == 0 or jac_targets.atm.empty())
    arts::user_error("Radar onion peeling requires at least one finalized atmospheric Jacobian target");
  if (state_min > state_max) arts::user_error("state_min ({}) exceeds state_max ({})", state_min, state_max);
  if (max_step <= 0.0) arts::user_error("max_step must be positive, got {}", max_step);
  if (tolerance <= 0.0) arts::user_error("tolerance must be positive, got {}", tolerance);
  if (max_iterations < 1) arts::user_error("max_iterations must be positive, got {}", max_iterations);
  if (max_sweeps < 1) arts::user_error("max_sweeps must be positive, got {}", max_sweeps);

  // The state is deliberately obtained and applied through the normal target
  // transformations.  The peeled profile is therefore an ordinary ARTS model
  // state and can be passed directly to OEM (including transformed targets).
  model_state_vec.resize(nx);
  model_state_vec = 0.0;
  for (const auto& target : jac_targets.atm) target.update_state(model_state_vec, atm_field);

  // The paths are traced for every update, since the path agenda may depend
  // on the atmosphere that is being peeled.
  const auto calculate = [&]() {
    const auto simulations = prepare_simulations(
        ws, measurement_sensor, radar_range_limits, settings.range_mode, surf_field, ray_path_observer_agenda);
    RadarResult result  = radar_forward(ws,
                                       measurement_sensor,
                                       simulations,
                                       atm_field,
                                       scattering_species,
                                       spectral_propmat_agenda,
                                       settings,
                                        {},
                                       jac_targets);
    measurement_vec_fit = std::move(result.measurement);
    measurement_jac     = std::move(result.jacobian);
    transform_jacobian(measurement_jac, jac_targets, atm_field);
  };

  calculate();

  // Onion order is distance from the sensor.  Distance and round-trip time
  // already increase away from it.  For altitude gates, use the altitude of
  // the observation geometry so both downward- and upward-looking radars are
  // ordered correctly.  This ordering is local to every sensor geometry;
  // interleaving independent profiles is harmless because their Jacobian
  // blocks have zero cross-sensitivity.
  std::vector<Size> order(ny);
  stdr::iota(order, Size{0});
  const auto distance_key = [&](const Size iy) {
    const Numeric center = std::midpoint(radar_range_limits[iy, 0], radar_range_limits[iy, 1]);
    if (settings.range_mode != RadarRangeMode::Altitude) return center;
    const auto& poslos = measurement_sensor[iy].poslos_grid();
    if (poslos.empty()) arts::user_error("Radar observation {} has no position/LOS geometry", iy);
    return std::abs(center - poslos.front().pos[0]);
  };
  stdr::stable_sort(order, [&](const Size a, const Size b) { return distance_key(a) < distance_key(b); });

  std::vector<bool>                  peeled(nx, false);
  std::vector<std::pair<Size, Size>> gate_state;
  for (const Size iy : order) {
    if (not std::isfinite(measurement_vec[iy]) or measurement_vec[iy] <= measurement_noise_floor) continue;

    // Associate this range gate with its strongest not-yet-retrieved local
    // atmospheric state sensitivity.  This uses the same mapped Jacobian as
    // OEM and avoids imposing a legacy pressure-grid or two-species layout.
    Size    ix_best = nx;
    Numeric best    = 0.0;
    for (Size ix = 0; ix < nx; ++ix) {
      const Numeric sensitivity = std::abs(measurement_jac[iy, ix]);
      if (not peeled[ix] and std::isfinite(sensitivity) and sensitivity > best) {
        best    = sensitivity;
        ix_best = ix;
      }
    }
    if (ix_best == nx) continue;

    bool converged = false;
    for (Index iteration = 0; iteration < max_iterations; ++iteration) {
      const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
      if (std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]))) {
        converged = true;
        break;
      }

      const Numeric derivative = measurement_jac[iy, ix_best];
      if (not std::isfinite(derivative) or derivative == 0.0)
        arts::user_error(
            "Radar onion peeling encountered a zero/non-finite derivative for gate {} and state {}", iy, ix_best);
      Numeric step             = residual / derivative;
      step                     = std::clamp(step, -max_step, max_step);
      model_state_vec[ix_best] = std::clamp(model_state_vec[ix_best] + step, state_min, state_max);
      for (const auto& target : jac_targets.atm) target.update_model(atm_field, model_state_vec);
      calculate();
    }
    if (not converged) {
      const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
      converged = std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]));
    }
    if (not converged)
      arts::user_error("Radar onion peeling did not converge for gate {} (observed {}, fitted {}) after {} iterations",
                       iy,
                       measurement_vec[iy],
                       measurement_vec_fit[iy],
                       max_iterations);
    peeled[ix_best] = true;
    gate_state.emplace_back(iy, ix_best);
  }

  // GeodeticField3 interpolation means adjacent model-state coordinates can
  // weakly affect a gate already visited.  Repeated outward sweeps are the
  // continuous-grid counterpart of ARTS2's layer-by-layer peeling and recover
  // mutual consistency without replacing the analytical Jacobian by a table.
  auto fit_is_converged = [&]() {
    return std::ranges::all_of(gate_state, [&](const auto& pair) {
      const Size    iy       = pair.first;
      const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
      return std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]));
    });
  };
  bool all_converged = fit_is_converged();
  for (Index sweep = 1; sweep < max_sweeps and not all_converged; ++sweep) {
    for (const auto [iy, ix] : gate_state) {
      for (Index iteration = 0; iteration < max_iterations; ++iteration) {
        const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
        if (std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]))) break;
        const Numeric derivative = measurement_jac[iy, ix];
        if (not std::isfinite(derivative) or derivative == 0.0)
          arts::user_error(
              "Radar onion peeling encountered a zero/non-finite derivative for gate {} and state {}", iy, ix);
        Numeric step        = std::clamp(residual / derivative, -max_step, max_step);
        model_state_vec[ix] = std::clamp(model_state_vec[ix] + step, state_min, state_max);
        for (const auto& target : jac_targets.atm) target.update_model(atm_field, model_state_vec);
        calculate();
      }
    }
    all_converged = fit_is_converged();
  }
  if (not all_converged)
    arts::user_error("Radar onion peeling did not produce a mutually consistent fit after {} sweeps", max_sweeps);

  // Ensure the returned field, fit, and Jacobian all represent precisely the
  // returned state, ready for a subsequent OEM calculation.
  for (const auto& target : jac_targets.atm) target.update_model(atm_field, model_state_vec);
  calculate();
}
ARTS_METHOD_ERROR_CATCH
