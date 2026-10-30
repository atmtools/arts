#include <arts_constants.h>
#include <mc_antenna.h>
#include <path_point.h>
#include <radar.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

void MCRadar(const Workspace&                ws,
             StokvecVector&                  radar_signal,
             StokvecVector&                  radar_error,
             Vector&                         radar_time,
             const AtmField&                 atm_field,
             const SurfaceField&             surf_field,
             const ArrayOfScatteringSpecies& scattering_species,
             const MCAntenna&                mc_antenna,
             const Agenda&                   ray_path_observer_agenda,
             const Agenda&                   spectral_propmat_agenda,
             const TransmittanceOption&      rte_option,
             const Numeric&                  frequency,
             const Vector3&                  sensor_pos,
             const Vector2&                  sensor_los,
             const Stokvec&                  mc_y_tx,
             const AscendingGrid&            range_bins,
             const Index&                    mc_seed,
             const Index&                    mc_max_iter,
             const Index&                    mc_max_scatorder,
             const Numeric&                  k2,
             const RadarUnit&                unit,
             const RadarRangeMode&           range_mode,
             const Numeric&                  ze_tref,
             const Numeric&                  pext_scaling) try {
  ARTS_TIME_REPORT

  const RadarRangeMode mode = radar::resolve(range_mode, range_bins.empty() ? 0.0 : range_bins.back());

  if (frequency <= 0.0) arts::user_error("frequency must be positive");
  if (range_bins.size() < 2) arts::user_error("range_bins must contain at least two edges");
  if (mode != RadarRangeMode::Altitude and range_bins.front() < 0.0)
    arts::user_error("range_bins cannot contain negative distances or times");
  if (mc_max_iter <= 0) arts::user_error("mc_max_iter must be positive");
  if (mc_max_scatorder != 1)
    arts::user_error("This first ARTS3 MCRadar implementation supports mc_max_scatorder=1; got {}", mc_max_scatorder);
  if (scattering_species.species.empty()) arts::user_error("MCRadar requires at least one scattering species");
  if (unit == RadarUnit::dBZe)
    arts::user_error("MCRadar returns a mean and its standard error, which are linear quantities; use Ze, not dBZe");
  if (pext_scaling < 0.0 or pext_scaling > 2.0)
    arts::user_error("pext_scaling must be between 0 and 2; got {}", pext_scaling);

  const Size            nbins = range_bins.size() - 1;
  const AscendingGrid   freq_grid{frequency};
  const JacobianTargets jac_targets{};
  const Numeric         factor = radar::unit_factor(unit, frequency, ze_tref, k2);

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
                                   freq_grid,
                                   Vector3{0.0, 0.0, 0.0},
                                   jac_targets,
                                   {},
                                   point,
                                   atm_point,
                                   spectral_propmat_agenda);
  };

  // The mean is weighted by the antenna gain of every sampled direction, so a
  // uniformly filled beam returns the same value as a pencil beam.  Weighted
  // sums of the return, and of its square, give the standard error.
  Numeric           weight_sum         = 0.0;
  Numeric           squared_weight_sum = 0.0;
  StokvecVector     sum(nbins, Stokvec{});
  StokvecVector     squared_weight_signal_sum(nbins, Stokvec{});
  StokvecVector     squared_sum(nbins, Stokvec{});
  Vector            time_sum(nbins, 0.0);
  Vector            coverage_sum(nbins, 0.0);
  std::vector<bool> covered(nbins, false);

  RandomNumberGenerator<> rng(mc_seed);
  // The general RNG avoids reusing seeds process-wide.  A forward model with
  // an explicit seed instead promises reproducibility across repeated calls.
  rng.force_seed(mc_seed);
  const Matrix33 ant_to_enu = rotmat_enu(sensor_los);
  Matrix33       enu_to_ant;
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j) enu_to_ant[i, j] = ant_to_enu[j, i];

  StokvecMatrix  signal;
  StokvecTensor4 signal_jac;
  for (Index iter = 0; iter < mc_max_iter; ++iter) {
    const auto [sampled_los, ray_rotation]  = mc_antenna.draw_los(rng, ant_to_enu, sensor_los);
    const Numeric antenna_weight            = mc_antenna.return_los(ray_rotation, enu_to_ant);
    const Muelmat tx_rotation               = rotmat_stokes(1.0, 1.0, ant_to_enu, ray_rotation);
    const Muelmat rx_rotation               = rotmat_stokes(-1.0, 1.0, ray_rotation, ant_to_enu);
    weight_sum                             += antenna_weight;
    squared_weight_sum                     += antenna_weight * antenna_weight;

    ArrayOfPropagationPathPoint ray_path;
    ray_path_observer_agendaExecute(ws, ray_path, sensor_pos, sampled_los, ray_path_observer_agenda);

    const radar::PathRange             range      = radar::path_range(ray_path, sensor_pos, surf_field.ellipsoid);
    const Vector                       coordinate = radar::range_coordinate(range, ray_path, mode);
    std::vector<rtepack::radar_sample> samples;
    for (Size ib = 0; ib < nbins; ++ib)
      radar::add_gate_samples(samples, coordinate, range_bins[ib], range_bins[ib + 1], ib);
    const Size np = radar::sort_samples(samples);
    if (samples.empty()) continue;

    const radar::Optics optics =
        radar::optics(ray_path, np, atm_field, freq_grid, scattering_species, jac_targets, pext_scaling, gas);
    radar::single_scattering(
        signal, signal_jac, tx_rotation * mc_y_tx, optics, range.layer, samples, nbins, rte_option);

    const Vector coverage = radar::gate_coverage(samples, nbins);
    const Vector time     = radar::gate_sum(samples, range.time, nbins);
    for (Size ib = 0; ib < nbins; ++ib) {
      if (coverage[ib] == 0.0) continue;

      covered[ib]       = true;
      time_sum[ib]     += antenna_weight * time[ib];
      coverage_sum[ib] += antenna_weight * coverage[ib];

      const Stokvec value = rx_rotation * signal[ib, 0];
      for (Size s = 0; s < 4; ++s) {
        sum[ib][s]                       += antenna_weight * value[s];
        squared_weight_signal_sum[ib][s] += antenna_weight * antenna_weight * value[s];
        squared_sum[ib][s]               += antenna_weight * antenna_weight * value[s] * value[s];
      }
    }
  }

  if (not(weight_sum > 0.0)) arts::user_error("No sampled antenna direction has a positive antenna weight");

  constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();
  radar_signal.resize(nbins);
  radar_error.resize(nbins);
  radar_time.resize(nbins);
  for (Size ib = 0; ib < nbins; ++ib) {
    if (not covered[ib]) {
      radar_signal[ib] = Stokvec{nan, nan, nan, nan};
      radar_error[ib]  = Stokvec{nan, nan, nan, nan};
      radar_time[ib]   = nan;
      continue;
    }

    radar_time[ib] = coverage_sum[ib] > 0.0 ? time_sum[ib] / coverage_sum[ib] : nan;
    for (Size s = 0; s < 4; ++s) {
      const Numeric mean     = sum[ib][s] / weight_sum;
      const Numeric variance = std::max(0.0,
                                        squared_sum[ib][s] - 2.0 * mean * squared_weight_signal_sum[ib][s] +
                                            mean * mean * squared_weight_sum) /
                               (weight_sum * weight_sum);
      radar_signal[ib][s] = factor * mean;
      radar_error[ib][s]  = factor * std::sqrt(variance);
    }
  }
}
ARTS_METHOD_ERROR_CATCH
