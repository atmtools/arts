#include <arts_omp.h>
#include <atm.h>
#include <jacobian.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>
#include <exception>

void spectral_rad_jacEmpty(StokvecMatrix &spectral_rad_jac,
                           const AscendingGrid &freq_grid,
                           const JacobianTargets &jac_targets) try {
  ARTS_TIME_REPORT

  spectral_rad_jac.resize(jac_targets.x_size(), freq_grid.size());
  spectral_rad_jac = Stokvec{0.0, 0.0, 0.0, 0.0};
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_jacFromBackground(
    StokvecMatrix &spectral_rad_jac,
    const StokvecMatrix &spectral_rad_bkg_jac,
    const TransmittanceMatrix &spectral_tramat) try {
  ARTS_TIME_REPORT

  const auto [nf, np, nq] = spectral_tramat.shape();

  //! The radiance derivative shape is the background shape
  spectral_rad_jac.resize(spectral_rad_bkg_jac.shape());

  if (nq == 0 or np == 0) return;

  const auto &&background_transmittance = spectral_tramat.P[joker, np - 1];

  ARTS_USER_ERROR_IF(
      background_transmittance.ncols() != spectral_rad_bkg_jac.ncols(),
      "Bad size of spectral_rad_bkg_jac, its inner dimension should match the frequency size of spectral_tramat. Sizes: {} != {}",
      spectral_rad_bkg_jac.ncols(),
      background_transmittance.ncols());

  //! Set the background radiance derivative as that which is seen after "this" swath
  for (Index i = 0; i < spectral_rad_jac.nrows(); i++) {
    const auto b = spectral_rad_bkg_jac[i];
    auto s       = spectral_rad_jac[i];
    std::transform(background_transmittance.begin(),
                   background_transmittance.end(),
                   b.begin(),
                   s.begin(),
                   std::multiplies<>());
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_jacAddPathPropagation(
    StokvecMatrix &spectral_rad_jac,
    const StokvecTensor3 &spectral_rad_jac_path,
    const JacobianTargets &jac_targets,
    const AtmField &atm_field,
    const ArrayOfPropagationPathPoint &ray_path) try {
  ARTS_TIME_REPORT

  const Size nf = spectral_rad_jac_path.npages();
  const Size np = spectral_rad_jac_path.nrows();
  const Size nt = jac_targets.target_count();
  const Size nx = jac_targets.x_size();

  jac_targets.throwing_check(nx);

  ARTS_USER_ERROR_IF(not same_shape({nx, nf}, spectral_rad_jac) or
                         ray_path.size() != np or
                         not same_shape({nf, np, nt}, spectral_rad_jac_path),
                     R"(Mismatched input sizes:

nf : {}
np : {}
nt : {}
nx : {}

spectral_rad_jac.shape()      : {:B,} [nx, nf]
spectral_rad_jac_path.shape() : {:B,} [nf, np, nt]
ray_path.size()               : [{}] [np]
)",
                     nf,
                     np,
                     nt,
                     nx,
                     spectral_rad_jac.shape(),
                     spectral_rad_jac_path.shape(),
                     ray_path.size())

  if (nt == 0) return;

  //! The derivative part from the atmosphere
  for (auto &atm_block : jac_targets.atm) {
    ARTS_USER_ERROR_IF(not atm_field.contains(atm_block.type),
                       "No {} in atm_field but in jac_targets",
                       atm_block.type)
    const auto &data = atm_field[atm_block.type];
    for (Size ip = 0; ip < np; ip++) {
      const auto weights = data.flat_weight(ray_path[ip].pos);
      const auto local = spectral_rad_jac_path[joker, ip, atm_block.target_pos];

      for (auto &w : weights) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.x_start;
          assert(i < static_cast<Size>(nx));
          auto sr = spectral_rad_jac[i];
          std::transform(
              local.begin(),
              local.end(),
              sr.begin(),
              sr.begin(),
              [x = w.second](const Stokvec &a, const Stokvec &b) -> Stokvec {
                return fma(x, a, b);
              });
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_transform_operatorSet(
    SpectralRadianceTransformOperator &spectral_rad_transform_operator,
    const SpectralRadianceUnitType &x) try {
  ARTS_TIME_REPORT

  spectral_rad_transform_operator = SpectralRadianceTransformOperator(x);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radApplyUnit(StokvecVector &spectral_rad,
                           StokvecMatrix &spectral_rad_jac,
                           const AscendingGrid &freq_grid,
                           const PropagationPathPoint &ray_point,
                           const SpectralRadianceTransformOperator
                               &spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  if (spectral_rad_jac.empty()) spectral_rad_jac.resize(0, freq_grid.size());

  spectral_rad_transform_operator(
      spectral_rad, spectral_rad_jac, freq_grid, ray_point);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radApplyForwardUnit(StokvecVector &spectral_rad,
                                  const AscendingGrid &freq_grid,
                                  const PropagationPathPoint &ray_point,
                                  const SpectralRadianceTransformOperator
                                      &spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecMatrix spectral_rad_jac(0, freq_grid.size());

  spectral_rad_transform_operator(
      spectral_rad, spectral_rad_jac, freq_grid, ray_point);
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_jacAddSensorJacobianPerturbations(
    const Workspace &ws,
    StokvecMatrix &spectral_rad_jac,
    const StokvecVector &spectral_rad,
    const ArrayOfSensorObsel &measurement_sensor,
    const AscendingGrid &freq_grid,
    const JacobianTargets &jac_targets,
    const Vector3 &pos,
    const Vector2 &los,
    const AtmField &atm_field,
    const SurfaceField &surf_field,
    const SubsurfaceField &subsurf_field,
    const Agenda &spectral_rad_observer_agenda) try {
  ARTS_TIME_REPORT

  /*
  
  This method likely calls itself "recursively" for sensor parameters

  However, the flag for bailing and stopping this recursion
  is that there are no more jacobian targets.  So it calls
  itself always with an empty jac_targets.  This is how
  it bails out of the recursion.

  As this method is useless unless there are sensor elements,
  we use empty Sensor targets rather than empty jac_targets
  as a mini-optimization.

  */
  if (jac_targets.sensor.empty()) return;

  ARTS_USER_ERROR_IF(
      spectral_rad.size() != freq_grid.size(),
      R"(spectral_rad must have same size as element frequency grid

spectral_rad.size() = {},
freq_grid.size()    = {}
)",
      spectral_rad.size(),
      freq_grid.size())

  ARTS_USER_ERROR_IF(not same_shape({jac_targets.x_size(), freq_grid.size()},
                                    spectral_rad_jac),
                     R"(spectral_rad_jac must be x-grid times frequency grid

spectral_rad_jac.shape() = {:B,},
jac_targets.x_size()     = {},
freq_grid.size()         = {}
)",
                     spectral_rad_jac.shape(),
                     jac_targets.x_size(),
                     freq_grid.size())

  const JacobianTargets jac_targets_empty{};
  StokvecMatrix spectral_rad_jac_empty{};
  ArrayOfPropagationPathPoint ray_path{};

  StokvecVector dsrad;
  auto call = [&](const AscendingGrid &freq_grid_2,
                  const Vector3 &pos2,
                  const Vector2 &los2,
                  const Numeric d) {
    spectral_rad_observer_agendaExecute(ws,
                                        dsrad,
                                        spectral_rad_jac_empty,
                                        ray_path,
                                        freq_grid_2,
                                        jac_targets_empty,
                                        pos2,
                                        los2,
                                        atm_field,
                                        surf_field,
                                        subsurf_field,
                                        spectral_rad_observer_agenda);

    // Convert to perturbed Jacobian
    dsrad -= spectral_rad;
    dsrad /= d;
  };

  const auto &x = freq_grid;
  const auto *b = x.begin();
  const auto *e = x.end();

  bool find_any = false;
  for (auto &target : jac_targets.sensor) {
    ARTS_USER_ERROR_IF(measurement_sensor.size() <=
                           static_cast<Size>(target.type.measurement_elem),
                       "Sensor element out of bounds");

    auto &elem      = measurement_sensor[target.type.measurement_elem];
    auto m          = spectral_rad_jac[target.target_pos];
    const Numeric d = target.d;

    // Check that the Jacobian targets are represented by this frequency grid and this pos-los pair
    const Index iposlos = elem.find(pos, los);
    if (iposlos == SensorObsel::dont_have) continue;
    if (elem.find(freq_grid) == SensorObsel::dont_have) continue;

    find_any = true;

    using enum SensorKeyType;
    switch (target.type.type) {
      case freq:
        call({b, e, [d](auto x) { return x + d; }}, pos, los, d);
        break;
      case zen: call(x, pos, {los[0] + d, los[1]}, d); break;
      case azi: call(x, pos, {los[0], los[1] + d}, d); break;
      case alt: call(x, {pos[0] + d, pos[1], pos[2]}, los, d); break;
      case lat: call(x, {pos[0], pos[1] + d, pos[2]}, los, d); break;
      case lon: call(x, {pos[0], pos[1], pos[2] + d}, los, d); break;
    }

    m += dsrad;
  }

  ARTS_USER_ERROR_IF(not find_any,
                     R"(No sensor element found for pos-los/frequency grid pair

  freq_grid: {:Bs,}
  pos:       {:B,}
  los:       {:B,}

Note: It is not allowed to change the frequency grid or the pos-los pair in an agenda
that calls this function.  This is because the actual memory address is used to identify
a sensor element.  Modifying pos, los or freq_grid will copy the data to a new memory
location and the sensor element will not be found.
)",
                     freq_grid,
                     pos,
                     los);
}
ARTS_METHOD_ERROR_CATCH

void measurement_vecFromSensor(
    const Workspace &ws,
    Vector &measurement_vec,
    Matrix &measurement_jac,
    const ArrayOfSensorObsel &measurement_sensor,
    const JacobianTargets &jac_targets,
    const AtmField &atm_field,
    const SurfaceField &surf_field,
    const SubsurfaceField &subsurf_field,
    const SpectralRadianceTransformOperator &spectral_rad_transform_operator,
    const Agenda &spectral_rad_observer_agenda) try {
  ARTS_TIME_REPORT

  measurement_vec.resize(measurement_sensor.size());
  measurement_vec = 0.0;

  measurement_jac.resize(measurement_sensor.size(), jac_targets.x_size());
  measurement_jac = 0.0;

  if (measurement_sensor.empty()) return;

  //! Check the observational elements that their dimensions are correct
  for (auto &obsel : measurement_sensor) obsel.check();

  const SensorSimulations simulations = collect_simulations(measurement_sensor);
  const Size N                        = simulations.size();

  std::string error{};

  //! Dynamic scheduling as some simulations may take much longer time than others
#pragma omp parallel for schedule(dynamic) if (arts_omp_parallel(-1, N > 1))
  for (Size i = 0; i < N; i++) {
    try {
      const Size ip         = simulations[i].iposlos;
      const auto &freq_grid = simulations[i].freq_grid;
      const auto &poslos    = simulations[i].poslos_grid[ip];

      StokvecVector spectral_rad;
      StokvecMatrix spectral_rad_jac;
      ArrayOfPropagationPathPoint ray_path;

      spectral_rad_observer_agendaExecute(ws,
                                          spectral_rad,
                                          spectral_rad_jac,
                                          ray_path,
                                          freq_grid,
                                          jac_targets,
                                          poslos.pos,
                                          poslos.los,
                                          atm_field,
                                          surf_field,
                                          subsurf_field,
                                          spectral_rad_observer_agenda);

      ARTS_USER_ERROR_IF(ray_path.empty(), "No ray path found");
      spectral_rad_transform_operator(
          spectral_rad, spectral_rad_jac, freq_grid, ray_path.front());

#pragma omp critical
      for (Size iv = 0; iv < measurement_sensor.size(); ++iv) {
        const SensorObsel &obsel = measurement_sensor[iv];
        if (obsel.same_freqs(freq_grid)) {
          measurement_vec[iv] += obsel.sumup(spectral_rad, ip);

          obsel.sumup(measurement_jac[iv], spectral_rad_jac, ip);
        }
      }
    } catch (const std::exception &e) {
#pragma omp critical
      if (error.empty()) {
        error = std::format(
            "Error in unflattening data for index {}: {}\n", i, e.what());
      }
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "Errors occurred:\n{:}", error);
}
ARTS_METHOD_ERROR_CATCH
