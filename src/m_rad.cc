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

void spectral_radiance_jacobianEmpty(
    StokvecMatrix &spectral_radiance_jacobian,
    const AscendingGrid &frequency_grid,
    const JacobianTargets &jacobian_targets) try {
  ARTS_TIME_REPORT

  spectral_radiance_jacobian.resize(jacobian_targets.x_size(),
                                    frequency_grid.size());
  spectral_radiance_jacobian = Stokvec{0.0, 0.0, 0.0, 0.0};
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_jacobianFromBackground(
    StokvecMatrix &spectral_radiance_jacobian,
    const StokvecMatrix &spectral_radiance_background_jacobian,
    const MuelmatVector &background_transmittance) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      static_cast<Size>(spectral_radiance_background_jacobian.ncols()) !=
          background_transmittance.size(),
      "spectral_radiance_background_jacobian must have same number of rows as the "
      "size of jacobian_targets")

  //! The radiance derivative shape is the background shape
  spectral_radiance_jacobian.resize(
      spectral_radiance_background_jacobian.shape());

  //! Set the background radiance derivative as that which is seen after "this" swath
  for (Index i = 0; i < spectral_radiance_jacobian.nrows(); i++) {
    std::transform(background_transmittance.begin(),
                   background_transmittance.end(),
                   spectral_radiance_background_jacobian[i].begin(),
                   spectral_radiance_jacobian[i].begin(),
                   std::multiplies<>());
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_jacobianAddPathPropagation(
    StokvecMatrix &spectral_radiance_jacobian,
    const ArrayOfStokvecMatrix &ray_path_spectral_radiance_jacobian,
    const JacobianTargets &jacobian_targets,
    const AtmField &atmospheric_field,
    const ArrayOfPropagationPathPoint &ray_path) try {
  ARTS_TIME_REPORT

  const auto np = ray_path_spectral_radiance_jacobian.size();
  const auto nj = spectral_radiance_jacobian.nrows();
  const auto nf = spectral_radiance_jacobian.ncols();
  const auto nt = jacobian_targets.target_count();

  if (nt == 0) return;

  ARTS_USER_ERROR_IF(
      static_cast<Size>(spectral_radiance_jacobian.nrows()) !=
          jacobian_targets.x_size(),
      "Bad size of spectral_radiance_jacobian, it's inner dimension should match the size of jacobian_targets. Sizes: "
      "{} != {}",
      spectral_radiance_jacobian.nrows(),
      jacobian_targets.x_size())

  ARTS_USER_ERROR_IF(
      ray_path.size() != np,
      "ray_path must have same size as the size of ray_path_spectral_radiance_jacobian.  Sizes: ",
      "{} != {}",
      ray_path.size(),
      np)

  for (auto &dr : ray_path_spectral_radiance_jacobian) {
    ARTS_USER_ERROR_IF(
        dr.ncols() != nf or dr.nrows() != static_cast<Index>(nt),
        "ray_path_spectral_radiance_jacobian elements must have same number of rows as the size of "
        "jacobian_targets.  Sizes: "
        "{:B,} != [{}, {}]",
        dr.shape(),
        nt,
        nf)
  }

  //! Checks that the jacobian_targets can be used and throws if not
  jacobian_targets.throwing_check(nj);

  //! The derivative part from the atmosphere
  for (auto &atm_block : jacobian_targets.atm) {
    ARTS_USER_ERROR_IF(not atmospheric_field.contains(atm_block.type),
                       "No {} in atmospheric_field but in jacobian_targets",
                       atm_block.type)
    const auto &data = atmospheric_field[atm_block.type];
    for (Size ip = 0; ip < np; ip++) {
      const auto weights = data.flat_weight(ray_path[ip].pos);
      const auto &local  = ray_path_spectral_radiance_jacobian[ip];

      for (auto &w : weights) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.x_start;
          assert(i < static_cast<Size>(nj));
          std::transform(
              local[atm_block.target_pos].begin(),
              local[atm_block.target_pos].end(),
              spectral_radiance_jacobian[i].begin(),
              spectral_radiance_jacobian[i].begin(),
              [x = w.second](const Stokvec &a, const Stokvec &b) -> Stokvec {
                return fma(x, a, b);
              });
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_transform_operatorSet(
    SpectralRadianceTransformOperator &spectral_radiance_transform_operator,
    const SpectralRadianceUnitType &x) try {
  ARTS_TIME_REPORT

  spectral_radiance_transform_operator = SpectralRadianceTransformOperator(x);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceApplyUnit(StokvecVector &spectral_radiance,
                                StokvecMatrix &spectral_radiance_jacobian,
                                const AscendingGrid &frequency_grid,
                                const PropagationPathPoint &ray_path_point,
                                const SpectralRadianceTransformOperator
                                    &spectral_radiance_transform_operator) try {
  ARTS_TIME_REPORT

  if (spectral_radiance_jacobian.empty()) {
    spectral_radiance_jacobian.resize(0, frequency_grid.size());
  }

  spectral_radiance_transform_operator(spectral_radiance,
                                       spectral_radiance_jacobian,
                                       frequency_grid,
                                       ray_path_point);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceApplyForwardUnit(
    StokvecVector &spectral_radiance,
    const AscendingGrid &frequency_grid,
    const PropagationPathPoint &ray_path_point,
    const SpectralRadianceTransformOperator
        &spectral_radiance_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecMatrix spectral_radiance_jacobian(0, frequency_grid.size());

  spectral_radiance_transform_operator(spectral_radiance,
                                       spectral_radiance_jacobian,
                                       frequency_grid,
                                       ray_path_point);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_jacobianAddSensorJacobianPerturbations(
    const Workspace &ws,
    StokvecMatrix &spectral_radiance_jacobian,
    const StokvecVector &spectral_radiance,
    const ArrayOfSensorObsel &measurement_sensor,
    const AscendingGrid &frequency_grid,
    const JacobianTargets &jacobian_targets,
    const Vector3 &pos,
    const Vector2 &los,
    const AtmField &atmospheric_field,
    const SurfaceField &surface_field,
    const SubsurfaceField &subsurface_field,
    const Agenda &spectral_radiance_observer_agenda) try {
  ARTS_TIME_REPORT

  /*
  
  This method likely calls itself "recursively" for sensor parameters

  However, the flag for bailing and stopping this recursion
  is that there are no more jacobian targets.  So it calls
  itself always with an empty jacobian_targets.  This is how
  it bails out of the recursion.

  As this method is useless unless there are sensor elements,
  we use empty Sensor targets rather than empty jacobian_targets
  as a mini-optimization.

  */
  if (jacobian_targets.sensor.empty()) return;

  ARTS_USER_ERROR_IF(
      spectral_radiance.size() != frequency_grid.size(),
      R"(spectral_radiance must have same size as element frequency grid

spectral_radiance.size() = {},
frequency_grid.size()    = {}
)",
      spectral_radiance.size(),
      frequency_grid.size())

  ARTS_USER_ERROR_IF(
      not same_shape<2>({jacobian_targets.x_size(), frequency_grid.size()},
                        spectral_radiance_jacobian),
      R"(spectral_radiance_jacobian must be x-grid times frequency grid

spectral_radiance_jacobian.shape() = {:B,},
jacobian_targets.x_size()          = {},
frequency_grid.size()              = {}
)",
      spectral_radiance_jacobian.shape(),
      jacobian_targets.x_size(),
      frequency_grid.size())

  const JacobianTargets jacobian_targets_empty{};
  StokvecMatrix spectral_radiance_jacobian_empty{};
  ArrayOfPropagationPathPoint ray_path{};

  StokvecVector dsrad;
  auto call = [&](const AscendingGrid &frequency_grid_2,
                  const Vector3 &pos2,
                  const Vector2 &los2,
                  const Numeric d) {
    spectral_radiance_observer_agendaExecute(ws,
                                             dsrad,
                                             spectral_radiance_jacobian_empty,
                                             ray_path,
                                             frequency_grid_2,
                                             jacobian_targets_empty,
                                             pos2,
                                             los2,
                                             atmospheric_field,
                                             surface_field,
                                             subsurface_field,
                                             spectral_radiance_observer_agenda);

    // Convert to perturbed Jacobian
    dsrad -= spectral_radiance;
    dsrad /= d;
  };

  const auto &x = frequency_grid;
  const auto *b = x.begin();
  const auto *e = x.end();

  bool find_any = false;
  for (auto &target : jacobian_targets.sensor) {
    ARTS_USER_ERROR_IF(measurement_sensor.size() <=
                           static_cast<Size>(target.type.measurement_elem),
                       "Sensor element out of bounds");

    auto &elem      = measurement_sensor[target.type.measurement_elem];
    auto m          = spectral_radiance_jacobian[target.target_pos];
    const Numeric d = target.d;

    // Check that the Jacobian targets are represented by this frequency grid and this pos-los pair
    const Index iposlos = elem.find(pos, los);
    if (iposlos == SensorObsel::dont_have) continue;
    if (elem.find(frequency_grid) == SensorObsel::dont_have) continue;

    find_any = true;

    using enum SensorKeyType;
    switch (target.type.type) {
      case f:   call({b, e, [d](auto x) { return x + d; }}, pos, los, d); break;
      case za:  call(x, pos, {los[0] + d, los[1]}, d); break;
      case aa:  call(x, pos, {los[0], los[1] + d}, d); break;
      case alt: call(x, {pos[0] + d, pos[1], pos[2]}, los, d); break;
      case lat: call(x, {pos[0], pos[1] + d, pos[2]}, los, d); break;
      case lon: call(x, {pos[0], pos[1], pos[2] + d}, los, d); break;
    }

    m += dsrad;
  }

  ARTS_USER_ERROR_IF(not find_any,
                     R"(No sensor element found for pos-los/frequency grid pair

  frequency_grid: {:Bs,}
  pos:            {:B,}
  los:            {:B,}

Note: It is not allowed to change the frequency grid or the pos-los pair in an agenda
that calls this function.  This is because the actual memory address is used to identify
a sensor element.  Modifying pos, los or frequency_grid will copy the data to a new memory
location and the sensor element will not be found.
)",
                     frequency_grid,
                     pos,
                     los);
}
ARTS_METHOD_ERROR_CATCH

void measurement_vectorFromSensor(
    const Workspace &ws,
    Vector &measurement_vector,
    Matrix &measurement_jacobian,
    const ArrayOfSensorObsel &measurement_sensor,
    const JacobianTargets &jacobian_targets,
    const AtmField &atmospheric_field,
    const SurfaceField &surface_field,
    const SubsurfaceField &subsurface_field,
    const SpectralRadianceTransformOperator
        &spectral_radiance_transform_operator,
    const Agenda &spectral_radiance_observer_agenda) try {
  ARTS_TIME_REPORT

  measurement_vector.resize(measurement_sensor.size());
  measurement_vector = 0.0;

  measurement_jacobian.resize(measurement_sensor.size(),
                              jacobian_targets.x_size());
  measurement_jacobian = 0.0;

  if (measurement_sensor.empty()) return;

  //! Check the observational elements that their dimensions are correct
  for (auto &obsel : measurement_sensor) obsel.check();

  const SensorSimulations simulations = collect_simulations(measurement_sensor);

  const auto flat_size = [](const SensorSimulations &simulations) {
    Size size = 0;
    for (const auto &[f_grid_ptr, poslos_set] : simulations) {
      for (const auto &poslos_gs : poslos_set) {
        size += poslos_gs->size();
      }
    }

    return size;
  };

  const Size N = flat_size(simulations);

  struct unflatten_data {
    const std::shared_ptr<const AscendingGrid> *f_grid_ptr{nullptr};
    const std::shared_ptr<const SensorPosLosVector> *poslos_ptr{nullptr};
    Size ip{std::numeric_limits<Size>::max()};

    unflatten_data(const SensorSimulations &simulations, Size n) {
      for (auto &[f_grid_ptr, poslos_set] : simulations) {
        for (auto &poslos_gs : poslos_set) {
          const Size np = poslos_gs->size();
          if (n < np) {
            this->f_grid_ptr = &f_grid_ptr;
            this->poslos_ptr = &poslos_gs;
            this->ip         = n;
            goto end;
          }
          n -= np;
        }
      }

    end:
      ARTS_USER_ERROR_IF(
          f_grid_ptr == nullptr or poslos_ptr == nullptr,
          "Failed to find f_grid_ptr and poslos_ptr for index {} in simulations",
          n)
      ARTS_USER_ERROR_IF(ip >= (**poslos_ptr).size(),
                         "Index {} out of bounds for poslos_ptr with size {}",
                         ip,
                         (**poslos_ptr).size());
    }
  };

  std::string error{};

#pragma omp parallel for if (not arts_omp_in_parallel() and N > 1)
  for (Size i = 0; i < N; i++) {
    try {
      const unflatten_data unflat(simulations, i);

      const Size ip    = unflat.ip;
      auto &f_grid_ptr = *unflat.f_grid_ptr;
      auto &poslos     = (**unflat.poslos_ptr)[ip];

      StokvecVector spectral_radiance;
      StokvecMatrix spectral_radiance_jacobian;
      ArrayOfPropagationPathPoint ray_path;

      spectral_radiance_observer_agendaExecute(
          ws,
          spectral_radiance,
          spectral_radiance_jacobian,
          ray_path,
          *f_grid_ptr,
          jacobian_targets,
          poslos.pos,
          poslos.los,
          atmospheric_field,
          surface_field,
          subsurface_field,
          spectral_radiance_observer_agenda);

      ARTS_USER_ERROR_IF(ray_path.empty(), "No ray path found");
      spectral_radiance_transform_operator(spectral_radiance,
                                           spectral_radiance_jacobian,
                                           *f_grid_ptr,
                                           ray_path.front());

#pragma omp critical
      for (Size iv = 0; iv < measurement_sensor.size(); ++iv) {
        const SensorObsel &obsel = measurement_sensor[iv];
        if (obsel.same_freqs(f_grid_ptr)) {
          measurement_vector[iv] += obsel.sumup(spectral_radiance, ip);

          obsel.sumup(measurement_jacobian[iv], spectral_radiance_jacobian, ip);
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
