#include <workspace.h>

void model_state_vecInit(Vector& model_state_vec,
                         const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  model_state_vec.resize(jac_targets.x_size());
  model_state_vec = 0.0;
}

void model_state_vecPerturbations(Vector& model_state_vec,
                                  const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  model_state_vecInit(model_state_vec, jac_targets);

  for (auto& target : jac_targets.atm) {
    const Range r(target.x_start, target.x_size);
    model_state_vec[r] = target.d;
  }

  for (auto& target : jac_targets.surf) {
    const Range r(target.x_start, target.x_size);
    model_state_vec[r] = target.d;
  }

  for (auto& target : jac_targets.subsurf) {
    const Range r(target.x_start, target.x_size);
    model_state_vec[r] = target.d;
  }

  for (auto& target : jac_targets.line) {
    const Range r(target.x_start, target.x_size);
    model_state_vec[r] = target.d;
  }

  for (auto& target : jac_targets.sensor) {
    const Range r(target.x_start, target.x_size);
    model_state_vec[r] = target.d;
  }

  for (auto& target : jac_targets.error) {
    const Range r(target.x_start, target.x_size);
    model_state_vec[r] = target.d;
  }
}

////// Update the fields from the model state vector

void atm_fieldFromModelState(AtmField& atm_field,
                             const Vector& model_state_vec,
                             const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.atm) {
    target.update_model(atm_field, model_state_vec);
  }
}

void surf_fieldFromModelState(SurfaceField& surf_field,
                              const Vector& model_state_vec,
                              const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.surf) {
    target.update_model(surf_field, model_state_vec);
  }
}

void subsurf_fieldFromModelState(SubsurfaceField& subsurf_field,
                                 const Vector& model_state_vec,
                                 const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.subsurf) {
    target.update_model(subsurf_field, model_state_vec);
  }
}

void abs_bandsFromModelState(AbsorptionBands& abs_bands,
                             const Vector& model_state_vec,
                             const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.line) {
    target.update_model(abs_bands, model_state_vec);
  }
}

void measurement_sensorFromModelState(ArrayOfSensorObsel& measurement_sensor,
                                      const Vector& model_state_vec,
                                      const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.sensor) {
    target.update_model(measurement_sensor, model_state_vec);
  }
}

////// Update the model state vector from the fields

void model_state_vecFromAtmosphere(Vector& model_state_vec,
                                   const AtmField& atm_field,
                                   const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.atm) {
    target.update_state(model_state_vec, atm_field);
  }
}

void model_state_vecFromSurface(Vector& model_state_vec,
                                const SurfaceField& surf_field,
                                const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.surf) {
    target.update_state(model_state_vec, surf_field);
  }
}

void model_state_vecFromSubsurface(Vector& model_state_vec,
                                   const SubsurfaceField& subsurf_field,
                                   const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.subsurf) {
    target.update_state(model_state_vec, subsurf_field);
  }
}

void model_state_vecFromBands(Vector& model_state_vec,
                              const AbsorptionBands& abs_bands,
                              const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.line) {
    target.update_state(model_state_vec, abs_bands);
  }
}

void model_state_vecFromSensor(Vector& model_state_vec,
                               const ArrayOfSensorObsel& measurement_sensor,
                               const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.sensor) {
    target.update_state(model_state_vec, measurement_sensor);
  }
}

///// Update the measurement Jacobian from the model state vector

void measurement_jacAtmosphereTransformation(
    Matrix& measurement_jac,
    const Vector& model_state_vec,
    const AtmField& field,
    const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.atm) {
    target.update_jac(measurement_jac, model_state_vec, field);
  }
}

void measurement_jacSurfaceTransformation(Matrix& measurement_jac,
                                          const Vector& model_state_vec,
                                          const SurfaceField& field,
                                          const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.surf) {
    target.update_jac(measurement_jac, model_state_vec, field);
  }
}

void measurement_jacSubsurfaceTransformation(
    Matrix& measurement_jac,
    const Vector& model_state_vec,
    const SubsurfaceField& field,
    const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.subsurf) {
    target.update_jac(measurement_jac, model_state_vec, field);
  }
}

void measurement_jacBandTransformation(Matrix& measurement_jac,
                                       const Vector& model_state_vec,
                                       const AbsorptionBands& field,
                                       const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.line) {
    target.update_jac(measurement_jac, model_state_vec, field);
  }
}

void measurement_jacSensorTransformation(Matrix& measurement_jac,
                                         const Vector& model_state_vec,
                                         const ArrayOfSensorObsel& field,
                                         const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jac_targets.sensor) {
    target.update_jac(measurement_jac, model_state_vec, field);
  }
}
