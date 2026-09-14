#include <enumsAtmKey.h>
#include <subsurf_field.h>
#include <workspace.h>

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void oemInit(JacobianTargets& jac_targets, OptimalEstimationData& oem) {
  ARTS_TIME_REPORT
  oem.clear();
  jac_targetsInit(jac_targets);
}

void oemAddSubsurface(JacobianTargets&                                                jac_targets,
                      OptimalEstimationData&                                          oem,
                      const Generic<const SubsurfaceKey, const SubsurfacePropertyTag> key,
                      const Numeric&                                                  d,
                      const BlockMatrix&                                              matrix,
                      const BlockMatrix&                                              inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddSubsurface(jac_targets, key, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.subsurf.back().type}] = {.first  = matrix,
                                                                                     .second = inverse};
}

void oemAddSurface(JacobianTargets&                                          jac_targets,
                   OptimalEstimationData&                                    oem,
                   const Generic<const SurfaceKey, const SurfacePropertyTag> key,
                   const Numeric&                                            d,
                   const BlockMatrix&                                        matrix,
                   const BlockMatrix&                                        inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddSurface(jac_targets, key, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.surf.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddAtmosphere(JacobianTargets&                    jac_targets,
                      OptimalEstimationData&              oem,
                      const Generic<const AtmKey,
                                    const QuantumLevelIdentifier,
                                    const ScatteringSpeciesProperty,
                                    const SpeciesEnum,
                                    const SpeciesIsotope> key,
                      const Numeric&                      d,
                      const BlockMatrix&                  matrix,
                      const BlockMatrix&                  inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddAtmosphere(jac_targets, key, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddSpeciesVMR(JacobianTargets&       jac_targets,
                      OptimalEstimationData& oem,
                      const SpeciesEnum&     species,
                      const Numeric&         d,
                      const BlockMatrix&     matrix,
                      const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddSpeciesVMR(jac_targets, species, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddSpeciesIsotopologueRatio(JacobianTargets&       jac_targets,
                                    OptimalEstimationData& oem,
                                    const SpeciesIsotope&  species,
                                    const Numeric&         d,
                                    const BlockMatrix&     matrix,
                                    const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddSpeciesIsotopologueRatio(jac_targets, species, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddMagneticField(JacobianTargets&       jac_targets,
                         OptimalEstimationData& oem,
                         const String&          component,
                         const Numeric&         d,
                         const BlockMatrix&     matrix,
                         const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddMagneticField(jac_targets, component, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddOverlappingMagneticField(JacobianTargets&       jac_targets,
                                    OptimalEstimationData& oem,
                                    const BlockMatrix&     matrix,
                                    const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddOverlappingMagneticField(jac_targets);

  const Size N = jac_targets.atm.size();

  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm[N - 2].type}] = {.first = matrix, .second = inverse};

  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm[N - 1].type}] = {.first = matrix, .second = inverse};
}

void oemAddOverlappingWindField(JacobianTargets&       jac_targets,
                                OptimalEstimationData& oem,
                                const BlockMatrix&     matrix,
                                const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddOverlappingWindField(jac_targets);

  const Size N = jac_targets.atm.size();

  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm[N - 2].type}] = {.first = matrix, .second = inverse};

  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm[N - 1].type}] = {.first = matrix, .second = inverse};
}

void oemAddWindField(JacobianTargets&       jac_targets,
                     OptimalEstimationData& oem,
                     const String&          component,
                     const Numeric&         d,
                     const BlockMatrix&     matrix,
                     const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddWindField(jac_targets, component, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddTemperature(JacobianTargets&       jac_targets,
                       OptimalEstimationData& oem,
                       const Numeric&         d,
                       const BlockMatrix&     matrix,
                       const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddTemperature(jac_targets, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddPressure(JacobianTargets&       jac_targets,
                    OptimalEstimationData& oem,
                    const Numeric&         d,
                    const BlockMatrix&     matrix,
                    const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddPressure(jac_targets, d);
  oem.covmat_diagonal_blocks[JacobianTargetType{jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void oemAddSensorFrequencyPolyOffset(JacobianTargets&          jac_targets,
                                     OptimalEstimationData&    oem,
                                     const ArrayOfSensorObsel& measurement_sensor,
                                     const Numeric&            d,
                                     const Index&              sensor_elem,
                                     const Index&              polyorder,
                                     const BlockMatrix&        matrix,
                                     const BlockMatrix&        inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddSensorFrequencyPolyOffset(jac_targets, measurement_sensor, d, sensor_elem, polyorder);
  auto keyk                        = JacobianTargetType{jac_targets.sensor.back().type};
  oem.covmat_diagonal_blocks[keyk] = {.first = matrix, .second = inverse};
}

void oemAddErrorPolyFit(JacobianTargets&          jac_targets,
                        OptimalEstimationData&    oem,
                        const ArrayOfSensorObsel& measurement_sensor,
                        const Vector&             t,
                        const Index&              sensor_elem,
                        const Index&              polyorder,
                        const BlockMatrix&        matrix,
                        const BlockMatrix&        inverse) {
  ARTS_TIME_REPORT
  oem.uncheck();

  jac_targetsAddErrorPolyFit(jac_targets, measurement_sensor, t, sensor_elem, polyorder);
  auto keyk                        = JacobianTargetType{jac_targets.error.back().type};
  oem.covmat_diagonal_blocks[keyk] = {.first = matrix, .second = inverse};
}
