#include <enumsAtmKey.h>
#include <subsurf_field.h>
#include <workspace.h>

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void RetrievalInit(JacobianTargets& jac_targets,
                   CovarianceMatrix& model_state_covariance_matrix,
                   JacobianTargetsDiagonalCovarianceMatrixMap&
                       covariance_matrix_diagonal_blocks) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixInit(model_state_covariance_matrix);
  jac_targetsInit(jac_targets);
  covariance_matrix_diagonal_blocks.clear();
}

void RetrievalAddSubsurface(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SubsurfaceKey& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSubsurface(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.subsurf.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddSubsurface(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SubsurfacePropertyTag& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSubsurface(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.subsurf.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddSurface(JacobianTargets& jac_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfaceKey& key,
                         const Numeric& d,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSurface(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.surf.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddSurface(JacobianTargets& jac_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfacePropertyTag& key,
                         const Numeric& d,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSurface(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.surf.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const AtmKey& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddAtmosphere(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddAtmosphere(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesIsotope& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddAtmosphere(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const QuantumLevelIdentifier& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddAtmosphere(jac_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddSpeciesVMR(JacobianTargets& jac_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& species,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSpeciesVMR(jac_targets, species, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddSpeciesIsotopologueRatio(
    JacobianTargets& jac_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const SpeciesIsotope& species,
    const Numeric& d,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSpeciesIsotopologueRatio(jac_targets, species, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddMagneticField(JacobianTargets& jac_targets,
                               JacobianTargetsDiagonalCovarianceMatrixMap&
                                   covariance_matrix_diagonal_blocks,
                               const String& component,
                               const Numeric& d,
                               const BlockMatrix& matrix,
                               const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddMagneticField(jac_targets, component, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddOverlappingMagneticField(
    JacobianTargets& jac_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddOverlappingMagneticField(jac_targets);

  const Size N = jac_targets.atm.size();

  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm[N - 2].type}] = {.first = matrix, .second = inverse};

  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm[N - 1].type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddOverlappingWindField(
    JacobianTargets& jac_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddOverlappingWindField(jac_targets);

  const Size N = jac_targets.atm.size();

  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm[N - 2].type}] = {.first = matrix, .second = inverse};

  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm[N - 1].type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddWindField(JacobianTargets& jac_targets,
                           JacobianTargetsDiagonalCovarianceMatrixMap&
                               covariance_matrix_diagonal_blocks,
                           const String& component,
                           const Numeric& d,
                           const BlockMatrix& matrix,
                           const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddWindField(jac_targets, component, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddTemperature(JacobianTargets& jac_targets,
                             JacobianTargetsDiagonalCovarianceMatrixMap&
                                 covariance_matrix_diagonal_blocks,
                             const Numeric& d,
                             const BlockMatrix& matrix,
                             const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddTemperature(jac_targets, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddPressure(JacobianTargets& jac_targets,
                          JacobianTargetsDiagonalCovarianceMatrixMap&
                              covariance_matrix_diagonal_blocks,
                          const Numeric& d,
                          const BlockMatrix& matrix,
                          const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddPressure(jac_targets, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jac_targets.atm.back().type}] = {.first = matrix, .second = inverse};
}

void RetrievalAddSensorFrequencyPolyOffset(
    JacobianTargets& jac_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const ArrayOfSensorObsel& measurement_sensor,
    const Numeric& d,
    const Index& sensor_elem,
    const Index& polyorder,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddSensorFrequencyPolyOffset(
      jac_targets, measurement_sensor, d, sensor_elem, polyorder);
  auto keyk = JacobianTargetType{jac_targets.sensor.back().type};
  covariance_matrix_diagonal_blocks[keyk] = {.first  = matrix,
                                             .second = inverse};
}

void RetrievalAddErrorPolyFit(JacobianTargets& jac_targets,
                              JacobianTargetsDiagonalCovarianceMatrixMap&
                                  covariance_matrix_diagonal_blocks,
                              const ArrayOfSensorObsel& measurement_sensor,
                              const Vector& t,
                              const Index& sensor_elem,
                              const Index& polyorder,
                              const BlockMatrix& matrix,
                              const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  jac_targetsAddErrorPolyFit(
      jac_targets, measurement_sensor, t, sensor_elem, polyorder);
  auto keyk = JacobianTargetType{jac_targets.error.back().type};
  covariance_matrix_diagonal_blocks[keyk] = {.first  = matrix,
                                             .second = inverse};
}
