#include <workspace.h>

#include "enumsAtmKey.h"
#include "enumsFieldComponent.h"

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void RetrievalInit(JacobianTargets& jacobian_targets,
                   CovarianceMatrix& model_state_covariance_matrix,
                   JacobianTargetsDiagonalCovarianceMatrixMap&
                       covariance_matrix_diagonal_blocks) {
  model_state_covariance_matrixInit(model_state_covariance_matrix);
  jacobian_targetsInit(jacobian_targets);
  covariance_matrix_diagonal_blocks.clear();
}

void RetrievalAddSurface(JacobianTargets& jacobian_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfaceKey& key,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse,
                         const Numeric& d) {
  jacobian_targetsAddSurface(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = SurfaceKeyVal{key}}, matrix, inverse);
}

void RetrievalAddSurface(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SurfaceTypeTag& key,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& d) {
  jacobian_targetsAddSurface(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = SurfaceKeyVal{key}}, matrix, inverse);
}

void RetrievalAddSurface(JacobianTargets& jacobian_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfacePropertyTag& key,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse,
                         const Numeric& d) {
  jacobian_targetsAddSurface(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = SurfaceKeyVal{key}}, matrix, inverse);
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const AtmKey& key,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& d) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{key}}, matrix, inverse);
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& key,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& d) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{key}}, matrix, inverse);
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesIsotope& key,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& d) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{key}}, matrix, inverse);
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const QuantumIdentifier& key,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& d) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{key}}, matrix, inverse);
}


void RetrievalAddSpeciesVMR(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& species,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& d) {
  jacobian_targetsAddSpeciesVMR(jacobian_targets, species, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{species}}, matrix, inverse);
}

void RetrievalAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const SpeciesIsotope& species,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse,
    const Numeric& d) {
  jacobian_targetsAddSpeciesIsotopologueRatio(jacobian_targets, species, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{species}}, matrix, inverse);
}

void RetrievalAddMagneticField(JacobianTargets& jacobian_targets,
                               JacobianTargetsDiagonalCovarianceMatrixMap&
                                   covariance_matrix_diagonal_blocks,
                               const String& component,
                               const BlockMatrix& matrix,
                               const BlockMatrix& inverse,
                               const Numeric& d) {
  jacobian_targetsAddMagneticField(jacobian_targets, component, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = to_mag(component)}, matrix, inverse);
}

void RetrievalAddWindField(JacobianTargets& jacobian_targets,
                           JacobianTargetsDiagonalCovarianceMatrixMap&
                               covariance_matrix_diagonal_blocks,
                           const String& component,
                           const BlockMatrix& matrix,
                           const BlockMatrix& inverse,
                           const Numeric& d) {
  jacobian_targetsAddWindField(jacobian_targets, component, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = to_wind(component)}, matrix, inverse);
}

void RetrievalAddTemperature(JacobianTargets& jacobian_targets,
                             JacobianTargetsDiagonalCovarianceMatrixMap&
                                 covariance_matrix_diagonal_blocks,
                             const BlockMatrix& matrix,
                             const BlockMatrix& inverse,
                             const Numeric& d) {
  jacobian_targetsAddTemperature(jacobian_targets, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{AtmKey::t}}, matrix, inverse);
}

void RetrievalAddPressure(JacobianTargets& jacobian_targets,
                          JacobianTargetsDiagonalCovarianceMatrixMap&
                              covariance_matrix_diagonal_blocks,
                          const BlockMatrix& matrix,
                          const BlockMatrix& inverse,
                          const Numeric& d) {
  jacobian_targetsAddPressure(jacobian_targets, d);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{AtmKey::p}}, matrix, inverse);
}
