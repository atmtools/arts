#include <covariance_matrix.h>
#include <retrieval_target.h>
#include <workspace.h>

#include "atm.h"
#include "jacobian.h"

////////////////////////////////////////////////////////////////////////////////
// Measurement vector error covariance matrix
////////////////////////////////////////////////////////////////////////////////

template <Jacobian::target_type T>
void add_diagonal_covmat(CovarianceMatrix& covmat,
                         const T& target,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  Range colrow(target.x_start, target.x_size);

  ARTS_USER_ERROR_IF(
      matrix.ncols() != colrow.extent or matrix.nrows() != colrow.extent,
      "The matrix must be square.  It must also have the same size as the target.");

  covmat.add_correlation({colrow,
                          colrow,
                          IndexPair{target.target_pos, target.target_pos},
                          matrix});

  if (inverse.not_null()) {
    ARTS_USER_ERROR_IF(
        inverse.ncols() != colrow.extent or inverse.nrows() != colrow.extent,
        "The inverse matrix must be square.  It must also have the same size as the target.");

    covmat.add_correlation({colrow,
                            colrow,
                            IndexPair{target.target_pos, target.target_pos},
                            matrix});
  }
}

void model_state_covariance_matrixInit(
    CovarianceMatrix& model_state_covariance_matrix) {
  model_state_covariance_matrix = CovarianceMatrix{};
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const AtmKeyVal& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.atm()) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for atmospheric target : ", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SurfaceKeyVal& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.surf()) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : ", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const LblLineKey& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.line()) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : ", new_target);
}

void model_state_covariance_matrixAddSpeciesVMR(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& species,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKeyVal{species},
                                   matrix,
                                   inverse);
}

////////////////////////////////////////////////////////////////////////////////
// Measurement vector error covariance matrix
////////////////////////////////////////////////////////////////////////////////

void measurement_vector_error_covariance_matrixConstant(
    CovarianceMatrix& measurement_vector_error_covariance_matrix,
    const ArrayOfSensorObsel& measurement_sensor,
    const Numeric& x) {
  measurement_vector_error_covariance_matrix = CovarianceMatrix{};

  const Size N = measurement_sensor.size();

  measurement_vector_error_covariance_matrix.add_correlation(
      {Range(0, N),
       Range(0, N),
       IndexPair{0, 0},
       std::make_shared<Sparse>(Sparse::diagonal(Vector(N, x)))});
}

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void RetrievalInit(JacobianTargets& jacobian_targets,
                   CovarianceMatrix& model_state_covariance_matrix,
                   JacobianTargetsDiagonalCovarianceMatrixMap&
                       covariance_matrix_diagonal_blocks) {
  model_state_covariance_matrixInit(model_state_covariance_matrix);
  jacobian_targets.clear();
  covariance_matrix_diagonal_blocks.clear();
}

void RetrievalAddSpeciesVMR(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& species,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse,
                            const Numeric& dvmr) {
  jacobian_targetsAddSpeciesVMR(jacobian_targets, species, dvmr);
  covariance_matrix_diagonal_blocks.set(
      {.target = AtmKeyVal{species}}, matrix, inverse);
}

void RetrievalFinalizeDiagonal(CovarianceMatrix& model_state_covariance_matrix,
                               JacobianTargets& jacobian_targets,
                               const JacobianTargetsDiagonalCovarianceMatrixMap&
                                   covariance_matrix_diagonal_blocks,
                               const AtmField& atmospheric_field,
                               const SurfaceField& surface_field,
                               const ArrayOfAbsorptionBand& absorption_bands) {
  jacobian_targetsFinalize(
      jacobian_targets, atmospheric_field, surface_field, absorption_bands);

  for (auto& key_data : covariance_matrix_diagonal_blocks) {
    std::visit(
        [&](auto& k) {
          model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                           jacobian_targets,
                                           k,
                                           key_data.second.first,
                                           key_data.second.second);
        },
        key_data.first.target);
  }
}