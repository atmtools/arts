#include <covariance_matrix.h>
#include <retrieval_target.h>
#include <workspace.h>

#include <concepts>

template <typename T>
concept CovarMatrixType = std::same_as<T, Matrix> or std::same_as<T, Sparse>;

template <Jacobian::target_type T>
void add_diagonal_covmat(CovarianceMatrix& covmat,
                         const T& target,
                         const BlockMatrix& forward,
                         const BlockMatrix& inverse) {
  Range colrow(target.x_start, target.x_size);

  ARTS_USER_ERROR_IF(
      forward.ncols() != colrow.extent or forward.nrows() != colrow.extent,
      "The forward matrix must be square.  It must also have the same size as the target.");

  covmat.add_correlation({colrow,
                          colrow,
                          IndexPair{target.target_pos, target.target_pos},
                          forward});

  if (inverse.not_null()) {
    ARTS_USER_ERROR_IF(
        inverse.ncols() != colrow.extent or inverse.nrows() != colrow.extent,
        "The inverse matrix must be square.  It must also have the same size as the target.");

    covmat.add_correlation({colrow,
                            colrow,
                            IndexPair{target.target_pos, target.target_pos},
                            forward});
  }
}

void model_state_covariance_matrixInit(CovarianceMatrix& covmat) {
  covmat = CovarianceMatrix{};
}

void model_state_covariance_matrixAdd(CovarianceMatrix& covmat,
                                      const JacobianTargets& jacobian_targets,
                                      const AtmKeyVal& new_target,
                                      const BlockMatrix& forward,
                                      const BlockMatrix& inverse) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.atm()) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(covmat, target, forward, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for atmospheric target : ", new_target);
}

void model_state_covariance_matrixAdd(CovarianceMatrix& covmat,
                                      const JacobianTargets& jacobian_targets,
                                      const SurfaceKeyVal& new_target,
                                      const BlockMatrix& forward,
                                      const BlockMatrix& inverse) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.surf()) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(covmat, target, forward, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : ", new_target);
}

void model_state_covariance_matrixAdd(CovarianceMatrix& covmat,
                                      const JacobianTargets& jacobian_targets,
                                      const LblLineKey& new_target,
                                      const BlockMatrix& forward,
                                      const BlockMatrix& inverse) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.line()) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(covmat, target, forward, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : ", new_target);
}

void model_state_covariance_matrixAddSpeciesVMR(
    CovarianceMatrix& covmat,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& species,
    const BlockMatrix& forward,
    const BlockMatrix& inverse) {
  model_state_covariance_matrixAdd(
      covmat, jacobian_targets, AtmKeyVal{species}, forward, inverse);
}

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

void RetrievalFinalizeDiagonal(
    CovarianceMatrix& covmat,
    JacobianTargets& jacobian_targets,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const ArrayOfAbsorptionBand& absorption_bands,
    const JacobianTargetsDiagonalCovarianceMatrixMap& covariance_matrix_diagonal_blocks) {
  jacobian_targetsFinalize(
      jacobian_targets, atmospheric_field, surface_field, absorption_bands);

  for (auto& key_data : covariance_matrix_diagonal_blocks.map) {
    std::visit(
        [&](auto& k) {
          model_state_covariance_matrixAdd(
              covmat, jacobian_targets, k, key_data.second.first, key_data.second.second);
        },
        key_data.first.target);
  }
}