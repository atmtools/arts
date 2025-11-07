#include <enumsSpeciesEnum.h>
#include <isotopologues.h>
#include <workspace.h>

////////////////////////////////////////////////////////////////////////////////
// Measurement vector error covariance matrix
////////////////////////////////////////////////////////////////////////////////

namespace {
template <Jacobian::target_type T>
void add_diagonal_covmat(CovarianceMatrix& covmat,
                         const T& target,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  const Range colrow(target.x_start, target.x_size);

  ARTS_USER_ERROR_IF(
      matrix.ncols() != colrow.nelem or matrix.nrows() != colrow.nelem,
      R"(The matrix must be square.  It must also have the same size as the target.
     shape(matrix) = {:B,},
     shape(target) = [{}, {}]
Target: {}
)",
      matrix.shape(),
      colrow.nelem,
      colrow.nelem,
      target.type);

  if (not target.overlap) {
    covmat.add_correlation({colrow,
                            colrow,
                            IndexPair{target.target_pos, target.target_pos},
                            matrix});
  }

  if (inverse.not_null()) {
    ARTS_USER_ERROR_IF(
        inverse.ncols() != colrow.nelem or inverse.nrows() != colrow.nelem,
        R"(The inverse matrix must be square.  It must also have the same size as the target.
     shape(matrix) = {:B,},
     shape(target) = [{}, {}]
Target: {}
)",
        inverse.shape(),
        colrow.nelem,
        colrow.nelem,
        target.type);

    if (not target.overlap) {
      covmat.add_correlation({colrow,
                              colrow,
                              IndexPair{target.target_pos, target.target_pos},
                              matrix});
    }
  }
}
}  // namespace

void model_state_covariance_matrixInit(
    CovarianceMatrix& model_state_covariance_matrix) {
  ARTS_TIME_REPORT

  model_state_covariance_matrix = CovarianceMatrix{};
}

namespace {
void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const AtmKeyVal& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.atm) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for atmospheric target : {}", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SurfaceKeyVal& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.surf) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : {}", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SubsurfaceKeyVal& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.subsurf) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : {}", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const LblLineKey& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.line) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for surface target : {}", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SensorKey& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.sensor) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for sensor target : {}", new_target);
}

void model_state_covariance_matrixAdd(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const ErrorKey& new_target,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.error) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(
          model_state_covariance_matrix, target, matrix, inverse);
    }
  }

  ARTS_USER_ERROR_IF(
      not found, "No target found for sensor target : {}", new_target);
}
}  // namespace

void model_state_covariance_matrixAddSpeciesVMR(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& species,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

  measurement_vector_error_covariance_matrix = CovarianceMatrix{};

  const Size N = measurement_sensor.size();

  measurement_vector_error_covariance_matrix.add_correlation(
      {Range(0, N),
       Range(0, N),
       IndexPair{0, 0},
       std::make_shared<Sparse>(Sparse::diagonal(Vector(N, x)))});
  measurement_vector_error_covariance_matrix.add_correlation_inverse(
      {Range(0, N),
       Range(0, N),
       IndexPair{0, 0},
       std::make_shared<Sparse>(Sparse::diagonal(Vector(N, 1.0 / x)))});
}

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void RetrievalFinalizeDiagonal(CovarianceMatrix& model_state_covariance_matrix,
                               JacobianTargets& jacobian_targets,
                               const JacobianTargetsDiagonalCovarianceMatrixMap&
                                   covariance_matrix_diagonal_blocks,
                               const AtmField& atmospheric_field,
                               const SurfaceField& surface_field,
                               const SubsurfaceField& subsurface_field,
                               const AbsorptionBands& abs_bands,
                               const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surface_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surface_field.ellipsoid)

  jacobian_targetsFinalize(jacobian_targets,
                           atmospheric_field,
                           surface_field,
                           subsurface_field,
                           abs_bands,
                           measurement_sensor);

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
