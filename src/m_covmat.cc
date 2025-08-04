#include <workspace.h>

#include "atm.h"
#include "enumsFieldComponent.h"
#include "enumsSpeciesEnum.h"
#include "isotopologues.h"

////////////////////////////////////////////////////////////////////////////////
// Measurement vector error covariance matrix
////////////////////////////////////////////////////////////////////////////////

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

void model_state_covariance_matrixInit(
    CovarianceMatrix& model_state_covariance_matrix) {
  ARTS_TIME_REPORT

  model_state_covariance_matrix = CovarianceMatrix{};
}

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

void model_state_covariance_matrixAddAtmosphere(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const AtmKey& key,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKeyVal{key},
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddAtmosphere(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& key,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKeyVal{key},
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddAtmosphere(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SpeciesIsotope& key,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKeyVal{key},
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddAtmosphere(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const QuantumLevelIdentifier& key,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKeyVal{key},
                                   matrix,
                                   inverse);
}

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

void model_state_covariance_matrixAddIsotopologueRatio(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const SpeciesIsotope& species,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKeyVal{species},
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddMagneticField(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const String& component,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   to_mag(component),
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddWindField(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const String& component,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   to_wind(component),
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddTemperature(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKey::t,
                                   matrix,
                                   inverse);
}

void model_state_covariance_matrixAddPressure(
    CovarianceMatrix& model_state_covariance_matrix,
    const JacobianTargets& jacobian_targets,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  ARTS_TIME_REPORT

  model_state_covariance_matrixAdd(model_state_covariance_matrix,
                                   jacobian_targets,
                                   AtmKey::p,
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
                               const AbsorptionBands& absorption_bands,
                               const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  jacobian_targetsFinalize(jacobian_targets,
                           atmospheric_field,
                           surface_field,
                           subsurface_field,
                           absorption_bands,
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
