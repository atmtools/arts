#include <covariance_matrix.h>

#include <concepts>

#include "debug.h"
#include "jacobian.h"
#include "matpack_data.h"
#include "matpack_sparse.h"
#include "obsel.h"

template <Jacobian::target_type T, typename M>
void add_diagonal_covmat(CovarianceMatrix& covmat,
                         const T& target,
                         const M& data)
  requires(std::same_as<M, Matrix> or std::same_as<M, Sparse>)
{
  Range colrow(target.x_start, target.x_size);

  ARTS_USER_ERROR_IF(
      data.ncols() != colrow.extent or data.nrows() != colrow.extent,
      "The data matrix must be square.  It must also have the same size as the target.");

  covmat.add_correlation({colrow,
                          colrow,
                          IndexPair{target.target_pos, target.target_pos},
                          std::make_shared<M>(data)});
}

void model_state_covariance_matrixInit(CovarianceMatrix& covmat) {
  covmat = CovarianceMatrix{};
}

template <typename M>
void _model_state_covariance_matrixAddSpeciesVMR(
    CovarianceMatrix& covmat,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& species,
    const M& data) {
  ARTS_USER_ERROR_IF(not jacobian_targets.finalized,
                     "Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jacobian_targets.atm()) {
    if (target.type == species) {
      found = true;
      add_diagonal_covmat(covmat, target, data);
    }
  }

  ARTS_USER_ERROR_IF(not found, "No target found for species : ", species);
}

void model_state_covariance_matrixAddSpeciesVMR(
    CovarianceMatrix& covmat,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& species,
    const Sparse& data) {
  _model_state_covariance_matrixAddSpeciesVMR(
      covmat, jacobian_targets, species, data);
}

void model_state_covariance_matrixAddSpeciesVMR(
    CovarianceMatrix& covmat,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& species,
    const Matrix& data) {
  _model_state_covariance_matrixAddSpeciesVMR(
      covmat, jacobian_targets, species, data);
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
