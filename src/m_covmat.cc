#include <enumsSpeciesEnum.h>
#include <isotopologues.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <utility>

////////////////////////////////////////////////////////////////////////////////
// Measurement vector error covariance matrix
////////////////////////////////////////////////////////////////////////////////

namespace {
template <Jacobian::target_type T> void add_diagonal_covmat(CovarianceMatrix&  covmat,
                                                            const T&           target,
                                                            const BlockMatrix& matrix,
                                                            const BlockMatrix& inverse) {
  const Range colrow(target.x_start, target.x_size);

  if (not matrix.not_null()) arts::user_error("The covariance matrix is null for target {}.", target.type);

  if (matrix.ncols() != colrow.nelem or matrix.nrows() != colrow.nelem)
    arts::user_error(R"(The matrix must be square.  It must also have the same size as the target.
     shape(matrix) = {:B,},
     shape(target) = [{}, {}]
Target: {}
)",
                     matrix.shape(),
                     colrow.nelem,
                     colrow.nelem,
                     target.type);

  if (inverse.not_null()) {
    if (inverse.ncols() != colrow.nelem or inverse.nrows() != colrow.nelem)
      arts::user_error(R"(The inverse matrix must be square.  It must also have the same size as the target.
     shape(matrix) = {:B,},
     shape(target) = [{}, {}]
Target: {}
)",
                       inverse.shape(),
                       colrow.nelem,
                       colrow.nelem,
                       target.type);
  }

  if (not target.overlap) {
    covmat.add_correlation({colrow, colrow, IndexPair{target.target_pos, target.target_pos}, matrix});
    if (inverse.not_null())
      covmat.add_correlation_inverse({colrow, colrow, IndexPair{target.target_pos, target.target_pos}, inverse});
  }
}
}  // namespace

void oemStateCovmatInit(OptimalEstimationData& data) {
  data.uncheck();
  auto& model_state_covmat = data.model_state_covmat;
  ARTS_TIME_REPORT

  model_state_covmat = CovarianceMatrix{};
}

namespace {
void model_state_covmatAdd(CovarianceMatrix&      model_state_covmat,
                           const JacobianTargets& jac_targets,
                           const AtmKeyVal&       new_target,
                           const BlockMatrix&     matrix,
                           const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT

  if (not jac_targets.finalized) arts::user_error("Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jac_targets.atm) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(model_state_covmat, target, matrix, inverse);
    }
  }

  if (not found) arts::user_error("No target found for atmospheric target : {}", new_target);
}

void model_state_covmatAdd(CovarianceMatrix&      model_state_covmat,
                           const JacobianTargets& jac_targets,
                           const SurfaceKeyVal&   new_target,
                           const BlockMatrix&     matrix,
                           const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT

  if (not jac_targets.finalized) arts::user_error("Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jac_targets.surf) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(model_state_covmat, target, matrix, inverse);
    }
  }

  if (not found) arts::user_error("No target found for surface target : {}", new_target);
}

void model_state_covmatAdd(CovarianceMatrix&       model_state_covmat,
                           const JacobianTargets&  jac_targets,
                           const SubsurfaceKeyVal& new_target,
                           const BlockMatrix&      matrix,
                           const BlockMatrix&      inverse) {
  ARTS_TIME_REPORT

  if (not jac_targets.finalized) arts::user_error("Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jac_targets.subsurf) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(model_state_covmat, target, matrix, inverse);
    }
  }

  if (not found) arts::user_error("No target found for surface target : {}", new_target);
}

void model_state_covmatAdd(CovarianceMatrix&      model_state_covmat,
                           const JacobianTargets& jac_targets,
                           const LblLineKey&      new_target,
                           const BlockMatrix&     matrix,
                           const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT

  if (not jac_targets.finalized) arts::user_error("Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jac_targets.line) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(model_state_covmat, target, matrix, inverse);
    }
  }

  if (not found) arts::user_error("No target found for surface target : {}", new_target);
}

void model_state_covmatAdd(CovarianceMatrix&      model_state_covmat,
                           const JacobianTargets& jac_targets,
                           const SensorKey&       new_target,
                           const BlockMatrix&     matrix,
                           const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT

  if (not jac_targets.finalized) arts::user_error("Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jac_targets.sensor) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(model_state_covmat, target, matrix, inverse);
    }
  }

  if (not found) arts::user_error("No target found for sensor target : {}", new_target);
}

void model_state_covmatAdd(CovarianceMatrix&      model_state_covmat,
                           const JacobianTargets& jac_targets,
                           const ErrorKey&        new_target,
                           const BlockMatrix&     matrix,
                           const BlockMatrix&     inverse) {
  ARTS_TIME_REPORT

  if (not jac_targets.finalized) arts::user_error("Jacobian targets not finalized.");

  bool found = false;

  for (const auto& target : jac_targets.error) {
    if (target.type == new_target) {
      found = true;
      add_diagonal_covmat(model_state_covmat, target, matrix, inverse);
    }
  }

  if (not found) arts::user_error("No target found for sensor target : {}", new_target);
}
}  // namespace

////////////////////////////////////////////////////////////////////////////////
// Measurement vector error covariance matrix
////////////////////////////////////////////////////////////////////////////////

namespace {
void measurement_covmat_constant(CovarianceMatrix& covariance, Size n, Numeric variance) {
  if (not std::isfinite(variance) or variance <= 0 or not std::isfinite(1.0 / variance))
    arts::user_error(
        "Measurement-error variance must be finite and strictly positive with a finite reciprocal, got {}.", variance);
  if (n == 0) arts::user_error("A nonempty measurement vector is required.");
  covariance = CovarianceMatrix{};
  covariance.add_correlation(
      {Range(0, n), Range(0, n), IndexPair{0, 0}, std::make_shared<Sparse>(Sparse::diagonal(Vector(n, variance)))});
  covariance.add_correlation_inverse({Range(0, n),
                                      Range(0, n),
                                      IndexPair{0, 0},
                                      std::make_shared<Sparse>(Sparse::diagonal(Vector(n, 1.0 / variance)))});
}
}  // namespace

void measurement_vec_error_covmatConstant(CovarianceMatrix& covariance,
                                          const Vector&     measurement_vec,
                                          const Numeric&    value) {
  measurement_covmat_constant(covariance, measurement_vec.size(), value);
}

void oemMeasurementCovmatConstant(OptimalEstimationData& data, const Numeric& value) {
  measurement_covmat_constant(data.measurement_vec_error_covmat, data.measurement_vec.size(), value);
}

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void oemFinalizeDiagonal(OptimalEstimationData&    oem,
                         JacobianTargets&          jac_targets,
                         const AtmField&           atm_field,
                         const SurfaceField&       surf_field,
                         const SubsurfaceField&    subsurf_field,
                         const AbsorptionBands&    abs_bands,
                         const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT
  oem.uncheck();
  CovarianceMatrix model_state_covmat;
  if (oem.covmat_diagonal_blocks.empty())
    arts::user_error("No pending covariance blocks; add targets before finalizing.");

  jac_targetsFinalize(jac_targets, atm_field, surf_field, subsurf_field, abs_bands, measurement_sensor);

  for (auto& key_data : oem.covmat_diagonal_blocks) {
    std::visit(
        [&](auto& k) {
          model_state_covmatAdd(model_state_covmat, jac_targets, k, key_data.second.first, key_data.second.second);
        },
        key_data.first.target);
  }
  oem.model_state_covmat = std::move(model_state_covmat);
  oem.check(&jac_targets);
}

namespace {
void correlate_atmosphere(CovarianceMatrix&      covariance,
                          const JacobianTargets& targets,
                          const AtmField&        atmosphere,
                          const AtmKeyVal&       key1,
                          const AtmKeyVal&       key2,
                          Numeric                correlation) {
  if (not std::isfinite(correlation) or std::abs(correlation) >= 1)
    arts::user_error("correlation must be finite and strictly between -1 and 1.");
  if (key1 == key2) arts::user_error("Choose two distinct atmospheric targets.");
  if (not targets.finalized) arts::user_error("Jacobian targets must be finalized.");
  const auto find_target = [&](const AtmKeyVal& key) -> const Jacobian::AtmTarget& {
    const Jacobian::AtmTarget* result = nullptr;
    for (const auto& target : targets.atm) {
      if (target.type == key) {
        if (result != nullptr) arts::user_error("Ambiguous atmospheric target: {}", key);
        result = &target;
      }
    }
    if (result == nullptr) arts::user_error("No atmospheric retrieval target: {}", key);
    if (result->overlap) arts::user_error("Overlapping targets are not supported: {}", key);
    return *result;
  };
  const auto& first  = find_target(key1);
  const auto& second = find_target(key2);
  const auto* grid1  = atmosphere[key1].get_if<GeodeticField3>();
  const auto* grid2  = atmosphere[key2].get_if<GeodeticField3>();
  if (not grid1 or not grid2) arts::user_error("Both targets must have gridded atmospheric fields.");
  if (not std::ranges::equal(grid1->grid<0>(), grid2->grid<0>()) or
      not std::ranges::equal(grid1->grid<1>(), grid2->grid<1>()) or
      not std::ranges::equal(grid1->grid<2>(), grid2->grid<2>()))
    arts::user_error("Targets must use identical altitude, latitude, and longitude grids.");
  if (first.x_size != second.x_size or first.x_size != grid1->data.size())
    arts::user_error("Targets must have one state element per grid point.");
  covariance.validate(targets.x_size());
  const auto& blocks    = std::as_const(covariance).get_blocks();
  const auto  variances = [&](const Jacobian::AtmTarget& target) {
    const auto it = std::ranges::find_if(blocks, [&](const Block& block) {
      return block.get_indices() == IndexPair{target.target_pos, target.target_pos};
    });
    if (it == blocks.end()) arts::user_error("Missing marginal covariance for target {}", target.type);
    if (it->is_sparse()) {
      const auto& matrix = it->get_sparse().matrix;
      for (Index row = 0; row < matrix.outerSize(); ++row)
        for (Eigen::SparseMatrix<Numeric, Eigen::RowMajor>::InnerIterator entry(matrix, row); entry; ++entry) {
          if (entry.row() != entry.col() and entry.value() != 0)
            arts::user_error("Marginal covariance must be diagonal for target {}", target.type);
        }
    } else {
      const auto& matrix = it->get_dense();
      for (Index i = 0; i < matrix.nrows(); ++i)
        for (Index j = 0; j < matrix.ncols(); ++j) {
          if ((i != j and matrix[i, j] != 0))
            arts::user_error("Marginal covariance must be diagonal for target {}", target.type);
        }
    }
    return it->diagonal();
  };
  const Vector v1 = variances(first), v2 = variances(second);
  Vector       cross(v1.size());
  for (Size i = 0; i < cross.size(); ++i) cross[i] = correlation * std::sqrt(v1[i]) * std::sqrt(v2[i]);
  const auto&     row = first.target_pos < second.target_pos ? first : second;
  const auto&     col = first.target_pos < second.target_pos ? second : first;
  const IndexPair indices{row.target_pos, col.target_pos};
  auto            replacement = blocks;
  std::erase_if(replacement, [&](const Block& block) { return block.get_indices() == indices; });
  if (correlation != 0) {
    replacement.emplace_back(
        Range(row.x_start, row.x_size), Range(col.x_start, col.x_size), indices, Sparse::diagonal(cross));
  }
  CovarianceMatrix candidate;
  candidate.set_blocks(std::move(replacement));
  // Other cross-correlations can make an individually valid pair inconsistent.
  // Validate the entire candidate before committing; inverse caches are discarded.
  candidate.validate(targets.x_size());
  covariance = std::move(candidate);
}
}  // namespace

using GenericAtmKey = Generic<const AtmKey,
                              const QuantumLevelIdentifier,
                              const ScatteringSpeciesProperty,
                              const SpeciesEnum,
                              const SpeciesIsotope>;

void oemStateCovmatCorrelateConstant(OptimalEstimationData& data,
                                     const JacobianTargets& targets,
                                     const AtmField&        atmosphere,
                                     const GenericAtmKey    target1,
                                     const GenericAtmKey    target2,
                                     const Numeric&         correlation) {
  auto& covariance = data.model_state_covmat;
  if (not data.model_state_vec_apriori.empty() and
      covariance.nrows() != static_cast<Index>(data.model_state_vec_apriori.size()))
    arts::user_error("State covariance size does not match the OEM prior.");
  const auto key = [](const auto& target) { return std::visit([](const auto& p) -> AtmKeyVal { return *p; }, target); };
  correlate_atmosphere(covariance, targets, atmosphere, key(target1), key(target2), correlation);
}

void oemMeasurementCovmatInit(OptimalEstimationData& data) {
  data.uncheck();
  data.measurement_vec_error_covmat = CovarianceMatrix{};
}

void oemMeasurementCovmatAdd(OptimalEstimationData& data, const BlockMatrix& matrix, const BlockMatrix& inverse) {
  auto&       covariance = data.measurement_vec_error_covmat;
  const Index start      = covariance.nrows();
  const Index size       = matrix.nrows();
  if (not matrix.not_null() or size == 0 or matrix.ncols() != size)
    arts::user_error("Measurement covariance block must be nonempty and square.");
  if (not data.measurement_vec.empty() and start + size > static_cast<Index>(data.measurement_vec.size()))
    arts::user_error(
        "Adding {} covariance rows exceeds the {} initialized measurements.", size, data.measurement_vec.size());
  if (inverse.not_null() and (inverse.nrows() != size or inverse.ncols() != size))
    arts::user_error("Inverse block dimensions must match the covariance block.");
  const Index block = covariance.ndiagblocks();
  const Range range(start, size);
  data.uncheck();
  covariance.add_correlation(Block(range, range, {block, block}, matrix));
  if (inverse.not_null()) covariance.add_correlation_inverse(Block(range, range, {block, block}, inverse));
}
