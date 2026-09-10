#include <workspace.h>

#include <cmath>
#include <thread>
#include <atomic>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string_view>

namespace {

void require(bool condition, std::string_view message) {
  if (not condition) throw std::runtime_error(std::string(message));
}

void rejects(const std::function<void()>& operation, std::string_view name) {
  try {
    operation();
  } catch (const std::exception& error) {
    require(not std::string_view(error.what()).empty(), "Validation error has no explanation");
    return;
  }
  throw std::runtime_error(std::format("Accepted invalid covariance: {}", name));
}

Matrix matrix(Index rows, Index cols, std::initializer_list<Numeric> values) {
  Matrix result(rows, cols);
  auto   value = values.begin();
  for (Index i = 0; i < rows; ++i)
    for (Index j = 0; j < cols; ++j) result[i, j] = *value++;
  return result;
}

Block block(Index offset, Index index, Numeric variance) {
  return {Range(offset, 1), Range(offset, 1), {index, index}, matrix(1, 1, {variance})};
}

CovarianceMatrix covariance(Matrix value) {
  const Index      n = value.nrows();
  CovarianceMatrix result;
  result.add_correlation({Range(0, n), Range(0, n), {0, 0}, std::move(value)});
  return result;
}

const std::vector<Block>& inverse_blocks(const CovarianceMatrix& value) { return value.get_inverse_blocks(); }

void close(Numeric actual, Numeric expected, std::string_view name) {
  require(std::isfinite(actual) and std::abs(actual - expected) <= 2e-12 * std::max(std::abs(expected), 1e-100), name);
}

CovarianceMatrix correlated_blocks() {
  CovarianceMatrix result;
  result.add_correlation(block(0, 0, 4));
  result.add_correlation(block(1, 1, 2));
  result.add_correlation({Range(0, 1), Range(1, 1), {0, 1}, matrix(1, 1, {1})});
  return result;
}

void structural() {
  rejects([] { CovarianceMatrix{}.validate(); }, "empty matrix");
  auto valid = covariance(matrix(2, 2, {4, 1, 1, 2}));
  valid.validate(2);
  rejects([&] { valid.validate(3); }, "expected dimension mismatch");

  for (const auto& malformed : std::vector<std::vector<Block>>{
           {Block{}},
           {{Range(0, 2), Range(0, 2), {0, 0}, matrix(1, 1, {1})}},
           {{Range(-1, 1), Range(-1, 1), {0, 0}, matrix(1, 1, {1})}},
           {{Range(0, -1), Range(0, -1), {0, 0}, matrix(1, 1, {1})}},
           {{Range(std::numeric_limits<Index>::max(), 2), Range(0, 2), {0, 0}, matrix(2, 2, {1, 0, 0, 1})}},
           {block(1, 0, 1)},
           {block(0, 0, 1), block(2, 1, 1)},
           {block(0, 0, 1), block(0, 1, 1)},
           {block(0, 0, 1), block(0, 0, 1)},
           {block(0, -1, 1)},
           {block(0, 0, 1), block(1, 1, 1), {Range(1, 1), Range(1, 1), {0, 1}, matrix(1, 1, {0.1})}},
           {block(0, 0, 1), {Range(0, 1), Range(1, 1), {0, 1}, matrix(1, 1, {0.1})}},
           {block(0, 0, 1), block(1, 1, 1), {Range(1, 1), Range(0, 1), {1, 0}, matrix(1, 1, {0.1})}},
       }) {
    CovarianceMatrix candidate;
    candidate.set_blocks(malformed);
    rejects([&] { candidate.validate(); }, "malformed block storage");
    rejects([&] { candidate.prepared(); }, "malformed prepared block storage");
  }

  auto duplicate = correlated_blocks();
  duplicate.add_correlation({Range(0, 1), Range(1, 1), {0, 1}, matrix(1, 1, {1})});
  rejects([&] { duplicate.validate(); }, "duplicate off-diagonal block");
  for (Numeric tolerance :
       {0.0, -1.0, std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::quiet_NaN()})
    rejects([&] { valid.validate(-1, tolerance); }, "invalid relative tolerance");
}

void numerical() {
  for (const auto& invalid : std::vector<Matrix>{
           matrix(2, 2, {-1, 0, 0, -2}),
           matrix(2, 2, {1, 2, 2, 1}),
           matrix(2, 2, {1, 1, 1, 1}),
           matrix(2, 2, {1, 0.1, 0.2, 1}),
           matrix(2, 2, {0, 0, 0, 1}),
           matrix(2, 2, {1, std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::infinity(), 1}),
           matrix(1, 1, {std::numeric_limits<Numeric>::quiet_NaN()}),
       }) {
    auto candidate = covariance(invalid);
    rejects([&] { candidate.validate(); }, "invalid covariance values");
    rejects([&] { candidate.compute_inverse(); }, "invalid covariance inversion");
    rejects([&] { candidate.prepared(); }, "invalid prepared covariance");
  }

  // A positive variance can be tiny in physical units without losing information.
  auto tiny = covariance(matrix(2, 2, {4e-40, 1e-20, 1e-20, 2}));
  tiny.validate(2);
  require(inverse_blocks(tiny).empty(), "Validation computed an inverse");
  tiny.compute_inverse();
  tiny.validate(2);
  close(tiny.get_inverse()[0, 0], 2e40 / 7, "Small physical units broke inversion");

  auto connected = correlated_blocks();
  connected.validate(2);
  rejects([&] { connected.validate(2, 1e-10, 1); }, "dense component exceeds allocation limit");
  connected.compute_inverse();
  connected.validate(2);
  close(connected.get_inverse()[0, 1], -1.0 / 7, "Connected covariance inverse lost correlation");

  // Block IDs identify targets, independently of their coordinate order.
  // The stored upper block triangle can occupy the lower physical triangle.
  CovarianceMatrix reversed;
  reversed.add_correlation(block(1, 0, 2));
  reversed.add_correlation(block(0, 2, 4));
  reversed.add_correlation({Range(1, 1), Range(0, 1), {0, 2}, matrix(1, 1, {1})});
  reversed.validate(2);
  reversed.compute_inverse();
  reversed.validate(2);
  const Matrix reversed_inverse = reversed.get_inverse();
  close(reversed_inverse[0, 0], 2.0 / 7, "Permuted block IDs changed first precision diagonal");
  close(reversed_inverse[1, 1], 4.0 / 7, "Permuted block IDs changed second precision diagonal");
  close(reversed_inverse[0, 1], -1.0 / 7, "Permuted block IDs lost precision correlation");

  auto indefinite                                  = correlated_blocks();
  indefinite.get_blocks().back().get_dense()[0, 0] = 3;
  rejects([&] { indefinite.validate(); }, "globally indefinite connected blocks");

  // Direct diagonal solves do not need a reciprocal cache. In particular,
  // this finite quotient must work even though 1/variance overflows.
  CovarianceMatrix tiny_diagonal;
  tiny_diagonal.add_correlation({Range(0, 2), Range(0, 2), {0, 0},
                                Sparse::diagonal(Vector{1e-310, 2.})});
  Vector rhs{1e-310, 6.}, solution(2);
  solve(solution, tiny_diagonal, rhs);
  close(solution[0], 1., "Diagonal division overflowed via reciprocal");
  close(solution[1], 3., "Diagonal vector solve");
  Matrix rhs_matrix(2, 3, 0.), result(2, 3);
  for (Index j = 0; j < 3; ++j) {
    rhs_matrix[0,j] = 1e-310;
    rhs_matrix[1,j] = 6.;
  }
  mult_inv(result, tiny_diagonal, rhs_matrix);
  for (Index j = 0; j < 3; ++j) {
    close(result[0,j], 1., "Left diagonal solve");
    close(result[1,j], 3., "Left diagonal solve");
  }
  Matrix right(3, 2);
  mult_inv(right, transpose(rhs_matrix), tiny_diagonal);
  for (Index j = 0; j < 3; ++j) {
    close(right[j,0], 1., "Right diagonal solve");
    close(right[j,1], 3., "Right diagonal solve");
  }
  require(inverse_blocks(tiny_diagonal).empty(), "Diagonal solve constructed inverse blocks");

  // This path must preserve sparse diagonal scaling, rather than allocate a
  // dense 10000 by 10000 matrix or run a cubic factorization.
  constexpr Index  n = 10000;
  CovarianceMatrix diagonal;
  diagonal.add_correlation({Range(0, n), Range(0, n), {0, 0}, Sparse::diagonal(Vector(n, 1e-30))});
  diagonal.validate(n);
  diagonal.validate(n, 1e-10, 1);
  diagonal.compute_inverse();
  diagonal.validate(n);
  require(inverse_blocks(diagonal).size() == 1 and inverse_blocks(diagonal)[0].is_sparse(),
          "Densified diagonal inverse");
  close(inverse_blocks(diagonal)[0].get_sparse().ro(0, 0), 1e30, "Sparse diagonal precision");
}

void structured_solves() {
  CovarianceMatrix a;
  a.add_correlation(block(0,0,2));
  a.add_correlation({Range(1,1),Range(1,1),{1,1},Sparse::diagonal(Vector{4})});
  a.add_correlation(block(2,2,2));
  a.add_correlation(block(3,3,5));
  Vector rhs{2,4,2,5}, out(4);
  solve(out,a,rhs);
  for(Index i=0;i<4;++i) close(out[i],1,"Initial independent diagonal solve");
  a.add_correlation({Range(1,1),Range(2,1),{1,2},matrix(1,1,{1})});
  rhs=Vector{2,5,3,5};
  solve(out,a,rhs);
  for(Index i=0;i<4;++i) close(out[i],1,"Joined component solve");
  require(inverse_blocks(a).empty(),"Structured solve formed inverse");
  // Keep a mutable reference across factor creation, then mutate through it.
  auto& blocks=a.get_blocks();
  solve(out,a,rhs);
  blocks.back().get_dense()[0,0]=0.5;
  rhs=Vector{2,4.5,2.5,5};
  solve(out,a,rhs);
  for(Index i=0;i<4;++i) close(out[i],1,"Retained reference invalidation");
  auto copy=a;
  blocks.pop_back();
  rhs=Vector{2,4,2,5};
  solve(out,a,rhs);
  for(Index i=0;i<4;++i) close(out[i],1,"Split component solve");
  rhs=Vector{2,4.5,2.5,5};
  solve(out,copy,rhs);
  for(Index i=0;i<4;++i) close(out[i],1,"Copied factor cache");
  Matrix r(4,2), x(4,2), right(2,4);
  for(Index i=0;i<4;++i) {r[i,0]=rhs[i];r[i,1]=2*rhs[i];}
  mult_inv(x,copy,r);
  mult_inv(right,transpose(r),copy);
  for(Index i=0;i<4;++i) for(Index j=0;j<2;++j) {
    close(x[i,j],Numeric(j+1),"Component matrix solve");
    close(right[j,i],Numeric(j+1),"Component right solve");
  }
  auto shared=std::make_shared<Matrix>(matrix(2,2,{2,0,0,3}));
  CovarianceMatrix external;
  external.add_correlation({Range(0,2),Range(0,2),{0,0},shared});
  Vector b{2,3}, y(2);
  solve(y,external,b);
  (*shared)[0,1]=(*shared)[1,0]=1;
  b=Vector{3,4};
  solve(y,external,b);
  close(y[0],1,"Shared matrix changed structure");
  close(y[1],1,"Shared matrix changed structure");
  (*shared)[0,0]=-1;
  rejects([&]{solve(y,external,b);},"Invalid mutation after cached solve");
}

void prepared_solves() {
  auto supplied = covariance(matrix(2, 2, {2, 0, 0, 4}));
  supplied.compute_inverse();
  supplied.prepared(true);
  supplied.get_inverse_blocks()[0].get_dense()[0, 0] = 3;
  rejects([&] { supplied.prepared(true); }, "Changed inverse reused validation proof");
  auto             storage = std::make_shared<Matrix>(matrix(2, 2, {2, 0, 0, 4}));
  CovarianceMatrix source;
  source.add_correlation({Range(0, 2), Range(0, 2), {0, 0}, storage});
  std::atomic<bool>         consistent{true};
  std::vector<std::jthread> readers;
  for (int i = 0; i < 8; ++i)
    readers.emplace_back([&] {
      const auto snapshot = source.prepared();
      Vector     rhs{2, 4}, out(2);
      for (int j = 0; j < 20; ++j) {
        solve(out, *snapshot, rhs);
        if (out[0] != 1 or out[1] != 1) consistent = false;
      }
    });
  readers.clear();
  require(consistent, "Concurrent preparation/solve mismatch");
  const auto first = source.prepared();
  require(first == source.prepared(), "Unchanged preparation was not reused");
  Vector rhs{2, 4}, out(2);
  solve(out, *first, rhs);
  close(out[0], 1, "Prepared diagonal solve");
  (*storage)[0, 0] = 4;
  solve(out, *first, rhs);
  close(out[0], 1, "Source alias changed immutable snapshot");
  const auto second = source.prepared();
  require(second != first, "Source mutation reused stale preparation");
  solve(out, *second, rhs);
  close(out[0], 0.5, "Updated preparation");
  auto copy                              = *second;
  copy.get_blocks()[0].get_dense()[0, 0] = 8;
  solve(out, *second, rhs);
  close(out[0], 0.5, "Mutable copy changed immutable snapshot");
  const auto precision = source.prepared(true);
  solve(out, *precision, rhs);
  close(out[0], 0.5, "Prepared explicit precision");
  require(precision == source.prepared(true), "Precision preparation not reused");
  storage->resize(1, 4);
  *storage = matrix(1, 4, {4, 0, 0, 4});
  rejects([&] { source.prepared(true); }, "Changed storage shape reused preparation");
}

void precision() {
  prepared_solves();
  structured_solves();
  auto valid = correlated_blocks();
  valid.add_correlation_inverse(block(0, 0, 2.0 / 7));
  valid.add_correlation_inverse(block(1, 1, 4.0 / 7));
  valid.add_correlation_inverse({Range(0, 1), Range(1, 1), {0, 1}, matrix(1, 1, {-1.0 / 7})});
  valid.validate(2);
  const Matrix before = valid.get_inverse();
  valid.compute_inverse();
  close(valid.get_inverse()[0, 1], before[0, 1], "Changed valid supplied precision");

  auto incorrect = correlated_blocks();
  incorrect.add_correlation_inverse(block(0, 0, 0.25));
  incorrect.add_correlation_inverse(block(1, 1, 0.5));
  rejects([&] { incorrect.validate(); }, "componentwise inversion ignored correlations");
  rejects([&] { incorrect.compute_inverse(); }, "inconsistent supplied precision");

  auto partial = correlated_blocks();
  partial.add_correlation_inverse(block(0, 0, 2.0 / 7));
  rejects([&] { partial.validate(); }, "incomplete correlated precision");
  rejects([&] { partial.compute_inverse(); }, "partial correlated precision during inversion");

  CovarianceMatrix independent;
  independent.add_correlation(block(0, 0, 4));
  independent.add_correlation(block(1, 1, 2));
  independent.add_correlation_inverse(block(0, 0, 0.25));
  independent.validate();
  require(inverse_blocks(independent).size() == 1, "Validation filled the uncached independent component");
  independent.compute_inverse();
  independent.validate();
  close(independent.get_inverse()[1, 1], 0.5, "Missing independent precision was not computed");

  // Adding a cross correlation makes both old diagonal precisions stale.
  independent.add_correlation({Range(0, 1), Range(1, 1), {0, 1}, matrix(1, 1, {1})});
  require(inverse_blocks(independent).empty(), "Correlation insertion kept stale precision");
  independent.compute_inverse();
  independent.validate();
  close(independent.get_inverse()[0, 1], -1.0 / 7, "Correlation insertion retained old inverse");

  auto edited = covariance(matrix(1, 1, {4}));
  edited.compute_inverse();
  edited.get_blocks()[0].get_dense()[0, 0] = 9;
  require(inverse_blocks(edited).empty(), "Mutable block access retained stale inverse");
  edited.compute_inverse();
  close(edited.get_inverse()[0, 0], 1.0 / 9, "Mutable block update was ignored");
  edited.set_blocks({block(0, 0, 16)});
  require(inverse_blocks(edited).empty(), "Block replacement retained stale inverse");
  edited.compute_inverse();
  close(edited.get_inverse()[0, 0], 1.0 / 16, "Replacement covariance update was ignored");

  // Public C++ matrix aliases can outlive the cache calculation. Validation
  // must catch edits through such aliases even without a mutating accessor.
  auto             shared = std::make_shared<Matrix>(matrix(1, 1, {4}));
  CovarianceMatrix aliased;
  aliased.add_correlation({Range(0, 1), Range(0, 1), {0, 0}, shared});
  aliased.compute_inverse();
  (*shared)[0, 0] = 9;
  require(not inverse_blocks(aliased).empty(), "Const inverse access unexpectedly cleared the cache");
  rejects([&] { aliased.validate(); }, "stale precision through retained matrix alias");
  rejects([&] { aliased.compute_inverse(); }, "stale precision through retained alias during inversion");

  auto malformed = covariance(matrix(1, 1, {4}));
  malformed.add_correlation_inverse(Block{});
  rejects([&] { malformed.validate(); }, "null precision block");
  auto negative = covariance(matrix(1, 1, {4}));
  negative.add_correlation_inverse(block(0, 0, -0.25));
  rejects([&] { negative.validate(); }, "negative precision");
}

void workspace_helpers() {
  JacobianTargets targets;
  targets.atm       = {Jacobian::AtmTarget{.type = SpeciesEnum::Water, .target_pos = 0, .x_start = 0, .x_size = 2}};
  targets.finalized = true;
  CovarianceMatrix prior;
  model_state_covmatAddSpeciesVMR(prior,
                                  targets,
                                  SpeciesEnum::Water,
                                  BlockMatrix{matrix(2, 2, {4, 1, 1, 2})},
                                  BlockMatrix{matrix(2, 2, {2.0 / 7, -1.0 / 7, -1.0 / 7, 4.0 / 7})});
  prior.validate(2);
  require(prior.nblocks() == 1 and inverse_blocks(prior).size() == 1,
          "Target helper duplicated covariance instead of storing inverse");
  close(prior.get_inverse()[0, 1], -1.0 / 7, "Target helper stored wrong supplied inverse");

  ArrayOfSensorObsel sensors(2);
  CovarianceMatrix   error;
  measurement_vec_error_covmatConstant(error, sensors, 0.25);
  error.validate(2);
  close(error.get_inverse()[1, 1], 4, "Constant measurement error precision");
  for (Numeric variance :
       {0.0, -1.0, std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::quiet_NaN()})
    rejects([&] { measurement_vec_error_covmatConstant(error, sensors, variance); },
            "invalid constant measurement variance");
}

// Same-executable comparison with the pre-refactor mult_inv implementation.
}  // namespace

int main(int argc, char** argv) try {
  require(argc == 2, "Specify structural, numerical, precision, or workspace");
  const std::string_view test = argv[1];
  if (test == "structural")
    structural();
  else if (test == "numerical")
    numerical();
  else if (test == "precision")
    precision();
  else if (test == "workspace")
    workspace_helpers();
  else
    throw std::runtime_error("Unknown covariance validation test");
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return EXIT_FAILURE;
}
