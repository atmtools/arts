#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <map>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "covariance_matrix.h"

namespace {
using Dense = Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic>;

// Iterate only stored entries for sparse matrices. Implicit zero variances are
// checked separately, and never require expanding a diagonal sparse block.
template <typename F> void entries(const Block& block, F&& f) {
  if (block.is_dense()) {
    const auto& matrix = block.get_dense();
    for (Index i = 0; i < matrix.nrows(); ++i)
      for (Index j = 0; j < matrix.ncols(); ++j) f(i, j, matrix[i, j]);
  } else {
    const auto& matrix = block.get_sparse().matrix;
    for (Index i = 0; i < matrix.outerSize(); ++i)
      for (Eigen::SparseMatrix<Numeric, Eigen::RowMajor>::InnerIterator it(matrix, i); it; ++it)
        f(it.row(), it.col(), it.value());
  }
}

Numeric element(const Block& block, Index i, Index j) {
  return block.is_dense() ? block.get_dense()[i, j] : block.get_sparse()[i, j];
}

struct Layout {
  std::map<Index, const Block*>     diagonal;
  std::map<IndexPair, const Block*> blocks;
  Index                             size{};
};

bool same_range(Range a, Range b) { return a.offset == b.offset and a.nelem == b.nelem; }

Layout layout(const std::vector<Block>& blocks, std::string_view name, const Layout* covariance = nullptr) {
  Layout result;
  for (const auto& block : blocks) {
    const auto [i, j] = block.get_indices();
    const auto row    = block.get_row_range();
    const auto col    = block.get_column_range();
    const auto fail   = [&](std::string_view reason) {
      throw std::runtime_error(std::format("{} block ({}, {}): {}.", name, i, j, reason));
    };
    if (i < 0 or j < i) fail("indices must identify an upper triangular block");
    if (not block.not_null()) fail("matrix is null");
    for (const auto range : {row, col}) {
      if (range.offset < 0 or range.nelem <= 0 or range.nelem > std::numeric_limits<Index>::max() - range.offset)
        fail("range must have a nonnegative offset and positive, nonoverflowing extent");
    }
    if (block.matrix_.nrows() != row.nelem or block.matrix_.ncols() != col.nelem)
      fail("matrix shape does not match its row and column ranges");
    if (not result.blocks.emplace(IndexPair{i, j}, &block).second) fail("duplicate block indices");
    if (i == j) {
      if (not same_range(row, col)) fail("diagonal row and column ranges differ");
      result.diagonal.emplace(i, &block);
    }
    entries(block, [&](Index r, Index c, Numeric value) {
      if (not std::isfinite(value))
        throw std::runtime_error(std::format("{} block ({}, {}), coordinate ({}, {}): non-finite value {}.",
                                             name,
                                             i,
                                             j,
                                             row.offset + r,
                                             col.offset + c,
                                             value));
    });
  }
  if (result.diagonal.empty()) throw std::runtime_error(std::format("{} needs diagonal blocks.", name));

  std::vector<std::pair<Index, const Block*>> ordered(result.diagonal.begin(), result.diagonal.end());
  std::ranges::sort(ordered, {}, [](const auto& item) { return item.second->get_row_range().offset; });
  Index end = 0;
  for (const auto& [index, block] : ordered) {
    const auto row = block->get_row_range();
    if (not covariance and row.offset != end)
      throw std::runtime_error(std::format(
          "{} diagonal block {}: diagonal ranges must cover [0, n) without gaps or overlaps; expected offset {}, got {}.",
          name,
          index,
          end,
          row.offset));
    end = row.offset + row.nelem;
  }
  result.size           = end;
  const auto& reference = covariance ? *covariance : result;
  for (const auto& [indices, block] : result.blocks) {
    const auto [i, j] = indices;
    const auto row    = reference.diagonal.find(i);
    const auto col    = reference.diagonal.find(j);
    if (row == reference.diagonal.end() or col == reference.diagonal.end() or
        not same_range(block->get_row_range(), row->second->get_row_range()) or
        not same_range(block->get_column_range(), col->second->get_column_range()))
      throw std::runtime_error(std::format(
          "{} block ({}, {}): ranges must match the corresponding covariance diagonal blocks.", name, i, j));
  }
  return result;
}

void check_diagonals(const Layout& layout, std::string_view name, Numeric tolerance) {
  for (const auto& [index, block] : layout.diagonal) {
    const auto start = block->get_row_range().offset;
    for (Index i = 0; i < block->nrows(); ++i) {
      const auto value = element(*block, i, i);
      if (value <= 0)
        throw std::runtime_error(
            std::format("{} diagonal block {}, coordinate ({}, {}): variance must be strictly positive, got {}.",
                        name,
                        index,
                        start + i,
                        start + i,
                        value));
    }
    entries(*block, [&](Index i, Index j, Numeric value) {
      if (i == j) return;
      const long double scale = std::sqrt(static_cast<long double>(element(*block, i, i))) *
                                std::sqrt(static_cast<long double>(element(*block, j, j)));
      if (std::abs(static_cast<long double>(value) - element(*block, j, i)) / scale > tolerance)
        throw std::runtime_error(std::format(
            "{} diagonal block {}, coordinate ({}, {}): matrix is not symmetric within relative_tolerance {} after variance scaling.",
            name,
            index,
            start + i,
            start + j,
            tolerance));
    });
  }
}

struct Component {
  std::vector<Index>        indices;
  std::vector<const Block*> covariance;
  std::vector<const Block*> precision;
};

std::vector<Component> components(const Layout& covariance, const Layout& precision) {
  std::map<Index, Index> parents;
  for (const auto& [index, block] : covariance.diagonal) parents[index] = index;
  const auto root = [&](Index i) {
    while (parents[i] != i) i = parents[i];
    return i;
  };
  // Precision edges are included so unexpected correlations between independent
  // covariance components are checked against the complete inverse relation.
  for (const auto* layout : {&covariance, &precision}) {
    for (const auto& [indices, block] : layout->blocks) {
      const auto [i, j] = indices;
      parents[root(i)]  = root(j);
    }
  }
  std::map<Index, Component> groups;
  for (const auto& [index, block] : covariance.diagonal) groups[root(index)].indices.push_back(index);
  for (const auto& [indices, block] : covariance.blocks) groups[root(indices.first)].covariance.push_back(block);
  for (const auto& [indices, block] : precision.blocks) groups[root(indices.first)].precision.push_back(block);
  std::vector<Component> result;
  for (auto& [index, component] : groups) result.push_back(std::move(component));
  return result;
}

bool diagonal_only(const std::vector<const Block*>& blocks) {
  bool diagonal = true;
  for (const auto* block : blocks) {
    const auto [i, j] = block->get_indices();
    entries(*block, [&](Index r, Index c, Numeric value) {
      if (value != 0 and (i != j or r != c)) diagonal = false;
    });
  }
  return diagonal;
}

void check_component(const Component& component,
                     const Layout&    covariance,
                     const Layout&    precision,
                     Numeric          tolerance,
                     Index            max_dense_elements) {
  const auto             first         = component.indices.front();
  const bool             has_precision = not component.precision.empty();
  std::map<Index, Index> starts;
  Index                  size = 0;
  for (const auto index : component.indices) {
    starts[index]  = size;
    size          += covariance.diagonal.at(index)->nrows();
    if (has_precision and not precision.diagonal.contains(index))
      throw std::runtime_error(std::format(
          "Inverse covariance is incomplete: missing diagonal block {} in the component starting at block {}.",
          index,
          first));
  }
  if (diagonal_only(component.covariance) and diagonal_only(component.precision)) {
    if (has_precision) {
      for (const auto index : component.indices) {
        const auto& c = *covariance.diagonal.at(index);
        const auto& p = *precision.diagonal.at(index);
        for (Index i = 0; i < c.nrows(); ++i) {
          const long double product = static_cast<long double>(element(c, i, i)) * element(p, i, i);
          if (not std::isfinite(product) or std::abs(product - 1) > tolerance)
            throw std::runtime_error(std::format(
                "Inverse covariance block {}, coordinate {}: covariance times inverse is {}, expected 1 within relative_tolerance {}.",
                index,
                c.get_row_range().offset + i,
                product,
                tolerance));
        }
      }
    }
    return;
  }

  if (size > max_dense_elements / size)
    throw std::runtime_error(std::format(
        "Covariance component starting at block {} needs a {} by {} dense factorization, exceeding max_dense_elements {}. Increase the explicit validation limit only if sufficient memory is available.",
        first,
        size,
        size,
        max_dense_elements));

  // Work in correlation units so physically small or large variances do not
  // create an arbitrary absolute positive-definiteness threshold.
  std::vector<long double> scales(size);
  for (const auto index : component.indices) {
    const auto& block = *covariance.diagonal.at(index);
    for (Index i = 0; i < block.nrows(); ++i)
      scales[starts[index] + i] = std::sqrt(static_cast<long double>(element(block, i, i)));
  }
  const auto assemble = [&](const std::vector<const Block*>& blocks, bool inverse) {
    Dense matrix = Dense::Zero(size, size);
    for (const auto* block : blocks) {
      const auto [i, j] = block->get_indices();
      entries(*block, [&](Index r, Index c, Numeric value) {
        const auto row        = starts[i] + r;
        const auto col        = starts[j] + c;
        const auto scale      = scales[row] * scales[col];
        const auto normalized = static_cast<Numeric>(inverse ? value * scale : value / scale);
        if (not std::isfinite(normalized))
          throw std::runtime_error(std::format(
              "{} covariance component starting at block {}: scaled coordinate ({}, {}) exceeds numerical range.",
              inverse ? "Inverse" : "",
              first,
              row,
              col));
        matrix(row, col) = normalized;
        if (i != j) matrix(col, row) = normalized;
      });
    }
    return matrix;
  };
  const Dense correlation = assemble(component.covariance, false);
  if (Eigen::LLT<Dense>(correlation).info() != Eigen::Success)
    throw std::runtime_error(std::format(
        "Covariance component starting at diagonal block {} is not positive definite after variance scaling (Cholesky factorization failed).",
        first));
  if (has_precision) {
    const Dense inverse = assemble(component.precision, true);
    if (Eigen::LLT<Dense>(inverse).info() != Eigen::Success)
      throw std::runtime_error(std::format(
          "Inverse covariance component starting at diagonal block {} is not positive definite (Cholesky factorization failed).",
          first));
    const Dense   residual = correlation * inverse - Dense::Identity(size, size);
    const Numeric error    = residual.cwiseAbs().maxCoeff();
    const Numeric limit    = tolerance * static_cast<Numeric>(size);
    if (not std::isfinite(error) or error > limit)
      throw std::runtime_error(std::format(
          "Inverse covariance component starting at block {} is inconsistent: maximum variance-scaled identity residual {} exceeds {}. The stored inverse may be stale or omit correlations.",
          first,
          error,
          limit));
  }
}
}  // namespace

void CovarianceMatrix::validate(Index expected_size, Numeric relative_tolerance, Index max_dense_elements) const {
  if (expected_size < -1) throw std::runtime_error("expected_size must be -1 or nonnegative.");
  if (max_dense_elements <= 0) throw std::runtime_error("max_dense_elements must be positive.");
  if (not std::isfinite(relative_tolerance) or relative_tolerance <= 0 or relative_tolerance >= 1)
    throw std::runtime_error("relative_tolerance must be finite and strictly between 0 and 1.");
  const auto covariance = layout(correlations_, "Covariance");
  if (expected_size >= 0 and covariance.size != expected_size)
    throw std::runtime_error(
        std::format("Covariance size {} does not match expected_size {}.", covariance.size, expected_size));
  check_diagonals(covariance, "Covariance", relative_tolerance);
  Layout precision;
  if (not inverses_.empty()) {
    precision = layout(inverses_, "Inverse covariance", &covariance);
    check_diagonals(precision, "Inverse covariance", relative_tolerance);
  }
  for (const auto& component : components(covariance, precision))
    check_component(component, covariance, precision, relative_tolerance, max_dense_elements);
}
