/*!
  \file   covmat.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2017-06-19

  \brief  Implementation of CovarianceMatrix class.
*/

#include "covariance_matrix.h"

#include <lin_alg.h>
#include <xml.h>

#include <Eigen/Cholesky>
#include <cmath>
#include <limits>
#include <mutex>
#include <optional>
#include <ostream>
#include <queue>
#include <set>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace {
struct CovarianceSignature {
  std::vector<Index>   layout;
  std::vector<Numeric> values;
  bool                 operator==(const CovarianceSignature &) const = default;
};
// One canonical traversal, shared by signature construction and comparison, so
// the two can never disagree about what identifies a covariance. A sink
// returning false stops the walk.
template <typename LayoutSink, typename ValueSink>
bool walk_covariance(const CovarianceMatrix &covariance, LayoutSink layout, ValueSink value) {
  const auto emit = [&](std::initializer_list<Index> items) {
    for (Index item : items)
      if (not layout(item)) return false;
    return true;
  };
  for (const auto *blocks : {&covariance.get_blocks(), &covariance.get_inverse_blocks()}) {
    if (not emit({static_cast<Index>(blocks->size())})) return false;
    for (const auto &b : *blocks) {
      if (not b.not_null()) throw std::runtime_error("Cannot prepare a null covariance block.");
      const auto [i, j] = b.get_indices();
      if (not emit({b.is_dense() ? b.get_dense().nrows() : b.get_sparse().nrows(),
                    b.is_dense() ? b.get_dense().ncols() : b.get_sparse().ncols(),
                    i,
                    j,
                    b.get_row_range().offset,
                    b.nrows(),
                    b.get_column_range().offset,
                    b.ncols(),
                    static_cast<Index>(b.is_dense())}))
        return false;
      if (b.is_dense()) {
        const auto &a = b.get_dense();
        for (auto it = a.elem_begin(); it != a.elem_end(); ++it)
          if (not value(*it)) return false;
      } else {
        const auto &a = b.get_sparse();
        if (not emit({static_cast<Index>(a.nnz())})) return false;
        for (const auto [row, col, entry] : a | by_elem) {
          if (not emit({row, col})) return false;
          if (not value(entry)) return false;
        }
      }
    }
  }
  return true;
}

CovarianceSignature covariance_signature(const CovarianceMatrix &covariance) {
  CovarianceSignature result;
  auto &[layout, values] = result;
  walk_covariance(
      covariance,
      [&layout](Index item) {
        layout.push_back(item);
        return true;
      },
      [&values](Numeric entry) {
        values.push_back(entry);
        return true;
      });
  return result;
}

// Compare against a stored signature without materializing a second one. The
// steady-state answer is "unchanged", and that case now costs no allocation and
// no copy; a change stops the walk at the first differing element. Detection is
// unchanged: every element is still read, so edits made through retained
// references or shared block storage are still caught.
bool signature_matches(const CovarianceMatrix     &covariance,
                       const std::vector<Index>   &layout,
                       const std::vector<Numeric> &values) {
  std::size_t layout_at = 0, value_at = 0;
  const bool  complete = walk_covariance(
      covariance,
      [&](Index item) { return layout_at < layout.size() and layout[layout_at++] == item; },
      [&](Numeric entry) { return value_at < values.size() and values[value_at++] == entry; });
  return complete and layout_at == layout.size() and value_at == values.size();
}

bool signature_matches(const CovarianceMatrix &covariance, const CovarianceSignature &expected) {
  return signature_matches(covariance, expected.layout, expected.values);
}
}  // namespace
struct CovariancePreparation {
  std::mutex                              mutex;
  std::vector<Index>                      layout;
  std::vector<Numeric>                    values;
  std::optional<CovarianceSignature>      validated_inverse;
  bool                                    precision = false;
  std::shared_ptr<const CovarianceMatrix> snapshot;
};

namespace {
std::vector<Block> detached_blocks(const std::vector<Block> &source) {
  std::vector<Block> out;
  out.reserve(source.size());
  for (const auto &b : source) {
    if (b.is_dense())
      out.emplace_back(b.get_row_range(), b.get_column_range(), b.get_indices(), b.get_dense());
    else
      out.emplace_back(b.get_row_range(), b.get_column_range(), b.get_indices(), b.get_sparse());
  }
  return out;
}
}  // namespace

BlockMatrix::BlockMatrix()                                   = default;
BlockMatrix::BlockMatrix(const BlockMatrix &)                = default;
BlockMatrix::BlockMatrix(BlockMatrix &&) noexcept            = default;
BlockMatrix &BlockMatrix::operator=(const BlockMatrix &)     = default;
BlockMatrix &BlockMatrix::operator=(BlockMatrix &&) noexcept = default;

Block::Block()                             = default;
Block::Block(const Block &)                = default;
Block::Block(Block &&) noexcept            = default;
Block &Block::operator=(const Block &)     = default;
Block &Block::operator=(Block &&) noexcept = default;
Block::~Block()                            = default;
CovarianceMatrix::CovarianceMatrix() : preparation_(std::make_shared<CovariancePreparation>()) {}
CovarianceMatrix::CovarianceMatrix(const CovarianceMatrix &other)
    : preparation_(std::make_shared<CovariancePreparation>()) {
  std::lock_guard lock(other.preparation_->mutex);
  solve_cache_  = other.solve_cache_;
  correlations_ = other.finalized_ ? detached_blocks(other.correlations_) : other.correlations_;
  inverses_     = other.finalized_ ? detached_blocks(other.inverses_) : other.inverses_;
}
CovarianceMatrix::CovarianceMatrix(CovarianceMatrix &&other) noexcept
    : solve_cache_(std::move(other.solve_cache_)),
      // Retain a valid preparation object in the empty source without allocating
      // in a noexcept move. Exact signatures distinguish its subsequent edits.
      preparation_(other.preparation_),
      finalized_(std::exchange(other.finalized_, false)),
      correlations_(std::exchange(other.correlations_, {})),
      inverses_(std::exchange(other.inverses_, {})) {}
CovarianceMatrix &CovarianceMatrix::operator=(const CovarianceMatrix &other) {
  if (this != &other) {
    CovarianceMatrix copy(other);
    *this = std::move(copy);
  }
  return *this;
}
CovarianceMatrix &CovarianceMatrix::operator=(CovarianceMatrix &&other) noexcept {
  if (this != &other) {
    solve_cache_  = std::move(other.solve_cache_);
    preparation_  = other.preparation_;
    finalized_    = std::exchange(other.finalized_, false);
    correlations_ = std::exchange(other.correlations_, {});
    inverses_     = std::exchange(other.inverses_, {});
  }
  return *this;
}
CovarianceMatrix::~CovarianceMatrix() = default;

void CovarianceMatrix::clear_cache() {
  auto fresh = std::make_shared<CovariancePreparation>();
  solve_cache_.reset();
  preparation_ = std::move(fresh);
  finalized_   = false;
  inverses_    = std::vector<Block>{};
}

void CovarianceMatrix::validate(Index expected_size, Numeric relative_tolerance, Index max_dense_elements) const {
  std::lock_guard lock(preparation_->mutex);
  validate_unlocked(expected_size, relative_tolerance, max_dense_elements);
}

Block::Block(Range row_range, Range column_range, IndexPair indices, BlockMatrix matrix)
    : row_range_(row_range), column_range_(column_range), indices_(std::move(indices)), matrix_(std::move(matrix)) {
  // Nothing to do here.
}

BlockMatrix::BlockMatrix(std::shared_ptr<Matrix> dense) : data(std::move(dense)) {}

BlockMatrix::BlockMatrix(std::shared_ptr<Sparse> sparse) : data(std::move(sparse)) {}

BlockMatrix::BlockMatrix(const Matrix &dense) : data(std::make_shared<Matrix>(dense)) {}

BlockMatrix::BlockMatrix(const Sparse &sparse) : data(std::make_shared<Sparse>(sparse)) {}

BlockMatrix &BlockMatrix::operator=(std::shared_ptr<Matrix> dense) {
  data = std::move(dense);
  return *this;
}

BlockMatrix &BlockMatrix::operator=(std::shared_ptr<Sparse> sparse) {
  data = std::move(sparse);
  return *this;
}

BlockMatrix &BlockMatrix::operator=(const Matrix &dense) {
  data = std::make_shared<Matrix>(dense);
  return *this;
}

BlockMatrix &BlockMatrix::operator=(const Sparse &sparse) {
  data = std::make_shared<Sparse>(sparse);
  return *this;
}

bool BlockMatrix::not_null() const {
  return std::visit([]<typename T>(const std::shared_ptr<T> &matrix) { return matrix != nullptr; }, data);
}

bool BlockMatrix::is_dense() const { return std::holds_alternative<std::shared_ptr<Matrix>>(data); }

bool BlockMatrix::is_sparse() const { return not is_dense(); }

Matrix &BlockMatrix::dense() {
  assert(is_dense());
  return *std::get<std::shared_ptr<Matrix>>(data);
}

const Matrix &BlockMatrix::dense() const {
  assert(is_dense());
  return *std::get<std::shared_ptr<Matrix>>(data);
}

Sparse &BlockMatrix::sparse() {
  assert(is_sparse());
  return *std::get<std::shared_ptr<Sparse>>(data);
}

const Sparse &BlockMatrix::sparse() const {
  assert(is_sparse());
  return *std::get<std::shared_ptr<Sparse>>(data);
}

Vector BlockMatrix::diagonal() const {
  return std::visit(
      []<typename T>(const std::shared_ptr<T> &matrix) -> Vector {
        if constexpr (std::same_as<T, Sparse>)
          return matrix->diagonal();
        else
          return Vector{matpack::diagonal(*matrix)};
      },
      data);
}

Index BlockMatrix::ncols() const {
  return std::visit([]<typename T>(const std::shared_ptr<T> &matrix) { return matrix ? matrix->ncols() : 0; }, data);
}

Index BlockMatrix::nrows() const {
  return std::visit([]<typename T>(const std::shared_ptr<T> &matrix) { return matrix ? matrix->nrows() : 0; }, data);
}

void BlockMatrix::multiply_left(StridedMatrixView out, StridedConstMatrixView rhs) const {
  std::visit([&]<typename T>(const std::shared_ptr<T> &matrix) { mult(out, *matrix, rhs); }, data);
}

void BlockMatrix::multiply_left(StridedVectorView out, StridedConstVectorView rhs) const {
  std::visit([&]<typename T>(const std::shared_ptr<T> &matrix) { mult(out, *matrix, rhs); }, data);
}

void BlockMatrix::multiply_right(StridedMatrixView out, StridedConstMatrixView lhs) const {
  std::visit([&]<typename T>(const std::shared_ptr<T> &matrix) { mult(out, lhs, *matrix); }, data);
}

bool BlockMatrix::is_finite() const {
  return std::visit(
      []<typename T>(const std::shared_ptr<T> &matrix) {
        if (not matrix) return false;
        if constexpr (std::same_as<T, Sparse>) {
          for (const auto [row, col, value] : *matrix | by_elem)
            if (not std::isfinite(value)) return false;
        } else {
          for (const auto value : *matrix | by_elem)
            if (not std::isfinite(value)) return false;
        }
        return true;
      },
      data);
}

bool BlockMatrix::is_identity() const {
  if (not not_null() or nrows() == 0 or nrows() != ncols()) return false;
  return std::visit(
      []<typename T>(const std::shared_ptr<T> &matrix) {
        if constexpr (std::same_as<T, Sparse>) {
          Index diagonal_count = 0;
          for (const auto [row, col, value] : *matrix | by_elem) {
            if (row == col) {
              if (value != 1) return false;
              ++diagonal_count;
            } else if (value != 0)
              return false;
          }
          return diagonal_count == matrix->nrows();
        } else {
          for (Index row = 0; row < matrix->nrows(); ++row)
            for (Index col = 0; col < matrix->ncols(); ++col)
              if ((*matrix)[row, col] != (row == col ? 1 : 0)) return false;
          return true;
        }
      },
      data);
}

void Block::set_matrix(std::shared_ptr<Sparse> sparse) { matrix_ = std::move(sparse); }
void Block::set_matrix(std::shared_ptr<Matrix> dense) { matrix_ = std::move(dense); }

std::array<Index, 2> BlockMatrix::shape() const {
  return std::visit(
      []<typename T>(const std::shared_ptr<T> &matrix) -> std::array<Index, 2> {
        return {matrix->nrows(), matrix->ncols()};
      },
      data);
}

//------------------------------------------------------------------------------
// Correlations
//------------------------------------------------------------------------------
void mult(StridedMatrixView C, StridedConstMatrixView A, const Block &B) {
  assert(B.not_null());

  StridedMatrixView      CView(C[joker, B.get_column_range()]);
  StridedMatrixView      CTView(C[joker, B.get_row_range()]);
  StridedConstMatrixView AView(A[joker, B.get_row_range()]);
  StridedConstMatrixView ATView(A[joker, B.get_column_range()]);

  Index i, j;
  std::tie(i, j) = B.get_indices();

  if (B.is_dense()) {
    mult(CView, AView, B.get_dense());
  } else {
    mult(CView, AView, B.get_sparse());
  }

  // Only upper blocks are stored, so if the correlation is between different RQs
  // we also need to account for the implicit lower block.
  if (i != j) {
    if (B.is_dense()) {
      mult(CTView, ATView, transpose(B.get_dense()));
    } else {
      Matrix D(CTView.ncols(), CTView.nrows());
      mult(D, B.get_sparse(), transpose(ATView));
      CTView += transpose(D);
    }
  }
}

void mult(StridedMatrixView C, const Block &A, StridedConstMatrixView B) {
  assert(A.not_null());

  StridedMatrixView      CView(C[A.get_row_range(), joker]);
  StridedMatrixView      CTView(C[A.get_column_range(), joker]);
  StridedConstMatrixView BView(B[A.get_column_range(), joker]);
  StridedConstMatrixView BTView(B[A.get_row_range(), joker]);

  if (A.is_dense()) {
    mult(CView, A.get_dense(), BView);
  } else {
    mult(CView, A.get_sparse(), BView);
  }

  Index i, j;
  std::tie(i, j) = A.get_indices();

  // Only upper blocks are stored, so if the correlation is between different RQs
  // we also need to account for the implicit lower block.
  if (i != j) {
    if (A.is_dense()) {
      mult(CTView, transpose(A.get_dense()), BTView);
    } else {
      Matrix D(CTView.ncols(), CTView.nrows());
      mult(D, transpose(BTView), A.get_sparse());
      CTView += transpose(D);
    }
  }
}

void mult(StridedVectorView w, const Block &A, StridedConstVectorView v) {
  StridedVectorView      wview(w[A.get_row_range()]), wtview(w[A.get_column_range()]);
  StridedConstVectorView vview(v[A.get_column_range()]), vtview(v[A.get_row_range()]);

  if (A.is_dense()) {
    mult(wview, A.get_dense(), vview);
  } else {
    mult(wview, A.get_sparse(), vview);
  }

  Index i, j;
  std::tie(i, j) = A.get_indices();

  if (i != j) {
    if (A.is_dense()) {
      mult(wtview, transpose(A.get_dense()), vtview);
    } else {
      transpose_mult(wtview, A.get_sparse(), vtview);
    }
  }
}

StridedMatrixView operator+=(StridedMatrixView A, const Block &B) {
  StridedMatrixView Aview(A[B.get_row_range(), B.get_column_range()]);
  StridedMatrixView ATview(A[B.get_column_range(), B.get_row_range()]);
  if (B.is_dense()) {
    Aview += B.get_dense();
  } else {
    Aview += static_cast<const Matrix>(B.get_sparse());
  }

  Index i, j;
  std::tie(i, j) = B.get_indices();

  if (i != j) {
    if (B.is_dense()) {
      ATview += transpose(B.get_dense());
    } else {
      ATview += transpose(static_cast<const Matrix>(B.get_sparse()));
    }
  }
  return A;
}

//------------------------------------------------------------------------------
// Covariance Matrix
//------------------------------------------------------------------------------
CovarianceMatrix::operator Matrix() const {
  Index  n = nrows();
  Matrix A(n, n);
  A = 0.0;

  for (const Block &c : correlations_) {
    StridedMatrixView Aview = A[c.get_row_range(), c.get_column_range()];
    if (c.is_dense()) {
      Aview = c.get_dense();
    } else {
      Aview = static_cast<const Matrix>(c.get_sparse());
    }

    Index ci, cj;
    std::tie(ci, cj) = c.get_indices();
    if (ci != cj) {
      StridedMatrixView ATview = A[c.get_column_range(), c.get_row_range()];
      if (c.is_dense()) {
        ATview = transpose(c.get_dense());
      } else {
        ATview = transpose(static_cast<const Matrix>(c.get_sparse()));
      }
    }
  }
  return A;
}

Matrix CovarianceMatrix::get_inverse() const {
  compute_inverse();
  Index  n = nrows();
  Matrix A(n, n);
  A = 0.0;

  for (const Block &c : inverses_) {
    StridedMatrixView Aview = A[c.get_row_range(), c.get_column_range()];
    if (c.is_dense()) {
      Aview = c.get_dense();
    } else {
      Aview = static_cast<const Matrix>(c.get_sparse());
    }

    Index ci, cj;
    std::tie(ci, cj) = c.get_indices();
    if (ci != cj) {
      StridedMatrixView ATview = A[c.get_column_range(), c.get_row_range()];
      if (c.is_dense()) {
        ATview = transpose(c.get_dense());
      } else {
        ATview = transpose(static_cast<const Matrix>(c.get_sparse()));
      }
    }
  }
  return A;
}

Index CovarianceMatrix::nrows() const {
  Index m1 = 0;

  for (const Block &c : correlations_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) { m1 += c.nrows(); }
  }

  Index m2 = 0;
  for (const Block &c : inverses_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) { m2 += c.nrows(); }
  }

  return std::max(m1, m2);
}

Index CovarianceMatrix::ncols() const { return nrows(); }

Index CovarianceMatrix::ndiagblocks() const {
  Index m = 0;

  for (const Block &c : correlations_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) { ++m; }
  }
  return m;
}

Index CovarianceMatrix::ninvdiagblocks() const {
  Index m = 0;

  for (const Block &c : inverses_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) { ++m; }
  }
  return m;
}

Index CovarianceMatrix::nblocks() const { return correlations_.size(); }

bool CovarianceMatrix::has_block(Index i, Index j) {
  if (i > j) { std::swap(i, j); }

  bool result = false;
  for (const Block &b : correlations_) { result |= b.get_indices() == std::make_pair(i, j); }

  return result;
}

const Block *CovarianceMatrix::get_block(Index i, Index j) {
  if (i > j) { std::swap(i, j); }
  Index bi, bj;
  for (const Block &b : correlations_) {
    std::tie(bi, bj) = b.get_indices();
    if (((i == bi) && (j == bj)) || ((i == -1) && (j == bj)) || ((i == -1) && (j == -1))) { return &b; }
  }
  return nullptr;
}

bool CovarianceMatrix::has_diagonal_blocks(const ArrayOfArrayOfIndex &jis) const {
  for (Index i = 0; i < static_cast<Index>(jis.size()); ++i) {
    Index n_blocks = 0;
    for (const Block &b : correlations_) {
      if (b.get_indices() == std::make_pair(i, i)) { ++n_blocks; }
    }
    if (n_blocks != 1) { return false; }
  }
  return true;
}

bool CovarianceMatrix::is_consistent(const ArrayOfArrayOfIndex &jis) const {
  auto pred = [&jis](const Block &b) {
    Index i, j;
    std::tie(i, j) = b.get_indices();

    Index row_start  = jis[i][0];
    Index row_extent = jis[i][1] - jis[i][0] + 1;
    Range row_range  = b.get_row_range();
    if ((row_range.offset != row_start) || (row_range.nelem != row_extent)) { return false; }

    Index column_start  = jis[j][0];
    Index column_extent = jis[j][1] - jis[j][0] + 1;
    Range column_range  = b.get_column_range();
    if ((column_range.offset != column_start) || (column_range.nelem != column_extent)) { return false; }
    return true;
  };

  if (!std::all_of(correlations_.begin(), correlations_.end(), pred)) { return false; }
  if (!std::all_of(inverses_.begin(), inverses_.end(), pred)) { return false; }
  return true;
}

bool CovarianceMatrix::is_consistent(const Block &b) const {
  Index i, j;
  std::tie(i, j) = b.get_indices();

  for (const Block &c : correlations_) {
    Index ii, jj;
    std::tie(ii, jj) = c.get_indices();

    if ((ii == i) && (c.nrows() != b.nrows())) { return false; }

    if ((jj == j) && (c.ncols() != b.ncols())) { return false; }

    if ((ii == i) && (jj == j)) { return false; }
  }
  return true;
}

bool CovarianceMatrix::has_inverse(IndexPair indices) const {
  for (const Block &b : inverses_) {
    if (indices == b.get_indices()) { return true; }
  }
  return false;
}

void CovarianceMatrix::generate_blocks(std::vector<std::vector<const Block *>> &corr_blocks) const {
  for (size_t i = 0; i < correlations_.size(); i++) {
    Index ci, cj;
    std::tie(ci, cj) = correlations_[i].get_indices();
  }

  std::vector<bool> has_blocks(correlations_.size(), false);
  std::queue<Index> rq_queue{};
  for (size_t i = 0; i < correlations_.size(); ++i) {
    if (!has_blocks[i]) {
      Index ci, cj;
      std::tie(ci, cj) = correlations_[i].get_indices();
      rq_queue.push(ci);
      if (ci != cj) { rq_queue.push(cj); }
      has_blocks[i] = true;
      corr_blocks.push_back(std::vector<const Block *>{&correlations_[i]});

      while (!rq_queue.empty()) {
        Index rq_index = rq_queue.front();
        rq_queue.pop();

        for (size_t j = 0; j < correlations_.size(); ++j) {
          if (!has_blocks[j]) {
            std::tie(ci, cj) = correlations_[j].get_indices();
            if ((ci == rq_index) || (cj == rq_index)) {
              if (ci != rq_index) { rq_queue.push(ci); }
              if (cj != rq_index) { rq_queue.push(cj); }
              corr_blocks.back().push_back(&correlations_[j]);
              has_blocks[j] = true;
            }
          }
        }
      }
    }
  }
}

void CovarianceMatrix::compute_inverse() const {
  if (finalized_) {
    if (inverses_.empty()) throw std::runtime_error("Prepared covariance did not request explicit precision.");
    return;
  }
  // Inverse-only matrices are also used internally as precision operators.
  if (correlations_.empty()) return;
  std::lock_guard lock(preparation_->mutex);
  // Reuse the complete, previously validated inverse. An exact signature also
  // catches edits through retained block references and shared matrix storage.
  if (preparation_->validated_inverse and signature_matches(*this, *preparation_->validated_inverse)) return;
  // Independent components may have supplied inverses while others still need
  // computing. A represented component must already have a complete inverse.
  // The optional confirmation API has a user-configurable allocation guard.
  // Existing inversion callers already request a dense component inverse and
  // must not inherit an unconfigurable size limit from that separate API.
  validate_unlocked(-1, 1e-10, std::numeric_limits<Index>::max());
  std::vector<std::vector<const Block *>> correlation_blocks{};
  generate_blocks(correlation_blocks);
  for (std::vector<const Block *> &cb : correlation_blocks) { invert_correlation_block(inverses_, cb); }
  // Remember completed validation, including the resulting inverse values.
  // Preparation still checks exact storage, so aliases cannot bypass validation.
  preparation_->validated_inverse = covariance_signature(*this);
}

void CovarianceMatrix::invert_correlation_block(std::vector<Block>         &inverses,
                                                std::vector<const Block *> &blocks) const {
  // Can't compute inverse of empty block.
  assert(blocks.size() > 0);

  // Sort blocks w.r.t. indices.
  auto comp = [](const Block *a, const Block *b) {
    Index a1, a2, b1, b2;
    std::tie(a1, a2) = a->get_indices();
    std::tie(b1, b2) = b->get_indices();
    return ((a1 < b1) || ((a1 == b1) && (a2 < b2)));
  };

  std::sort(blocks.begin(), blocks.end(), comp);

  // validate has checked any supplied/cached inverse for this complete
  // component. Its sparsity pattern need not equal that of the covariance.
  if (std::any_of(blocks.begin(), blocks.end(), [this](const Block *a) {
        const auto [i, j] = a->get_indices();
        return i == j and has_inverse({i, i});
      }))
    return;

  // The usual independent measurement covariance should remain sparse.
  if (blocks.size() == 1 and blocks.front()->is_sparse()) {
    const Block &block    = *blocks.front();
    const auto  &sparse   = block.get_sparse();
    bool         diagonal = true;
    for (const auto [row, col, value] : sparse | by_elem)
      if (row != col and value != 0) diagonal = false;
    if (diagonal) {
      Vector values = block.diagonal();
      for (auto &value : values) {
        value = 1 / value;
        if (not std::isfinite(value)) throw std::runtime_error("Covariance inverse diagonal exceeds numerical range.");
      }
      inverses.emplace_back(
          block.get_row_range(), block.get_column_range(), block.get_indices(), Sparse::diagonal(values));
      return;
    }
  }

  // Otherwise go on to precompute the inverse of a block consisting
  // of correlations between multiple retrieval quantities.

  // The single blocks corresponding to a set of correlated retrieval quantities
  // can be distributed freely over the covariance matrix, so we need to establish
  // a mapping to a continuous square matrix to compute the inverse. This is done by
  // mapping the coordinates of each block to a start row and extent in the continuous
  // matrix A. We also record which retrieval quantity indices belong to this
  // set of retrieval quantities.
  Index                  n = 0;
  std::map<Index, Index> block_start{};
  std::map<Index, Index> block_extent{};
  std::map<Index, Index> block_start_cont{};
  std::map<Index, Index> block_extent_cont{};
  std::vector<Index>     block_indices{};

  for (size_t i = 0; i < blocks.size(); ++i) {
    Index ci, cj;
    std::tie(ci, cj) = blocks[i]->get_indices();

    if (ci == cj) {
      Index extent = blocks[i]->get_row_range().nelem;
      block_start.insert(std::make_pair(ci, blocks[i]->get_row_range().offset));
      block_extent.insert(std::make_pair(ci, blocks[i]->get_row_range().nelem));
      block_start_cont.insert(std::make_pair(ci, n));
      block_extent_cont.insert(std::make_pair(ci, extent));
      block_indices.push_back(ci);
      n += extent;
    }
  }

  // Copy blocks into a single dense matrix.
  Matrix A(n, n);
  A = 0.0;

  for (size_t i = 0; i < blocks.size(); ++i) {
    Index ci, cj;
    std::tie(ci, cj) = blocks[i]->get_indices();
    Range             row_range(block_start_cont[ci], block_extent_cont[ci]);
    Range             column_range(block_start_cont[cj], block_extent_cont[cj]);
    StridedMatrixView A_view = A[row_range, column_range];

    if (blocks[i]->is_dense()) {
      A_view = blocks[i]->get_dense();
    } else {
      A_view = static_cast<const Matrix>(blocks[i]->get_sparse());
    }
    // Only off-diagonal blocks have an implicit transpose. A diagonal block
    // stores both triangles, already checked for symmetry before this step.
    if (ci != cj) A[column_range, row_range] = transpose(A_view);
  }
  inv(A, A);

  if (std::any_of(A.elem_begin(), A.elem_end(), [](Numeric value) { return not std::isfinite(value); }))
    throw std::runtime_error(
        "Covariance inverse contains non-finite values; the component cannot be inverted in the available numerical range.");

  // // Invert matrix using LAPACK.
  // char uplo = 'L';
  // int  ni, info1(0), info2(0);
  // ni = static_cast<int>(n);
  // lapack::dpotrf_(&uplo, &ni, A.get_raw_data(), &ni, &info1);
  // lapack::dpotri_(&uplo, &ni, A.get_raw_data(), &ni, &info2);
  // if ((info1 != 0) || info2 !=0) {
  //     throw std::runtime_error("Error inverting block of covariance matrix."
  //                              "Make sure that it is symmetric, positive definite"
  //                              "or provide the inverse manually.");
  // }

  // Now we need to disassemble the matrix inverse bach to the separate block in the
  // covariance matrix. Note, however, that blocks that previously were implicitly
  // zero are now non-zero, i.e. the inverse may contain more blocks than the covariance
  // matrix itself.
  for (Index bi : block_indices) {
    for (Index bj : block_indices) {
      if (bi <= bj) {
        Range             row_range_A(block_start_cont[bi], block_extent_cont[bi]);
        Range             column_range_A(block_start_cont[bj], block_extent_cont[bj]);
        Range             row_range(block_start[bi], block_extent[bi]);
        Range             column_range(block_start[bj], block_extent[bj]);
        StridedMatrixView A_view = A[row_range_A, column_range_A];
        inverses.push_back(Block(row_range, column_range, std::make_pair(bi, bj), std::make_shared<Matrix>(A_view)));
      }
    }
  }
}

void CovarianceMatrix::set_blocks(std::vector<Block> blocks) {
  correlations_ = std::move(blocks);
  inverses_.clear();
}

void CovarianceMatrix::add_correlation(Block c) {
  // A new edge can join previously independent covariance components. Their
  // old inverses no longer apply, but inverses of untouched components do.
  const auto [i, j] = c.get_indices();
  std::set<Index> touched{i, j};
  Size            previous;
  do {
    previous = touched.size();
    for (const auto &block : correlations_) {
      const auto [row, col] = block.get_indices();
      if (touched.contains(row) or touched.contains(col)) {
        touched.insert(row);
        touched.insert(col);
      }
    }
  } while (previous != touched.size());
  std::erase_if(inverses_, [&](const Block &block) {
    const auto [row, col] = block.get_indices();
    return touched.contains(row) or touched.contains(col);
  });
  correlations_.push_back(std::move(c));
}

void CovarianceMatrix::add_correlation_inverse(Block c) { inverses_.push_back(std::move(c)); }

Vector CovarianceMatrix::diagonal() const {
  Vector diag(nrows());
  for (const Block &b : correlations_) {
    Index i, j;
    std::tie(i, j) = b.get_indices();

    if (i == j) { diag[b.get_row_range()] = b.diagonal(); }
  }
  return diag;
}

Vector CovarianceMatrix::inverse_diagonal() const {
  compute_inverse();

  Vector diag(nrows());
  for (const Block &b : inverses_) {
    Index i, j;
    std::tie(i, j) = b.get_indices();

    if (i == j) { diag[b.get_row_range()] = b.diagonal(); }
  }
  return diag;
}

namespace {
// Discover exact diagonal structure without materializing a dense covariance.
// No persistent cache: block storage can be shared with external callers.
std::optional<Vector> diagonal_values(const std::vector<Block> &blocks, Index n) {
  if (blocks.empty()) return std::nullopt;
  Vector values(n, 0.);
  for (const Block &block : blocks) {
    const auto [i, j] = block.get_indices();
    if (i != j) return std::nullopt;
    if (block.is_dense()) {
      const auto &a = block.get_dense();
      for (Index r = 0; r < a.nrows(); ++r)
        for (Index c = 0; c < a.ncols(); ++c)
          if (r != c and a[r, c] != 0) return std::nullopt;
    } else {
      const auto &a = block.get_sparse();
      for (const auto [row, col, value] : a | by_elem)
        if (row != col and value != 0) return std::nullopt;
    }
    const Vector d      = block.diagonal();
    const Index  offset = block.get_row_range().offset;
    // A well-formed covariance has its diagonal blocks tiling [0, n), so n
    // covers every row. Block storage is settable without validation, so a
    // gapped layout can make n too small; such a layout is not a diagonal of
    // length n, and answering "not diagonal" leaves it to the generic path.
    if (offset < 0 or offset + static_cast<Index>(d.size()) > n) return std::nullopt;
    for (Index r = 0; r < static_cast<Index>(d.size()); ++r) values[offset + r] = d[r];
  }
  return values;
}
}  // namespace

// Structural types are internal: public inputs remain Matrix and Sparse.
std::optional<Vector> CovarianceMatrix::diagonal_if_diagonal() const { return diagonal_values(correlations_, nrows()); }

struct CovarianceSolveCache {
  struct Diagonal {
    Vector values;
  };
  struct Cholesky {
    Eigen::LLT<Eigen::MatrixXd> factor;
  };
  struct Component {
    std::vector<Index>               rows;
    std::variant<Diagonal, Cholesky> solver;
  };
  std::vector<Index>     layout;
  std::vector<Numeric>   values;
  std::vector<Component> components;
};

std::shared_ptr<const CovarianceMatrix> CovarianceMatrix::prepared(bool need_precision) const {
  std::lock_guard lock(preparation_->mutex);
  auto           &cache = *preparation_;
  // Answer the common "nothing changed" case before building a signature.
  if (cache.snapshot and cache.precision == need_precision and signature_matches(*this, cache.layout, cache.values))
    return cache.snapshot;
  auto signature         = covariance_signature(*this);
  auto &[layout, values] = signature;
  if (correlations_.empty() and inverses_.empty()) throw std::runtime_error("Cannot prepare an empty covariance.");
  auto snapshot           = std::make_shared<CovarianceMatrix>();
  snapshot->correlations_ = detached_blocks(correlations_);
  snapshot->inverses_     = detached_blocks(inverses_);
  if (need_precision or not inverses_.empty()) {
    if (not cache.validated_inverse or not signature_matches(*this, *cache.validated_inverse))
      snapshot->compute_inverse();
  } else {
    Matrix empty(0, 0);
    snapshot->solve_components(empty, empty);
  }
  snapshot->finalized_ = true;
  cache.layout         = std::move(layout);
  cache.values         = std::move(values);
  cache.precision      = need_precision;
  cache.snapshot       = std::move(snapshot);
  return cache.snapshot;
}

bool CovarianceMatrix::solve_components(StridedMatrixView out, StridedConstMatrixView rhs) const {
  // Only source cache construction/publication needs synchronization. A
  // published snapshot has immutable factors, while independent solves retain
  // their own cache handle and release the lock before doing any arithmetic.
  std::unique_lock lock(preparation_->mutex, std::defer_lock);
  if (not finalized_) lock.lock();
  // Preserve explicitly supplied (including precision-only) representations.
  if (correlations_.empty() or not inverses_.empty()) return false;
  if (not finalized_) {
    // Exact snapshots detect mutations through retained references/shared storage.
    if (not solve_cache_ or not signature_matches(*this, solve_cache_->layout, solve_cache_->values)) {
      auto [layout, values] = covariance_signature(*this);
      validate_unlocked(-1, 1e-10, std::numeric_limits<Index>::max());
      auto cache    = std::make_shared<CovarianceSolveCache>();
      cache->layout = std::move(layout);
      cache->values = std::move(values);
      std::vector<std::vector<const Block *>> groups;
      generate_blocks(groups);
      for (const auto &group : groups) {
        CovarianceSolveCache::Component component;
        std::map<Index, Index>          starts;
        for (const auto *b : group) {
          const auto [i, j] = b->get_indices();
          if (i != j) continue;
          starts[i] = component.rows.size();
          for (Index k = 0; k < b->nrows(); ++k) component.rows.push_back(b->get_row_range().offset + k);
        }
        if (group.size() == 1) {
          // Reuse the exact structure classifier; never infer independence from small values.
          if (auto d = diagonal_values(std::vector<Block>{*group.front()}, nrows())) {
            Vector local(component.rows.size());
            for (Index i = 0; i < static_cast<Index>(local.size()); ++i) local[i] = (*d)[component.rows[i]];
            component.solver = CovarianceSolveCache::Diagonal{std::move(local)};
            cache->components.push_back(std::move(component));
            continue;
          }
        }
        const Index     n     = component.rows.size();
        Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(n, n);
        for (const auto *b : group) {
          const auto [i, j] = b->get_indices();
          const Index r0 = starts.at(i), c0 = starts.at(j);
          auto        put = [&](Index r, Index c, Numeric v) {
            dense(r0 + r, c0 + c) = v;
            if (i != j) dense(c0 + c, r0 + r) = v;
          };
          if (b->is_dense()) {
            for (Index r = 0; r < b->nrows(); ++r)
              for (Index c = 0; c < b->ncols(); ++c) put(r, c, b->get_dense()[r, c]);
          } else {
            const auto &a = b->get_sparse();
            for (const auto [row, col, value] : a | by_elem) put(row, col, value);
          }
        }
        CovarianceSolveCache::Cholesky factor{Eigen::LLT<Eigen::MatrixXd>(dense)};
        if (factor.factor.info() != Eigen::Success)
          throw std::runtime_error("Covariance component Cholesky factorization failed.");
        component.solver = std::move(factor);
        cache->components.push_back(std::move(component));
      }
      solve_cache_ = std::move(cache);
    }
  }
  const auto                                 *cache = solve_cache_.get();
  std::shared_ptr<const CovarianceSolveCache> retained_cache;
  if (lock.owns_lock()) {
    retained_cache = solve_cache_;
    lock.unlock();
  }
  for (const auto &component : cache->components) {
    std::visit(
        [&]<typename T>(const T &solver) {
          if constexpr (std::same_as<T, CovarianceSolveCache::Diagonal>) {
            for (Index i = 0; i < static_cast<Index>(component.rows.size()); ++i)
              for (Index j = 0; j < rhs.ncols(); ++j)
                out[component.rows[i], j] = rhs[component.rows[i], j] / solver.values[i];
          } else {
            Eigen::MatrixXd local(component.rows.size(), rhs.ncols());
            for (Index i = 0; i < local.rows(); ++i)
              for (Index j = 0; j < local.cols(); ++j) local(i, j) = rhs[component.rows[i], j];
            Eigen::MatrixXd solved = solver.factor.solve(local);
            for (Index i = 0; i < local.rows(); ++i)
              for (Index j = 0; j < local.cols(); ++j) out[component.rows[i], j] = solved(i, j);
          }
        },
        component.solver);
  }
  return true;
}

CovarianceSquareRoot::CovarianceSquareRoot(const CovarianceMatrix &covariance) {
  covariance.validate(covariance.nrows(), 1e-10, std::numeric_limits<Index>::max());
  // Work on detached covariance blocks: supplied inverses do not define L.
  CovarianceMatrix local;
  local.correlations_ = detached_blocks(covariance.get_blocks());
  Matrix empty(0, 0);
  local.solve_components(empty, empty);
  factors_ = std::move(local.solve_cache_);
}

void CovarianceSquareRoot::multiply_left(StridedMatrixView out, StridedConstMatrixView rhs, bool transpose) const {
  apply(out, rhs, false, transpose);
}

void CovarianceSquareRoot::solve_left(StridedMatrixView out, StridedConstMatrixView rhs, bool transpose) const {
  apply(out, rhs, true, transpose);
}

void CovarianceSquareRoot::apply(StridedMatrixView      out,
                                 StridedConstMatrixView rhs,
                                 bool                   inverse,
                                 bool                   transpose) const {
  assert(out.shape() == rhs.shape());
  for (const auto &component : factors_->components) {
    std::visit(
        [&]<typename T>(const T &solver) {
          if constexpr (std::same_as<T, CovarianceSolveCache::Diagonal>) {
            for (Index i = 0; i < static_cast<Index>(component.rows.size()); ++i) {
              const Numeric scale = std::sqrt(solver.values[i]);
              for (Index j = 0; j < rhs.ncols(); ++j)
                out[component.rows[i], j] =
                    inverse ? rhs[component.rows[i], j] / scale : rhs[component.rows[i], j] * scale;
            }
          } else {
            Eigen::MatrixXd local(component.rows.size(), rhs.ncols());
            for (Index i = 0; i < local.rows(); ++i)
              for (Index j = 0; j < local.cols(); ++j) local(i, j) = rhs[component.rows[i], j];
            Eigen::MatrixXd result;
            if (inverse) {
              if (transpose)
                result = solver.factor.matrixU().solve(local);
              else
                result = solver.factor.matrixL().solve(local);
            } else {
              if (transpose)
                result = solver.factor.matrixU() * local;
              else
                result = solver.factor.matrixL() * local;
            }
            for (Index i = 0; i < local.rows(); ++i)
              for (Index j = 0; j < local.cols(); ++j) out[component.rows[i], j] = result(i, j);
          }
        },
        component.solver);
  }
}

void mult(StridedMatrixView C, StridedConstMatrixView A, const CovarianceMatrix &B) {
  if (auto d = diagonal_values(B.correlations_, B.nrows())) {
    for (Index i = 0; i < C.nrows(); ++i)
      for (Index j = 0; j < C.ncols(); ++j) C[i, j] = A[i, j] * (*d)[j];
    return;
  }
  C = 0.0;
  Matrix T(C);
  for (const Block &c : B.correlations_) {
    T = 0.0;
    mult(T, A, c);
    C += T;
  }
}

void mult(StridedMatrixView C, const CovarianceMatrix &A, StridedConstMatrixView B) {
  if (auto d = diagonal_values(A.correlations_, A.nrows())) {
    for (Index i = 0; i < C.nrows(); ++i)
      for (Index j = 0; j < C.ncols(); ++j) C[i, j] = (*d)[i] * B[i, j];
    return;
  }
  C = 0.0;
  Matrix T(C);
  for (const Block &c : A.correlations_) {
    T = 0.0;
    mult(T, c, B);
    C += T;
  }
}

void mult(StridedVectorView w, const CovarianceMatrix &A, StridedConstVectorView v) {
  w = 0.0;
  Vector t(w);
  for (const Block &c : A.correlations_) {
    t = 0.0;
    mult(t, c, v);
    w += t;
  }
}

void mult_inv(StridedMatrixView C, StridedConstMatrixView A, const CovarianceMatrix &B) {
  if (B.solve_components(transpose(C), transpose(A))) return;
  B.compute_inverse();
  C = 0.0;
  Matrix T(C);
  for (const Block &c : B.inverses_) {
    T = 0.0;
    mult(T, A, c);
    C += T;
  }
}

void mult_inv(StridedMatrixView C, const CovarianceMatrix &A, StridedConstMatrixView B) {
  if (A.solve_components(C, B)) return;
  A.compute_inverse();
  C = 0.0;
  Matrix T(C);
  for (const Block &c : A.inverses_) {
    T = 0.0;
    mult(T, c, B);
    C += T;
  }
}

void solve(StridedVectorView w, const CovarianceMatrix &A, StridedConstVectorView v) {
  // Views preserve arbitrary vector strides.
  Matrix input(v.size(), 1), output(v.size(), 1);
  for (Index i = 0; i < static_cast<Index>(v.size()); ++i) input[i, 0] = v[i];
  if (A.solve_components(output, input)) {
    for (Index i = 0; i < static_cast<Index>(w.size()); ++i) w[i] = output[i, 0];
    return;
  }
  A.compute_inverse();
  w = 0.0;
  Vector t(w);
  for (const Block &c : A.inverses_) {
    t = 0.0;
    mult(t, c, v);
    w += t;
  }
}

StridedMatrixView operator+=(StridedMatrixView A, const CovarianceMatrix &B) {
  for (const Block &c : B.correlations_) { A += c; }
  return A;
}

void add_inv(StridedMatrixView A, const CovarianceMatrix &B) {
  B.compute_inverse();
  for (const Block &c : B.inverses_) { A += c; }
}

std::ostream &operator<<(std::ostream &os, const CovarianceMatrix &covmat) {
  os << "Covariance Matrix, ";
  os << "\tDimensions: [" << covmat.nrows() << " x " << covmat.ncols() << "]" << '\n';
  os << "Blocks:" << '\n';
  for (const Block &b : covmat.correlations_) {
    Index i, j;
    std::tie(i, j) = b.get_indices();
    os << "\ti = " << i << ", j = " << j << ": " << b.get_row_range().nelem;
    os << " x " << b.get_column_range().nelem;
    os << ", has inverse: " << (covmat.has_inverse(std::make_pair(i, j)) ? "yes" : "no");
    os << '\n';
  }
  return os;
}

void xml_io_stream<Block>::write(std::ostream &os, const Block &x, bofstream *pbofs, std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.row_range_, pbofs);
  xml_write_to_stream(os, x.column_range_, pbofs);
  xml_write_to_stream(os, x.indices_, pbofs);
  xml_write_to_stream(os, x.matrix_, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Block>::read(std::istream &is, Block &x, bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.row_range_, pbifs);
  xml_read_from_stream(is, x.column_range_, pbifs);
  xml_read_from_stream(is, x.indices_, pbifs);
  xml_read_from_stream(is, x.matrix_, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<BlockMatrix>::write(std::ostream      &os,
                                       const BlockMatrix &x,
                                       bofstream         *pbofs,
                                       std::string_view   name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.data, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<BlockMatrix>::read(std::istream &is, BlockMatrix &x, bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.data, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<CovarianceMatrix>::write(std::ostream           &os,
                                            const CovarianceMatrix &x,
                                            bofstream              *pbofs,
                                            std::string_view        name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.get_blocks(), pbofs);
  xml_write_to_stream(os, x.get_inverse_blocks(), pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<CovarianceMatrix>::read(std::istream &is, CovarianceMatrix &x, bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.get_blocks(), pbifs);
  xml_read_from_stream(is, x.get_inverse_blocks(), pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
