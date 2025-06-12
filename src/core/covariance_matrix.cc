/*!
  \file   covariance_matrix.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2017-06-19

  \brief  Implementation of CovarianceMatrix class.
*/

#include "covariance_matrix.h"

#include <lin_alg.h>

#include <ostream>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

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
  if (is_dense()) return std::get<std::shared_ptr<Matrix>>(data) != nullptr;
  return std::get<std::shared_ptr<Sparse>>(data) != nullptr;
}

bool BlockMatrix::is_dense() const {
  return std::holds_alternative<std::shared_ptr<Matrix>>(data);
}

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
  if (is_dense())
    return Vector{matpack::diagonal(*std::get<std::shared_ptr<Matrix>>(data))};
  return std::get<std::shared_ptr<Sparse>>(data)->diagonal();
}

Index BlockMatrix::ncols() const {
  if (is_dense()) return dense().ncols();
  return sparse().ncols();
}

Index BlockMatrix::nrows() const {
  if (is_dense()) return dense().nrows();
  return sparse().nrows();
}

std::ostream &operator<<(std::ostream &os, const BlockMatrix &m) {
  if (m.is_dense()) {
    os << m.dense();
  } else {
    os << m.sparse();
  }
  return os;
}

std::array<Index, 2> BlockMatrix::shape() const {
  if (is_dense()) {
    return dense().shape();
  }
  return {sparse().nrows(), sparse().ncols()};
}

//------------------------------------------------------------------------------
// Correlations
//------------------------------------------------------------------------------
void mult(StridedMatrixView C, StridedConstMatrixView A, const Block &B) {
  assert(B.not_null());

  StridedMatrixView CView(C[joker, B.get_column_range()]);
  StridedMatrixView CTView(C[joker, B.get_row_range()]);
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

  StridedMatrixView CView(C[A.get_row_range(), joker]);
  StridedMatrixView CTView(C[A.get_column_range(), joker]);
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
  StridedVectorView wview(w[A.get_row_range()]),
      wtview(w[A.get_column_range()]);
  StridedConstVectorView vview(v[A.get_column_range()]),
      vtview(v[A.get_row_range()]);

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
  Index n = nrows();
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
  Index n = nrows();
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
    if (i == j) {
      m1 += c.nrows();
    }
  }

  Index m2 = 0;
  for (const Block &c : inverses_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) {
      m2 += c.nrows();
    }
  }

  return std::max(m1, m2);
}

Index CovarianceMatrix::ncols() const { return nrows(); }

Index CovarianceMatrix::ndiagblocks() const {
  Index m = 0;

  for (const Block &c : correlations_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) {
      ++m;
    }
  }
  return m;
}

Index CovarianceMatrix::ninvdiagblocks() const {
  Index m = 0;

  for (const Block &c : inverses_) {
    Index i, j;
    std::tie(i, j) = c.get_indices();
    if (i == j) {
      ++m;
    }
  }
  return m;
}

Index CovarianceMatrix::nblocks() const { return correlations_.size(); }

bool CovarianceMatrix::has_block(Index i, Index j) {
  if (i > j) {
    std::swap(i, j);
  }

  bool result = false;
  for (const Block &b : correlations_) {
    result |= b.get_indices() == std::make_pair(i, j);
  }

  return result;
}

const Block *CovarianceMatrix::get_block(Index i, Index j) {
  if (i > j) {
    std::swap(i, j);
  }
  Index bi, bj;
  for (const Block &b : correlations_) {
    std::tie(bi, bj) = b.get_indices();
    if (((i == bi) && (j == bj)) || ((i == -1) && (j == bj)) ||
        ((i == -1) && (j == -1))) {
      return &b;
    }
  }
  return nullptr;
}

bool CovarianceMatrix::has_diagonal_blocks(
    const ArrayOfArrayOfIndex &jis) const {
  for (Index i = 0; i < static_cast<Index>(jis.size()); ++i) {
    Index n_blocks = 0;
    for (const Block &b : correlations_) {
      if (b.get_indices() == std::make_pair(i, i)) {
        ++n_blocks;
      }
    }
    if (n_blocks != 1) {
      return false;
    }
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
    if ((row_range.offset != row_start) || (row_range.nelem != row_extent)) {
      return false;
    }

    Index column_start  = jis[j][0];
    Index column_extent = jis[j][1] - jis[j][0] + 1;
    Range column_range  = b.get_column_range();
    if ((column_range.offset != column_start) ||
        (column_range.nelem != column_extent)) {
      return false;
    }
    return true;
  };

  if (!std::all_of(correlations_.begin(), correlations_.end(), pred)) {
    return false;
  }
  if (!std::all_of(inverses_.begin(), inverses_.end(), pred)) {
    return false;
  }
  return true;
}

bool CovarianceMatrix::is_consistent(const Block &b) const {
  Index i, j;
  std::tie(i, j) = b.get_indices();

  for (const Block &c : correlations_) {
    Index ii, jj;
    std::tie(ii, jj) = c.get_indices();

    if ((ii == i) && (c.nrows() != b.nrows())) {
      return false;
    }

    if ((jj == j) && (c.ncols() != b.ncols())) {
      return false;
    }

    if ((ii == i) && (jj == j)) {
      return false;
    }
  }
  return true;
}

bool CovarianceMatrix::has_inverse(IndexPair indices) const {
  for (const Block &b : inverses_) {
    if (indices == b.get_indices()) {
      return true;
    }
  }
  return false;
}

void CovarianceMatrix::generate_blocks(
    std::vector<std::vector<const Block *>> &corr_blocks) const {
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
      if (ci != cj) {
        rq_queue.push(cj);
      }
      has_blocks[i] = true;
      corr_blocks.push_back(std::vector<const Block *>{&correlations_[i]});

      while (!rq_queue.empty()) {
        Index rq_index = rq_queue.front();
        rq_queue.pop();

        for (size_t j = 0; j < correlations_.size(); ++j) {
          if (!has_blocks[j]) {
            std::tie(ci, cj) = correlations_[j].get_indices();
            if ((ci == rq_index) || (cj == rq_index)) {
              if (ci != rq_index) {
                rq_queue.push(ci);
              }
              if (cj != rq_index) {
                rq_queue.push(cj);
              }
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
  std::vector<std::vector<const Block *>> correlation_blocks{};
  generate_blocks(correlation_blocks);
  for (std::vector<const Block *> &cb : correlation_blocks) {
    invert_correlation_block(inverses_, cb);
  }
}

void CovarianceMatrix::invert_correlation_block(
    std::vector<Block> &inverses, std::vector<const Block *> &blocks) const {
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

  auto block_has_inverse = [this](const Block *a) {
    return has_inverse(a->get_indices());
  };
  if (std::all_of(blocks.begin(), blocks.end(), block_has_inverse)) return;

  // Otherwise go on to precompute the inverse of a block consisting
  // of correlations between multiple retrieval quantities.

  // The single blocks corresponding to a set of correlated retrieval quantities
  // can be distributed freely over the covariance matrix, so we need to establish
  // a mapping to a continuous square matrix to compute the inverse. This is done by
  // mapping the coordinates of each block to a start row and extent in the continuous
  // matrix A. We also record which retrieval quantity indices belong to this
  // set of retrieval quantities.
  Index n = 0;
  std::map<Index, Index> block_start{};
  std::map<Index, Index> block_extent{};
  std::map<Index, Index> block_start_cont{};
  std::map<Index, Index> block_extent_cont{};
  std::vector<Index> block_indices{};

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
    Range row_range(block_start_cont[ci], block_extent_cont[ci]);
    Range column_range(block_start_cont[cj], block_extent_cont[cj]);
    StridedMatrixView A_view = A[row_range, column_range];

    if (blocks[i]->is_dense()) {
      A_view = blocks[i]->get_dense();
    } else {
      A_view = static_cast<const Matrix>(blocks[i]->get_sparse());
    }
  }

  for (Index i = 0; i < n; ++i) {
    for (Index j = i + 1; j < n; ++j) {
      A[j, i] = A[i, j];
    }
  }
  inv(A, A);

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
        Range row_range_A(block_start_cont[bi], block_extent_cont[bi]);
        Range column_range_A(block_start_cont[bj], block_extent_cont[bj]);
        Range row_range(block_start[bi], block_extent[bi]);
        Range column_range(block_start[bj], block_extent[bj]);
        StridedMatrixView A_view = A[row_range_A, column_range_A];
        inverses.push_back(Block(row_range,
                                 column_range,
                                 std::make_pair(bi, bj),
                                 std::make_shared<Matrix>(A_view)));
      }
    }
  }
}

void CovarianceMatrix::add_correlation(Block c) { correlations_.push_back(c); }

void CovarianceMatrix::add_correlation_inverse(Block c) {
  inverses_.push_back(c);
}

Vector CovarianceMatrix::diagonal() const {
  Vector diag(nrows());
  for (const Block &b : correlations_) {
    Index i, j;
    std::tie(i, j) = b.get_indices();

    if (i == j) {
      diag[b.get_row_range()] = b.diagonal();
    }
  }
  return diag;
}

Vector CovarianceMatrix::inverse_diagonal() const {
  compute_inverse();

  Vector diag(nrows());
  for (const Block &b : inverses_) {
    Index i, j;
    std::tie(i, j) = b.get_indices();

    if (i == j) {
      diag[b.get_row_range()] = b.diagonal();
    }
  }
  return diag;
}

void mult(StridedMatrixView C,
          StridedConstMatrixView A,
          const CovarianceMatrix &B) {
  C = 0.0;
  Matrix T(C);
  for (const Block &c : B.correlations_) {
    T = 0.0;
    mult(T, A, c);
    C += T;
  }
}

void mult(StridedMatrixView C,
          const CovarianceMatrix &A,
          StridedConstMatrixView B) {
  C = 0.0;
  Matrix T(C);
  for (const Block &c : A.correlations_) {
    T = 0.0;
    mult(T, c, B);
    C += T;
  }
}

void mult(StridedVectorView w,
          const CovarianceMatrix &A,
          StridedConstVectorView v) {
  w = 0.0;
  Vector t(w);
  for (const Block &c : A.correlations_) {
    t = 0.0;
    mult(t, c, v);
    w += t;
  }
}

void mult_inv(StridedMatrixView C,
              StridedConstMatrixView A,
              const CovarianceMatrix &B) {
  C = 0.0;
  Matrix T(C);
  for (const Block &c : B.inverses_) {
    T = 0.0;
    mult(T, A, c);
    C += T;
  }
}

void mult_inv(StridedMatrixView C,
              const CovarianceMatrix &A,
              StridedConstMatrixView B) {
  C = 0.0;
  Matrix T(C);
  for (const Block &c : A.inverses_) {
    T = 0.0;
    mult(T, c, B);
    C += T;
  }
}

void solve(StridedVectorView w,
           const CovarianceMatrix &A,
           StridedConstVectorView v) {
  w = 0.0;
  Vector t(w);
  for (const Block &c : A.inverses_) {
    t = 0.0;
    mult(t, c, v);
    w += t;
  }
}

StridedMatrixView operator+=(StridedMatrixView A, const CovarianceMatrix &B) {
  for (const Block &c : B.correlations_) {
    A += c;
  }
  return A;
}

void add_inv(StridedMatrixView A, const CovarianceMatrix &B) {
  for (const Block &c : B.inverses_) {
    A += c;
  }
}

std::ostream &operator<<(std::ostream &os, const CovarianceMatrix &covmat) {
  os << "Covariance Matrix, ";
  os << "\tDimensions: [" << covmat.nrows() << " x " << covmat.ncols() << "]"
     << '\n';
  os << "Blocks:" << '\n';
  for (const Block &b : covmat.correlations_) {
    Index i, j;
    std::tie(i, j) = b.get_indices();
    os << "\ti = " << i << ", j = " << j << ": " << b.get_row_range().nelem;
    os << " x " << b.get_column_range().nelem;
    os << ", has inverse: "
       << (covmat.has_inverse(std::make_pair(i, j)) ? "yes" : "no");
    os << '\n';
  }
  return os;
}
