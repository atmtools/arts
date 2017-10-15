/* Copyright (C) 2017
   Simon Pfreundschuh <simonpf@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   covariance_matrix.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2017-06-19

  \brief  Implementation of CovarianceMatrix class.
*/

#include <utility>
#include <tuple>
#include <queue>
#include <vector>

#include "covariance_matrix.h"
#include "lapack.h"

//------------------------------------------------------------------------------
// Correlations
//------------------------------------------------------------------------------
void mult(MatrixView C, ConstMatrixView A, const Block &B)
{
    assert((B.sparse_ != nullptr) || (B.dense_ != nullptr)) ;

    MatrixView CView(C(joker, B.get_column_range()));
    MatrixView CTView(C(joker, B.get_row_range()));
    ConstMatrixView AView(A(joker, B.get_row_range()));
    ConstMatrixView ATView(A(joker, B.get_column_range()));

    Index i,j;
    std::tie(i,j) = B.get_indices();

    if (B.matrix_type_ == Block::MatrixType::dense) {
        mult(CView, AView, *B.dense_);
    } else {
        mult(CView, AView, *B.sparse_);
    }

    // Only upper blocks are stored, so if the correlation is between different RQs
    // we also need to account for the implicit lower block.
    if (i != j) {
        if (B.matrix_type_ == Block::MatrixType::dense) {
            mult(CTView, ATView, transpose(*B.dense_));
        } else {
            Matrix D(CTView.ncols(), CTView.nrows());
            mult(D, *B.sparse_, transpose(ATView));
            CTView += transpose(D);
        }
    }
}

void mult(MatrixView C, const Block &A, ConstMatrixView B)
{
    assert((A.sparse_ != nullptr) || (A.dense_ != nullptr)) ;

    MatrixView CView(C(A.get_row_range(), joker));
    MatrixView CTView(C(A.get_column_range(),joker));
    ConstMatrixView BView(B(A.get_column_range(), joker));
    ConstMatrixView BTView(B(A.get_row_range(), joker));

    if (A.matrix_type_ == Block::MatrixType::dense) {
        mult(CView, *A.dense_, BView);
    } else {
        mult(CView, *A.sparse_, BView);
    }

    Index i,j;
    std::tie(i,j) = A.get_indices();

    // Only upper blocks are stored, so if the correlation is between different RQs
    // we also need to account for the implicit lower block.
    if (i != j) {
        if (A.matrix_type_ == Block::MatrixType::dense) {
            mult(CTView, transpose(*A.dense_), BTView);
        } else {
            Matrix D(CTView.ncols(), CTView.nrows());
            mult(D, transpose(BTView), *A.sparse_);
            CTView += transpose(D);
        }
    }
}

void mult(VectorView w, const Block &A, ConstVectorView v)
{
    VectorView wview(w[A.get_row_range()]), wtview(w[A.get_column_range()]);
    ConstVectorView vview(v[A.get_column_range()]), vtview(v[A.get_row_range()]);

    if (A.matrix_type_ == Block::MatrixType::dense) {
        mult(wview, *A.dense_, vview);
    } else {
        mult(wview, *A.sparse_, vview);
    }

    Index i,j;
    std::tie(i,j) = A.get_indices();

    if (i != j) {
        if (A.matrix_type_ == Block::MatrixType::dense) {
            mult(wtview, transpose(*A.dense_), vtview);
        } else {
            transpose_mult(wtview, *A.sparse_, vtview);
        }
    }
}

MatrixView & operator+=(MatrixView &A, const Block &B)
{
    MatrixView Aview(A(B.get_row_range(), B.get_column_range()));
    MatrixView ATview(A(B.get_column_range(), B.get_row_range()));
    if (B.get_matrix_type() == Block::MatrixType::dense) {
        Aview += B.get_dense();
    } else {
        Aview += static_cast<Matrix>(B.get_sparse());
    }

    Index i,j;
    std::tie(i,j) = B.get_indices();

    if (i != j) {
        if (B.get_matrix_type() == Block::MatrixType::dense) {
            ATview += transpose(B.get_dense());
        } else {
            ATview += transpose(static_cast<Matrix>(B.get_sparse()));
        }
    }
    return A;
}

//------------------------------------------------------------------------------
// Covariance Matrix
//------------------------------------------------------------------------------
CovarianceMatrix::operator Matrix() const {

    Index n = nrows();
    Matrix A(n,n);
    A = 0.0;

    for (const Block & c : correlations_) {
        MatrixView Aview = A(c.get_row_range(), c.get_column_range());
        if (c.get_matrix_type() == Block::MatrixType::dense) {
            Aview = c.get_dense();
        } else {
            Aview = c.get_sparse();
        }

        Index ci, cj;
        std::tie(ci, cj) = c.get_indices();
        if (ci != cj) {
            MatrixView ATview = A(c.get_column_range(), c.get_row_range());
            if (c.get_matrix_type() == Block::MatrixType::dense) {
                ATview = transpose(c.get_dense());
            } else {
                ATview = transpose(static_cast<Matrix>(c.get_sparse()));
            }
        }
    }
    return A;
}

Matrix CovarianceMatrix::get_inverse() const {

    Index n = nrows();
    Matrix A(n,n);
    A = 0.0;

    for (const Block & c : inverses_) {
        MatrixView Aview = A(c.get_row_range(), c.get_column_range());
        if (c.get_matrix_type() == Block::MatrixType::dense) {
            Aview = c.get_dense();
        } else {
            Aview = c.get_sparse();
        }

        Index ci, cj;
        std::tie(ci, cj) = c.get_indices();
        if (ci != cj) {
            MatrixView ATview = A(c.get_column_range(), c.get_row_range());
            if (c.get_matrix_type() == Block::MatrixType::dense) {
                ATview = transpose(c.get_dense());
            } else {
                ATview = transpose(static_cast<Matrix>(c.get_sparse()));
            }
        }
    }
    return A;
}

Index CovarianceMatrix::nrows() const
{
    Index m = 0;

    for (const Block & c : correlations_) {
        Index i,j;
        std::tie(i,j) = c.get_indices();
        if (i == j) {
            m += c.nrows();
        }
    }
    return m;
}

Index CovarianceMatrix::ncols() const
{
    return nrows();
}

Index CovarianceMatrix::ndiagblocks() const
{
    Index m = 0;

    for (const Block & c : correlations_) {
        Index i,j;
        std::tie(i,j) = c.get_indices();
        if (i == j) {
            ++m;
        }
    }
    return m;
}

bool CovarianceMatrix::has_diagonal_blocks(const ArrayOfArrayOfIndex &jis) const
{
    for (Index i = 0; i < static_cast<Index>(jis.size()); ++i) {
        Index n_blocks = 0;
        for (const Block &b : correlations_) {
            if (b.get_indices() == std::make_pair(i,i)) {
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
        Index i,j;
        std::tie(i,j) = b.get_indices();

        Index row_start  = jis[i][0];
        Index row_extent = jis[i][1] - jis[i][0] + 1;
        Range row_range  = b.get_row_range();
        if ((row_range.get_start() != row_start)
            || (row_range.get_extent() != row_extent)) {
            return false;
        }

        Index column_start  = jis[i][0];
        Index column_extent = jis[i][1] - jis[i][0] + 1;
        Range column_range  = b.get_column_range();
        if ((column_range.get_start() != column_start)
            || (column_range.get_extent() != column_extent)) {
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

bool CovarianceMatrix::has_inverse(std::pair<Index, Index> indices) const
{
    for (const Block &b : inverses_) {
        if (indices == b.get_indices()) {
            return true;
        }
    }
    return false;
}

void CovarianceMatrix::generate_blocks(std::vector<std::vector<const Block *>> &corr_blocks) const
{
    for (size_t i = 0; i < correlations_.size(); i++) {
        Index ci, cj;
        std::tie(ci, cj) = correlations_[i].get_indices();
    }

    std::vector<bool> has_block(correlations_.size(), false);
    std::queue<Index> rq_queue{};
    for (size_t i = 0; i < correlations_.size(); ++i) {
        if (!has_block[i]) {
            Index ci, cj;
            std::tie(ci, cj) = correlations_[i].get_indices();
            rq_queue.push(ci);
            if (ci != cj) {
                rq_queue.push(cj);
            }
            has_block[i] = true;
            corr_blocks.push_back(std::vector<const Block *>{&correlations_[i]});

            while (!rq_queue.empty()) {
                Index rq_index = rq_queue.front();
                rq_queue.pop();

                for (size_t j = 0; j < correlations_.size(); ++j) {
                    if (!has_block[j]) {
                        std::tie(ci, cj) = correlations_[j].get_indices();
                        if ((ci == rq_index) || (cj == rq_index)) {
                            if (ci != rq_index) {
                                rq_queue.push(ci);
                            }
                            if (cj != rq_index) {
                                rq_queue.push(cj);
                            }
                            corr_blocks.back().push_back(&correlations_[j]);
                            has_block[j] = true;
                        }
                    }
                }
            }
        }
    }
}

void CovarianceMatrix::compute_inverse() const
{
    std::vector<std::vector<const Block*>> correlation_blocks{};
    generate_blocks(correlation_blocks);
    for (std::vector<const Block *> & cb: correlation_blocks) {
        invert_correlation_block(inverses_, cb);
    }
}

void CovarianceMatrix::invert_correlation_block(std::vector<Block> &inverses,
                                                std::vector<const Block *> &blocks) const
{
    // Can't compute inverse of empty block.
    assert(blocks.size() > 0);

    // Sort blocks w.r.t. indices.
    auto comp = [](const Block *a, const Block *b)
        {
            Index a1, a2, b1, b2;
            std::tie(a1, a2) = a->get_indices();
            std::tie(b1, b2) = b->get_indices();
            return ((a1 < b1) || ((a1 == b1) && (a2 < b2)));
        };

    std::sort(blocks.begin(), blocks.end(), comp);


    // Nothing to do if single block which already has inverse.
    if ((blocks.size() == 1) && has_inverse(blocks[0]->get_indices())) {
        return;
    }

    // Otherwise go on to precompute the inverse of a block consisting
    // of correlations between multiple retrieval quantities.

    // The single blocks corresponding to a set of correlated retrieval quantities
    // can be distributed freely over the covariance matrix, so we need to establish
    // a mapping to a continues square matrix to compute the inverse. This is done by
    // mapping the coordinates of each block to a start row and extent in the continuous
    // matrix A. We also record which retrieval quantity indices belong to this
    // set of retrieval quantities.
    Index n = 0;
    std::map<Index,Index> block_start{};
    std::map<Index,Index> block_extent{};
    std::map<Index,Index> block_start_cont{};
    std::map<Index,Index> block_extent_cont{};
    std::vector<Index> block_indices{};

    for (size_t i = 0; i < blocks.size(); ++i) {
        Index ci,cj;
        std::tie(ci,cj) = blocks[i]->get_indices();

        if (ci == cj) {
            Index extent = blocks[i]->get_row_range().get_extent();
            block_start.insert(std::make_pair(ci, blocks[i]->get_row_range().get_start()));
            block_extent.insert(std::make_pair(ci, blocks[i]->get_row_range().get_extent()));
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
        Index ci,cj;
        std::tie(ci,cj) = blocks[i]->get_indices();
        Range row_range(block_start_cont[ci], block_extent_cont[ci]);
        Range column_range(block_start_cont[cj], block_extent_cont[cj]);
        MatrixView A_view = A(row_range, column_range);

        if (blocks[i]->get_matrix_type() == Block::MatrixType::dense) {
            A_view = blocks[i]->get_dense();
        } else {
            A_view = blocks[i]->get_sparse();
        }
    }

    for (Index i = 0; i < n; ++i) {
        for (Index j = i + 1; j < n; ++j) {
            A(j,i) = A(i,j);
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
                MatrixView A_view = A(row_range_A, column_range_A);
                inverses.push_back(Block(row_range, column_range,
                                         std::make_pair(bi, bj),
                                         std::make_shared<Matrix>(A_view)));
            }
        }
    }
}

void CovarianceMatrix::add_correlation(Block c)
{
    correlations_.push_back(c);
}

void CovarianceMatrix::add_correlation_inverse(Block c)
{
    inverses_.push_back(c);
}

template <typename MatrixType>
void CovarianceMatrix::add_correlation(Index i,
                                       Index j,
                                       std::shared_ptr<MatrixType> matrix)
{
    using CorrelationIterator = std::vector<Block>::const_iterator;
    CorrelationIterator i_block = std::find_if(
        correlations_.begin(),
        correlations_.end(),
        [i](const Block &b) {return b.get_indices() == std::make_pair(i,i);}
        );
    CorrelationIterator j_block = std::find_if(
        correlations_.begin(),
        correlations_.end(),
        [j](const Block &b) {return b.get_indices() == std::make_pair(j,j);}
        );

    if ((i_block == correlations_.end()) || (j_block == correlations_.end())) {
            runtime_error("Cannot add correlation between blocks that have not yet been added "
                          "to the covariance matrix.");
    }

    Range row_range(i_block->get_row_range());
    Range column_range(j_block->get_column_range());

    correlations_.emplace_back(row_range, column_range, std::make_pair(i,j), matrix);
}

template <typename MatrixType>
void CovarianceMatrix::add_correlation(Index i,
                                       std::shared_ptr<MatrixType> matrix)
{
    if (i == 0) {
        if (correlations_.size() > 0) {
            runtime_error("Trying to add a correlation block with index 0 to a covariance "
                          "matrix that already contains a correlation. This would lead to "
                          " and inconsistent state and is not supported.");
            }
            Index start = 0;
            Index extent = matrix->nrows();
            Range range(start, extent);
            correlations_.emplace_back(range, range, std::make_pair(i,i), matrix);
    } else {
        using CorrelationIterator = std::vector<Block>::const_iterator;
        auto pred = [i](const Block &b) {return b.get_indices() == std::make_pair(i,i);};

        CorrelationIterator i_block = std::find_if(correlations_.begin(),
                                                   correlations_.end(),
                                                   pred);
        if ((i_block != correlations_.end())) {
            runtime_error("Covariance matrix already contains a block with the same index.");
        }
        auto pred_prev = [i](const Block &b) {return b.get_indices() == std::make_pair(i-1,i-1);};
        CorrelationIterator prev_block = std::find_if(correlations_.begin(),
                                                      correlations_.end(),
                                                      pred_prev);
        if ((prev_block != correlations_.end())) {
            runtime_error("Could not find block preceding the block that should be added.");
        }


        Index start = prev_block->get_row_range().get_start()
                      + prev_block->get_row_range().get_extent();
        Index extent = matrix->nrows();
        Range range(start, extent);
        correlations_.emplace_back(range, range, std::make_pair(i,i), matrix);
    }
}

template <typename MatrixType>
void CovarianceMatrix::add_correlation(Index i,
                                       Index j,
                                       const ArrayOfArrayOfIndex &jacobian_indices,
                                       std::shared_ptr<MatrixType>    matrix)
{
    if (i > j) {
        std::swap(i,j);
    }

    Index row_extent = jacobian_indices[i][1] - jacobian_indices[i][0];
    Index col_extent = jacobian_indices[j][1] - jacobian_indices[j][0];

    Range row_range(jacobian_indices[i][0], row_extent);
    Range column_range(jacobian_indices[j][0], col_extent);


    correlations_.emplace_back(row_range, column_range, std::make_pair(i,j), matrix);
}

void mult(MatrixView C, ConstMatrixView A, const CovarianceMatrix &B)
{
    C = 0.0;
    Matrix T(C);
    for (const Block &c : B.correlations_) {
        T = 0.0;
        mult(T, A, c);
        C += T;
    }
}

void mult(MatrixView C, const CovarianceMatrix &A, ConstMatrixView B)
{
    C = 0.0;
    Matrix T(C);
    for (const Block &c : A.correlations_) {
        T = 0.0;
        mult(T, c, B);
        C += T;
    }
}

void mult(VectorView w, const CovarianceMatrix &A, ConstVectorView v)
{
    w = 0.0;
    Vector t(w);
    for (const Block &c : A.correlations_) {
        t = 0.0;
        mult(t, c, v);
        w += t;
    }
}

void mult_inv(MatrixView C, ConstMatrixView A, const CovarianceMatrix &B)
{
    C = 0.0;
    Matrix T(C);
    for (const Block &c : B.inverses_) {
        T = 0.0;
        mult(T, A, c);
        C += T;
    }
}

void mult_inv(MatrixView C, const CovarianceMatrix &A, ConstMatrixView B)
{
    C = 0.0;
    Matrix T(C);
    for (const Block &c : A.inverses_) {
        T = 0.0;
        mult(T, c, B);
        C += T;
    }
}

void solve(VectorView w, const CovarianceMatrix &A, ConstVectorView v)
{
    w = 0.0;
    Vector t(w);
    for (const Block &c : A.inverses_) {
        t = 0.0;
        mult(t, c, v);
        w += t;
    }
}

MatrixView & operator+=(MatrixView &A, const CovarianceMatrix &B)
{
    for (const Block &c : B.correlations_) {
        A += c;
    }
    return A;
}

void add_inv(MatrixView A, const CovarianceMatrix &B)
{
    for (const Block &c : B.inverses_) {
        A += c;
    }
}

std::ostream& operator<<(std::ostream &os, const CovarianceMatrix &covmat)
{
    os << "Covariance Matrix, ";
    os << "\tDimensions: [" << covmat.nrows() << " x " << covmat.ncols() << "]" << std::endl;
    os << "Blocks:" << std::endl;
    for (const Block &b : covmat.correlations_) {
        Index i,j;
        tie(i,j) = b.get_indices();
        os << "\ti = " << i << ", j = " << j << ": " << b.get_row_range().get_extent();
        os << " x " << b.get_column_range().get_extent();
        os << ", has inverse: " << (covmat.has_inverse(std::make_pair(i,j)) ? "yes" : "no");
        os << std::endl;
    }
    return os;
}

// Explicit instantiations.

template void CovarianceMatrix::add_correlation<Matrix>(Index, std::shared_ptr<Matrix>);
template void CovarianceMatrix::add_correlation<Sparse>(Index, std::shared_ptr<Sparse>);
template void CovarianceMatrix::add_correlation<Matrix>(Index, Index, std::shared_ptr<Matrix>);
template void CovarianceMatrix::add_correlation<Sparse>(Index, Index, std::shared_ptr<Sparse>);
template void CovarianceMatrix::add_correlation<Matrix>(Index, Index,
                                                        const ArrayOfArrayOfIndex &,
                                                        std::shared_ptr<Matrix>);
template void CovarianceMatrix::add_correlation<Sparse>(Index, Index,
                                                        const ArrayOfArrayOfIndex &,
                                                        std::shared_ptr<Sparse>);
