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
  \file   covariance_matrix.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2017-06-19

  \brief  Header files of CovarianceMatrix class.

  Notes:

  Defines the CovarianceMatrix class which implements the specific structure of
  covariance matrices and their inverse.

*/

#ifndef covariance_matrix_h
#define covariance_matrix_h

#include <memory>

#include "jacobian.h"
#include "matpackI.h"
#include "matpackII.h"

class CovarianceMatrix;

//------------------------------------------------------------------------------
// Type Aliases
//------------------------------------------------------------------------------

using IndexPair = std::pair<Index, Index>;

//------------------------------------------------------------------------------
// Block
//------------------------------------------------------------------------------
/*! A block in a covariance matrix.
 *
 * A block in a covariance matrix represents a correlation between two retrieval
 * quantities. Each block holds a pointer to either a dense matrix of type
 * Matrix or type Sparse. In addition to this it holds to two block
 * indices i and j and two range objects that  describe the position of the block
 * in terms of the rows and columns of blocks in the matrix and the row and
 * columns of elements, respectively.
*/
class Block {
public:

    /*
    * Enum class representing the type of the matrix representing the block, i.e.
    * whether its of type Matrix (dense) or of type Sparse(sparse).
    */
    enum class MatrixType {dense, sparse};

    /*
     * Create a correlation block from given row_range, column_range and shared_ptr
     * to a dense matrix.
     *
     * @param row_range The range of element-rows covered by the block
     * @param column_range The range of element-column covered by the block
     * @param indices The pair of block-indices that identifies the block in the
     *        covariance matrix.
     * @param dense A shared pointer to a den matrix of type Matrix
     */
    Block(Range row_range, Range column_range, IndexPair indices,
                std::shared_ptr<Matrix> dense)
        : row_range_(row_range), column_range_(column_range), indices_(indices),
        matrix_type_(MatrixType::dense), dense_(dense), sparse_(nullptr)
    {
        // Nothing to do here.
    }

    /*
     * Create a correlation block from given row_range, column_range and shared_ptr
     * to a sparse matrix.
     *
     * @param row_range The range of element-rows covered by the block
     * @param column_range The range of element-column covered by the block
     * @param indices The pair of block-indices that identifies the block in the
     *        covariance matrix.
     * @param dense A shared pointer to a den matrix of type Sparse
     */
    Block(Range row_range, Range column_range, IndexPair indices,
                std::shared_ptr<Sparse> sparse)
        : row_range_(row_range), column_range_(column_range), indices_(indices),
          matrix_type_(MatrixType::sparse), dense_(nullptr), sparse_(sparse)
    {
        // Nothing to do here.
    }

    Block(const Block &)  = default;
    Block(      Block &&) = default;
    Block & operator=(const Block &)  = default;
    Block & operator=(      Block &&) = default;

    ~Block() = default;

    /*! Number of rows of this block */
    Index nrows() const {return row_range_.get_extent();}
    /*! Number of columns of this block */
    Index ncols() const {return column_range_.get_extent();}

    /*! The row range of this block*/
    Range get_row_range() const {return row_range_;}
    /*! The column range of this block*/
    Range get_column_range() const {return column_range_;}

    void set_matrix(std::shared_ptr<Sparse> sparse)     {sparse_ = sparse;}
    void set_matrix(std::shared_ptr<Matrix> dense)      {dense_  = dense;}

    /*! Return the diagonal as a vector.*/
    Vector diagonal() const
    {
        if (dense_) {
            return dense_->diagonal();
        } else {
            return sparse_->diagonal();
        }
    }

    /*! Return the indices of the retrieval quantities correlated by this block as std::pair. */
    IndexPair  get_indices()         const {return indices_;}
    /*! Return the type of the matrix holding the correlation coefficients. */
    MatrixType get_matrix_type()     const {return matrix_type_;}

    const Matrix & get_dense()  const
    {
        assert(dense_);
        return *dense_;
    }

    Matrix & get_dense()
    {
        assert(dense_);
        return *dense_;
    }

    const Sparse & get_sparse() const
    {
        assert(sparse_);
        return *sparse_;
    }

    Sparse & get_sparse()
    {
        assert(sparse_);
        return *sparse_;
    }

    // Friend declarations.
    friend void mult(MatrixView, ConstMatrixView, const Block &);
    friend void mult(MatrixView, const Block &, ConstMatrixView);
    friend void mult(VectorView, const Block &, ConstVectorView);

    friend MatrixView & operator+=(MatrixView &, const Block &);

private:

    Range     row_range_, column_range_;
    IndexPair indices_;

    MatrixType matrix_type_;

    std::shared_ptr<Matrix> dense_;
    std::shared_ptr<Sparse> sparse_;
};

void mult(MatrixView, ConstMatrixView, const Block &);
void mult(MatrixView, const Block &, ConstMatrixView);
void mult(VectorView, const Block &, ConstVectorView);

MatrixView & operator+=(MatrixView &, const Block &);
void add_inv(MatrixView A, const Block &);

//------------------------------------------------------------------------------
// Covariance Matrices
//------------------------------------------------------------------------------
/*! A Covariance Matrix
 *
 * This class represents correlations between retrieval quantities in the form
 * of a covariance matrix. The covariance matrix is represented by a vector of
 * blocks that describe the correlations between the retrieval quantities.
 *
 * A block in a covariance matrix is identified by its block-row and
 * block-column indices. The block-row i is defined by the row of blocks that
 * contains the ith diagonal block. Similarly, the block-column j is defined
 * as  the column of block  that contains the jth diagonal block.
 *
 * The covariance matrix class implements overloads for the common
 * ARTS operations for multiplication by matrices and vectors as well as
 * adding covariance matrices to other matrices.
 *
 * Computing inverses of covariance matrices is handled indirectly by providing
 * mult_inv methods that multiply the inverse of the covariance matrix by a given
 * vector or matrix. This, however, requires previously having computed the inverse
 * of the matrix using the compute_inverse method.
 */
class CovarianceMatrix {

public:

    CovarianceMatrix() = default;
    CovarianceMatrix(const CovarianceMatrix &)  = default;
    CovarianceMatrix(      CovarianceMatrix &&) = default;
    CovarianceMatrix & operator=(const CovarianceMatrix &)  = default;
    CovarianceMatrix & operator=(      CovarianceMatrix &&) = default;

    ~CovarianceMatrix() = default;

    operator Matrix() const;
    Matrix get_inverse() const;

    Index nrows() const;
    Index ncols() const;

    /**! The number of diagonal blocks in the matrix excluding inverses.*/
    Index ndiagblocks() const;

    /**! The number of blocks in the matrix excluding inverses.*/
    Index nblocks() const;

    /**
     * Check if the block with indices (i, j) is contained in the
     * covariance matrix.
     *
     * @param i the block-row index
     * @param j the block-column index
     */
    bool has_block(Index i, Index j);

    /**
     * Return a pointer to the block with indices (i,j). If any
     * of i or j is less than zero, than the first block found in this
     * row or column is returned. This is useful for obtaining the element range
     * corresponding to a given block row or column range.
     *
     * @param i The block-row index of the block to return or -1.
     * @param j The block-column index of the block to return or -1.
     *
     * @return A pointer to the block at (i,j) or nullptr if
     *         this block doesn't exist.
     */
    const Block* get_block(Index i = -1, Index j = -1);

    /** Block in the covariance matrix.
     *
     * @return Reference to the std::vector holding the block
     * objects of this covariance matrix.
     */
    std::vector<Block>& get_blocks() {
        return correlations_;
    };

    /** Blocks of the inverse covariance matrix.
     *
     * @return Reference to the std::vector holding the blocks
     * objects of the inverse of the covariance matrix.
     */
    std::vector<Block>& get_inverse_blocks() {
        return inverses_;
    };

    /**
     * Checks that the covariance matrix contains one diagonal block per retrieval
     * quantity.
     *
     * TODO: This should be moved to m_retrieval
     *
     * @param jis The ArrayOfArrayOfIndex containing the first and last indices
     *        of each retrieval quantity in the state vector.
     * @return true if the covariance matrix contains a diagonal block for each
     *         retrieval quantity
     */
    bool has_diagonal_blocks(const ArrayOfArrayOfIndex &jis) const;

    /**
     * Checks that the dimensions of the blocks in the covariance matrix are
     * consistent with ranges occupied by the different retrieval quantities
     * in the state vector.
     *
     * TODO: This should be moved to m_retrieval
     *
     * @param jis The ArrayOfArrayOfIndex containing the first and last indices
     *        of each retrieval quantity in the state vector.
     * @return true if the covariance matrix contains a diagonal block for each
     *         retrieval quantity
     */
    bool is_consistent(const ArrayOfArrayOfIndex &jis) const;

    /**
     * This method checks whether a block is consistent with existing blocks
     * in the covariance matrix, i.e. that there is no other block in the
     * given row i (colum j) that has a different number of rows (columns)
     * than the given block.
     *
     * @param block The block for which to check consistent
     * @return true if the block is consistent with the other block
     *         contained in the matrix
     */
    bool is_consistent(const Block &block) const;

    /**
     * Compute the inverse of this correlation matrix. This function must be executed
     * after all block have been added to the covariance matrix and before any of the
     * mult_inv or add_inv methods is used.
     */
    void compute_inverse() const;

    /** Add block to covariance matrix.
     *
     * This function add a given block to the covariance matrix.
     * If this block is not consistent with other blocks in the matrix
     * an error will be thrown.
     *
     * @param c The block to add to the covariance matrix
     */
    void add_correlation(Block c);

    /** Add block inverse to covariance matrix.
     *
     * This function adds the inverse of a given block to a covariance
     * matrix. An error will be thrown if the corresponding non-inverse
     * block is not already in the covariance matrix.
     *
     * @param c The inverse of a block already in the matrix.
     */
    void add_correlation_inverse(Block c);

    /** Diagonal elements as vector
     *
     * Extracts the diagonal elements from the covariance matrix.
     * matrix.
     *
     * @return A vector containing the diagonal elements.
     */
    Vector diagonal() const;

    // Friend declarations.
    friend void mult(MatrixView, ConstMatrixView, const CovarianceMatrix &);
    friend void mult(MatrixView, const CovarianceMatrix &, ConstMatrixView);
    friend void mult(VectorView, const CovarianceMatrix &, ConstVectorView);

    friend void mult_inv(MatrixView, ConstMatrixView, const CovarianceMatrix &);
    friend void mult_inv(MatrixView, const CovarianceMatrix &, ConstMatrixView);
    friend void solve(VectorView, const CovarianceMatrix &, ConstVectorView);

    friend MatrixView & operator+=(MatrixView &, const CovarianceMatrix &);
    friend void add_inv(MatrixView, const CovarianceMatrix &);

    friend void xml_read_from_stream(istream &, CovarianceMatrix &, bifstream *,
                                     const Verbosity &);
    friend void xml_write_to_stream(ostream &, const CovarianceMatrix &, bofstream *,
                                    const String&, const Verbosity &);
    friend std::ostream& operator<<(std::ostream& os, const CovarianceMatrix& v);
private:

    void generate_blocks(std::vector<std::vector<const Block *>> &) const;
    void invert_correlation_block(std::vector<Block> &inverses,
                                  std::vector<const Block *> &blocks) const;
    bool has_inverse(IndexPair indices) const;

    std::vector<Block> correlations_;
    mutable std::vector<Block> inverses_;

};

void mult(MatrixView, ConstMatrixView, const CovarianceMatrix &);
void mult(MatrixView, const CovarianceMatrix &, ConstMatrixView);
void mult(VectorView, const CovarianceMatrix &, ConstVectorView);

void mult_inv(MatrixView, ConstMatrixView, const CovarianceMatrix &);
void mult_inv(MatrixView, const CovarianceMatrix &, ConstMatrixView);
void solve(VectorView, const CovarianceMatrix &, ConstVectorView);

MatrixView & operator+=(MatrixView &, const CovarianceMatrix &);
void add_inv(MatrixView, const CovarianceMatrix &);

std::ostream& operator<<(std::ostream& os, const ConstVectorView& v);

#endif // covariance_matrix_h
