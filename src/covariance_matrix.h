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
 * Matrix or type Sparse. In addition to this it holds to two indices i and j
 * and two range objects that  describe the position of the block in the
 * covariance matrix. The indices i and j are the indices of the retrieval
 * quantities that this block relates.
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
     * to matrix object.
     */
    Block(Range row_range, Range column_range, IndexPair indices,
                std::shared_ptr<Matrix> dense)
        : row_range_(row_range), column_range_(column_range), indices_(indices),
        matrix_type_(MatrixType::dense), dense_(dense), sparse_(nullptr)
    {
        // Nothing to do here.
    }

    /*
     * Same as above but constructs a Block holding a matrix of type Sparse.
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

    /*! Return the indices of the retrieval quantities correlated by this block as std::pair. */
    IndexPair  get_indices()         const {return indices_;}
    /*! Return the type of the matrix holding the correlation coefficients. */
    MatrixType get_matrix_type()     const {return matrix_type_;}

    const Matrix & get_dense()  const
    {
        assert(dense_);
        return *dense_;
    }
    const Sparse & get_sparse() const
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
 * The covariance matrix class implements provides overloads for the common
 * ARTS operations for multiplication by matrices and vectors as well as
 * adding covariance matrices to other matrices.
 *
 * Computing inverses of covariance matrices is handled indirectly by providing
 * mult_inv methods that multiply the inverse of the covariance matrix by a given
 * vector or matrix. This, however, requires previously having computed the inverse
 * of the matrix using the comput_inverse method.
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

    Index nrows() const;
    Index ncols() const;
    Index ndiagblocks() const;

    /**
     * Checks that the covariance matrix contains one diagonal block per retrieval
     * quantity.
     */
    bool has_diagonal_blocks(const ArrayOfArrayOfIndex &jis) const;

    /**
     * Checks if block positions and dimensions agree with the jacobian.
     */
    bool is_consistent(const ArrayOfArrayOfIndex &jis) const;

    /**
     * Checks if block positions and dimensions agree with the jacobian.
     */

    /**
     * Compute the inverse of this correlation matrix. This function must be executed
     * after all block have been added to the covariance matrix and before any of the
     * mult_inv or add_inv methods is used.
     */
    void compute_inverse() const;

    /** Add correlation in form of given block to this matrix. */
    void add_correlation(Block c);
    /** Add the inverse of a correlation block to this matrix. This
     *  will only work if the original block is also contained in the
     *  covariance matrix. */
    void add_correlation_inverse(Block c);
    /** Add a diagonal block describing the internal correlations of retrieval quantity
     * with index i to the matrix. A runtime error
     * i and j given by a dense matrix of type Matrix. This function requires that
     * that the diagonal blocks corresponding to both quantities i and j have been
     * added to the covariance matrix previously otherwise a runtime error will
     * be thrown
     * @param i The index of the first retrieval quantity.
     * @param j The index of the second retrieval quantity.
     * @param matrix A shared pointer to the matrix holding the correlation coefficients.
     */
    template<typename MatrixType>
    void add_correlation(Index i,
                         std::shared_ptr<MatrixType> matrix);
    /** Add a correlation between two retrieval quantities that have been previously
     * added to the matrix.
     * @param i The index of the first retrieval quantity.
     * @param j The index of the second retrieval quantity.
     * @param matrix A shared pointer to the matrix holding the correlation coefficients.
     */
    template<typename MatrixType>
    void add_correlation(Index i,
                         Index j,
                         std::shared_ptr<MatrixType> matrix);
    /** Add a correlation relation between the retrieval quantitites with indices
     * i and j represented by a dense or a sparse matrix.
     * @param i The index of the first retrieval quantity.
     * @param j The index of the second retrieval quantity.
     * @param jacobian_indices The ArrayOfArrayOfIndex containing the first and last
     * indices of all retrieval quantities.
     * @param matrix A shared pointer to the matrix holding the correlation coefficients.
     */
    template<typename MatrixType>
    void add_correlation(Index i,
                         Index j,
                         const ArrayOfArrayOfIndex & jacobian_indices,
                         std::shared_ptr<MatrixType> sparse);

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
    bool has_inverse(std::pair<Index, Index> indices) const;

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
