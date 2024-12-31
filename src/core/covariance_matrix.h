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

#include <matpack.h>
#include <mystring.h>

#include <iosfwd>
#include <memory>
#include <utility>

class CovarianceMatrix;

//------------------------------------------------------------------------------
// Type Aliases
//------------------------------------------------------------------------------

using IndexPair = std::pair<Index, Index>;

class BlockMatrix {
 public:
  using variant_t =
      std::variant<std::shared_ptr<Matrix>, std::shared_ptr<Sparse>>;

 private:
  variant_t data;

 public:
  BlockMatrix()                               = default;
  BlockMatrix(const BlockMatrix &)            = default;
  BlockMatrix(BlockMatrix &&)                 = default;
  BlockMatrix &operator=(const BlockMatrix &) = default;
  BlockMatrix &operator=(BlockMatrix &&)      = default;

  BlockMatrix(std::shared_ptr<Matrix> dense) : data(std::move(dense)) {}
  BlockMatrix(std::shared_ptr<Sparse> sparse) : data(std::move(sparse)) {}
  BlockMatrix(const Matrix &dense) : data(std::make_shared<Matrix>(dense)) {}
  BlockMatrix(const Sparse &sparse) : data(std::make_shared<Sparse>(sparse)) {}

  BlockMatrix &operator=(std::shared_ptr<Matrix> dense);

  BlockMatrix &operator=(std::shared_ptr<Sparse> sparse);

  BlockMatrix &operator=(const Matrix &dense);

  BlockMatrix &operator=(const Sparse &sparse);

  [[nodiscard]] bool not_null() const;

  [[nodiscard]] bool is_dense() const;

  [[nodiscard]] bool is_sparse() const;

  [[nodiscard]] Matrix &dense();

  [[nodiscard]] const Matrix &dense() const;

  [[nodiscard]] Sparse &sparse();

  [[nodiscard]] const Sparse &sparse() const;

  [[nodiscard]] Vector diagonal() const;

  [[nodiscard]] Index ncols() const;

  [[nodiscard]] Index nrows() const;

  friend std::ostream &operator<<(std::ostream &os, const BlockMatrix &m);

  [[nodiscard]] std::array<Index, 2> shape() const;
};

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
     * Create a correlation block from given row_range, column_range and shared_ptr
     * to a dense matrix.
     *
     * @param row_range The range of element-rows covered by the block
     * @param column_range The range of element-column covered by the block
     * @param indices The pair of block-indices that identifies the block in the
     *        covariance matrix.
     * @param dense A shared pointer to a den matrix of type Matrix
     */
  Block(Range row_range,
        Range column_range,
        IndexPair indices,
        BlockMatrix matrix)
      : row_range_(row_range),
        column_range_(column_range),
        indices_(std::move(indices)),
        matrix_(std::move(matrix)) {
    // Nothing to do here.
  }

  Block(const Block &)            = default;
  Block(Block &&)                 = default;
  Block &operator=(const Block &) = default;
  Block &operator=(Block &&)      = default;

  ~Block() = default;

  /*! Number of rows of this block */
  [[nodiscard]] Index nrows() const { return row_range_.nelem; }
  /*! Number of columns of this block */
  [[nodiscard]] Index ncols() const { return column_range_.nelem; }

  /*! The row range of this block*/
  [[nodiscard]] Range get_row_range() const { return row_range_; }
  /*! The column range of this block*/
  [[nodiscard]] Range get_column_range() const { return column_range_; }

  /*! The row range of this block*/
  Range &get_row_range() { return row_range_; }
  /*! The column range of this block*/
  Range &get_column_range() { return column_range_; }

  void set_matrix(std::shared_ptr<Sparse> sparse) {
    matrix_ = std::move(sparse);
  }
  void set_matrix(std::shared_ptr<Matrix> dense) { matrix_ = std::move(dense); }

  /*! Return the diagonal as a vector.*/
  [[nodiscard]] Vector diagonal() const { return matrix_.diagonal(); }

  /*! Return the indices of the retrieval quantities correlated by this block as std::pair. */
  [[nodiscard]] IndexPair get_indices() const { return indices_; }

  /*! Return the indices of the retrieval quantities correlated by this block as std::pair. */
  void set_indices(Index f, Index s) { indices_ = {f, s}; }

  [[nodiscard]] bool not_null() const { return matrix_.not_null(); }
  [[nodiscard]] bool is_dense() const { return matrix_.is_dense(); }
  [[nodiscard]] bool is_sparse() const { return matrix_.is_sparse(); }

  [[nodiscard]] const Matrix &get_dense() const { return matrix_.dense(); }
  Matrix &get_dense() { return matrix_.dense(); }

  [[nodiscard]] const Sparse &get_sparse() const { return matrix_.sparse(); }
  Sparse &get_sparse() { return matrix_.sparse(); }

  // Friend declarations.
  friend void mult(StridedMatrixView, StridedConstMatrixView, const Block &);
  friend void mult(StridedMatrixView, const Block &, StridedConstMatrixView);
  friend void mult(StridedVectorView, const Block &, StridedConstVectorView);

  friend StridedMatrixView operator+=(StridedMatrixView, const Block &);

 private:
  Range row_range_, column_range_;
  IndexPair indices_;

  BlockMatrix matrix_;
};

void mult(StridedMatrixView, StridedConstMatrixView, const Block &);
void mult(StridedMatrixView, const Block &, StridedConstMatrixView);
void mult(StridedVectorView, const Block &, StridedConstVectorView);

StridedMatrixView operator+=(StridedMatrixView, const Block &);
void add_inv(StridedMatrixView A, const Block &);

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
  CovarianceMatrix()                                    = default;
  CovarianceMatrix(const CovarianceMatrix &)            = default;
  CovarianceMatrix(CovarianceMatrix &&)                 = default;
  CovarianceMatrix &operator=(const CovarianceMatrix &) = default;
  CovarianceMatrix &operator=(CovarianceMatrix &&)      = default;

  ~CovarianceMatrix() = default;

  explicit operator Matrix() const;
  Matrix get_inverse() const;

  Index nrows() const;
  Index ncols() const;

  /**! The number of diagonal blocks in the matrix excluding inverses.*/
  Index ndiagblocks() const;

  /**! The number of inverse diagonal blocks in the matrix.*/
  Index ninvdiagblocks() const;

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
  const Block *get_block(Index i = -1, Index j = -1);

  /** Block in the covariance matrix.
     *
     * @return Reference to the std::vector holding the block
     * objects of this covariance matrix.
     */
  std::vector<Block> &get_blocks() { return correlations_; };

  /** Blocks of the inverse covariance matrix.
     *
     * @return Reference to the std::vector holding the blocks
     * objects of the inverse of the covariance matrix.
     */
  std::vector<Block> &get_inverse_blocks() { return inverses_; };

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

  /** Diagonal of the inverse of the covariance matrix as vector
     *
     * Extracts the diagonal elements from the inverse of the covariance matrix.
     * This can trigger the computation of the inverse of the matrix, if it has
     * not been provided by the user.
     *
     * @return A vector containing the diagonal elements.
     */
  Vector inverse_diagonal() const;

  // Friend declarations.
  friend void mult(StridedMatrixView, StridedConstMatrixView, const CovarianceMatrix &);
  friend void mult(StridedMatrixView, const CovarianceMatrix &, StridedConstMatrixView);
  friend void mult(StridedVectorView, const CovarianceMatrix &, StridedConstVectorView);

  friend void mult_inv(StridedMatrixView, StridedConstMatrixView, const CovarianceMatrix &);
  friend void mult_inv(StridedMatrixView, const CovarianceMatrix &, StridedConstMatrixView);
  friend void solve(StridedVectorView, const CovarianceMatrix &, StridedConstVectorView);

  friend StridedMatrixView operator+=(StridedMatrixView, const CovarianceMatrix &);
  friend void add_inv(StridedMatrixView, const CovarianceMatrix &);

  friend void xml_read_from_stream(std::istream &,
                                   CovarianceMatrix &,
                                   bifstream *);
  friend void xml_write_to_stream(std::ostream &,
                                  const CovarianceMatrix &,
                                  bofstream *,
                                  const String &);
  friend std::ostream &operator<<(std::ostream &os, const CovarianceMatrix &v);

 private:
  void generate_blocks(std::vector<std::vector<const Block *>> &) const;
  void invert_correlation_block(std::vector<Block> &inverses,
                                std::vector<const Block *> &blocks) const;
  bool has_inverse(IndexPair indices) const;

  std::vector<Block> correlations_;
  mutable std::vector<Block> inverses_;
};

void mult(StridedMatrixView, StridedConstMatrixView, const CovarianceMatrix &);
void mult(StridedMatrixView, const CovarianceMatrix &, StridedConstMatrixView);
void mult(StridedVectorView, const CovarianceMatrix &, StridedConstVectorView);

void mult_inv(StridedMatrixView, StridedConstMatrixView, const CovarianceMatrix &);
void mult_inv(StridedMatrixView, const CovarianceMatrix &, StridedConstMatrixView);
void solve(StridedVectorView, const CovarianceMatrix &, StridedConstVectorView);

StridedMatrixView operator+=(StridedMatrixView, const CovarianceMatrix &);
void add_inv(StridedMatrixView, const CovarianceMatrix &);

template <>
struct std::formatter<BlockMatrix> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const BlockMatrix &v, FmtContext &ctx) const {
    if (v.not_null())
      return v.is_dense() ? tags.format(ctx, v.dense())
                          : tags.format(ctx, v.sparse());
    tags.add_if_bracket(ctx, '[');
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<CovarianceMatrix> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const CovarianceMatrix &v,
                              FmtContext &ctx) const {
    const std::string_view sep = tags.sep(true);
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, Matrix{v}, sep, v.get_inverse());
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

#endif  // covariance_matrix_h
