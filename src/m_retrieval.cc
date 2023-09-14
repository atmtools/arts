#include <algorithm>

#include <workspace.h>
#include "covariance_matrix.h"
#include "jacobian.h"

////////////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////////////

void check_and_add_block(CovarianceMatrix& covmat,
                         const RetrievalQuantity& jq,
                         const Index rq_index,
                         const Index grid_dimensions,
                         const Sparse& covmat_block,
                         const Sparse& covmat_inv_block) {
  Index start = covmat.nrows();
  Index extent = covmat_block.nrows();
  Range range(start, extent);
  Index n_gps = 1;
  for (Index j = 0; j < grid_dimensions; ++j) {
    n_gps *= jq.Grids()[j].size();
  }

  if (!covmat_block.empty()) {
    if ((n_gps == extent) && (n_gps == covmat_block.ncols())) {
      std::shared_ptr<Sparse> mat = std::make_shared<Sparse>(covmat_block);
      covmat.add_correlation(
          Block(range, range, std::make_pair(rq_index, rq_index), mat));
    } else {
      ARTS_USER_ERROR (
        "The matrix in covmat_block was expected to have dimensions [",
        n_gps, ", ", n_gps, "] but found  to have dimensions [",
        covmat_block.nrows(), ", ", covmat_block.ncols(), "].")
    }
  }
  if (!covmat_inv_block.empty()) {
    if ((n_gps == covmat_inv_block.nrows()) &&
        (n_gps == covmat_inv_block.ncols())) {
      std::shared_ptr<Sparse> mat = std::make_shared<Sparse>(covmat_inv_block);
      covmat.add_correlation_inverse(
          Block(range, range, std::make_pair(rq_index, rq_index), mat));
    } else {
      ARTS_USER_ERROR (
        "The matrix in covmat_inv_block was expected to have dimensions [",
        n_gps, ", ", n_gps, "] but found  to have dimensions [",
        covmat_block.nrows(), ", ", covmat_block.ncols(), "].")
    }
  }
}

void add_scalar_variance(CovarianceMatrix& covmat,
                         ArrayOfRetrievalQuantity& jacobian_quantities,
                         Numeric var) {
  Index i = jacobian_quantities.size() - 1;
  Index start = covmat.nrows();
  Index extent = 1;
  Range range(start, extent);
  std::shared_ptr<Matrix> mat = std::make_shared<Matrix>(1, 1);
  mat->operator()(0, 0) = var;
  covmat.add_correlation(Block(range, range, std::make_pair(i, i), mat));
  mat = std::make_shared<Matrix>(1, 1);
  mat->operator()(0, 0) = 1.0 / var;
  covmat.add_correlation_inverse(
      Block(range, range, std::make_pair(i, i), mat));
}

////////////////////////////////////////////////////////////////////////////////
// Simple Pressure-Height Conversion
////////////////////////////////////////////////////////////////////////////////

void ZFromPSimple(Vector& z_grid, const Vector& p_grid) {
  z_grid = Vector(p_grid.size());

  for (Index i = 0; i < p_grid.size(); ++i) {
    ARTS_USER_ERROR_IF (p_grid[i] < 0.01,
                        "Pressures below 0.01 Pa are not accedpted.");
  }

  for (Index i = 0; i < p_grid.size(); ++i) {
    z_grid[i] = 16e3 * (5.0 - log10(p_grid[i]));
  }
}

void PFromZSimple(Vector& p_grid, const Vector& z_grid) {
  p_grid = Vector(z_grid.size());

  for (Index i = 0; i < p_grid.size(); ++i) {
    ARTS_USER_ERROR_IF (z_grid[i] > 120e3,
                        "Altitudes above 120 km are not accepted.");
  }

  for (Index i = 0; i < z_grid.size(); ++i) {
    p_grid[i] = pow(10.0, 5.0 - z_grid[i] / 16e3);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Creation of Correlation Blocks
////////////////////////////////////////////////////////////////////////////////

template <typename MatrixType>
void insert_elements(MatrixType& matrix,
                     Index /*m*/,
                     Index /*n*/,
                     const ArrayOfIndex& row_indices,
                     const ArrayOfIndex& column_indices,
                     const Vector& elements) {
  matrix.insert_elements(
      row_indices.size(), row_indices, column_indices, elements);
}

template <>
void insert_elements(Matrix& matrix,
                     Index m,
                     Index n,
                     const ArrayOfIndex& row_indices,
                     const ArrayOfIndex& column_indices,
                     const Vector& elements) {
  ARTS_ASSERT(row_indices.size() == column_indices.size());
  ARTS_ASSERT(column_indices.size() == static_cast<Size>(elements.size()));

  matrix.resize(m, n);

  for (Index i = 0; i < elements.size(); ++i) {
    matrix(row_indices[i], column_indices[i]) = elements[i];
  }
}

template <typename MatrixType>
void covmatDiagonal(MatrixType& block,
                    MatrixType& block_inv,
                    const Vector& vars) {
  ARTS_USER_ERROR_IF (vars.empty(),
        "Cannot pass empty vector of variances to covmat_blockSetDiagonal");
  Index n = vars.size();
  block = MatrixType(n, n);
  block_inv = MatrixType(n, n);

  ArrayOfIndex indices(n);
  Vector elements(n), elements_inv(n);
  for (Index i = 0; i < n; ++i) {
    indices[i] = i;
    elements[i] = vars[i];
    elements_inv[i] = 1.0 / vars[i];
  }

  insert_elements(block, n, n, indices, indices, elements);
  insert_elements(block_inv, n, n, indices, indices, elements_inv);
}

void covmatDiagonal(Matrix& block, Matrix& block_inv, const Vector& vars) {
  covmatDiagonal<Matrix>(block, block_inv, vars);
}

void covmatDiagonal(Sparse& block, Sparse& block_inv, const Vector& vars) {
  covmatDiagonal<Sparse>(block, block_inv, vars);
}

template <typename MatrixType>
void covmat1DMarkov(MatrixType& block,
                    MatrixType& block_inv,
                    const Vector& grid,
                    const Vector& sigma,
                    const Numeric& lc,
                    const Numeric& /*co*/) {
  Index n = grid.size();
  ArrayOfIndex row_indices, column_indices;
  Vector elements, elements_inv;

  ARTS_USER_ERROR_IF (n != sigma.size(),
                      "Size of grid incompatible with given variances.");

  elements = Vector(row_indices.size());
  for (size_t i = 0; i < row_indices.size(); ++i) {
    Numeric dz = abs(grid[row_indices[i]] - grid[column_indices[i]]);
    Numeric e =
        sigma[row_indices[i]] * sigma[column_indices[i]] * exp(-dz / lc);
    elements[i] = e;
  }

  block = MatrixType(n, n);
  insert_elements(block, n, n, row_indices, column_indices, elements);
  row_indices = ArrayOfIndex{};
  column_indices = ArrayOfIndex{};
  elements = Vector(3 * n - 2);

  Numeric dz = abs(grid[1] - grid[0]);
  Numeric alpha = exp(-dz / lc);
  Numeric c1 = -alpha / (1.0 - alpha * alpha);
  Numeric c2 = 1.0 / (1.0 - alpha * alpha);

  for (Index i = 0; i < n; ++i) {
    // Lower Tri-diagonal
    if (i > 0) {
      column_indices.push_back(i - 1);
      row_indices.push_back(i);
      elements[i * 3 - 1] = c1;
      elements[i * 3 - 1] /= (sigma[i] * sigma[i - 1]);
    }

    // Diagonal Elements
    column_indices.push_back(i);
    row_indices.push_back(i);
    if (i == 0 || i == n - 1) {
      elements[i * 3] = c2;
    } else {
      elements[i * 3] = c2 * (1.0 + alpha * alpha);
    }
    elements[i * 3] /= (sigma[i] * sigma[i]);

    // Upper Tri-diagonal
    if (i < n - 1) {
      column_indices.push_back(i + 1);
      row_indices.push_back(i);
      elements[i * 3 + 1] = c1;
      elements[i * 3 + 1] /= (sigma[i] * sigma[i + 1]);
    }
  }
  block_inv = MatrixType(n, n);
  insert_elements(block_inv, n, n, row_indices, column_indices, elements);
}

void covmat1DMarkov(Matrix& block,
                    Matrix& block_inv,
                    const Vector& grid,
                    const Vector& sigma,
                    const Numeric& lc,
                    const Numeric& /*co*/) {
  covmat1DMarkov<Matrix>(block, block_inv, grid, sigma, lc, 0.0);
}

void covmat1DMarkov(Sparse& block,
                    Sparse& block_inv,
                    const Vector& grid,
                    const Vector& sigma,
                    const Numeric& lc,
                    const Numeric& /*co*/) {
  covmat1DMarkov<Sparse>(block, block_inv, grid, sigma, lc, 0.0);
}

template <typename MatrixType>
void covmat1D(MatrixType& block,
              const Vector& grid1,
              const Vector& grid2,
              const Vector& sigma1,
              const Vector& sigma2,
              const Vector& lc1,
              const Vector& lc2,
              const Numeric& co,
              const String& fname) {
  Index m = grid1.size();
  Vector sigma1_copy(sigma1), lc1_copy(lc1);

  ARTS_ASSERT((sigma1.size() == m) || (sigma1.size() == 1));
  if (sigma1.size() == 1) {
    Numeric v = sigma1[0];
    sigma1_copy = Vector(m);
    sigma1_copy = v;
  }

  ARTS_ASSERT((lc1.size() == m) || (lc1.size() == 1));
  if (lc1.size() == 1) {
    Numeric v = lc1[0];
    lc1_copy = Vector(m);
    lc1_copy = v;
  };

  Index n = grid2.size();
  Vector sigma2_copy(sigma2), lc2_copy(lc2);

  ARTS_ASSERT((sigma2.size() == n) || (sigma2.size() == 1));
  if (sigma2.size() == 1) {
    Numeric v = sigma2[0];
    sigma2_copy = Vector(n);
    sigma2_copy = v;
  }

  ARTS_ASSERT((lc2.size() == n) || (lc2.size() == 1));
  if (lc2.size() == 1) {
    Numeric v = lc2[0];
    lc2_copy = Vector(n);
    lc2_copy = v;
  };

  ConstVectorView grid_view_1(grid1), sigma_view_1(sigma1_copy),
      cls_view_1(lc1_copy);
  ConstVectorView grid_view_2(grid2), sigma_view_2(sigma2_copy),
      cls_view_2(lc2_copy);

  if (n == 0) {
    n = m;
    grid_view_2.set(grid_view_1);
    sigma_view_2.set(sigma_view_1);
    cls_view_2.set(cls_view_1);
  }

  // Correlation Functions
  auto f_lin = [&](Index i, Index j) {
    Numeric d = abs(grid_view_1[i] - grid_view_2[j]);
    Numeric cl = 0.5 * (cls_view_1[i] + cls_view_2[j]);
    return 1.0 - (1.0 - exp(-1.0)) * (d / cl);
  };

  auto f_exp = [&](Index i, Index j) {
    Numeric d = abs(grid_view_1[i] - grid_view_2[j]);
    Numeric cl = 0.5 * (cls_view_1[i] + cls_view_2[j]);
    return exp(-d / cl);
  };

  auto f_gau = [&](Index i, Index j) {
    Numeric d = abs(grid_view_1[i] - grid_view_2[j]);
    Numeric cl = 0.5 * (cls_view_1[i] + cls_view_2[j]);
    Numeric x = d / cl;
    x *= x;
    return exp(-x);
  };

  ArrayOfIndex row_indices;
  ArrayOfIndex column_indices;

  row_indices.reserve(n * m);
  column_indices.reserve(n * m);

  std::function<Numeric(Index, Index)> f;

  if (fname == "exp") {
    f = f_exp;
  } else if (fname == "lin") {
    f = f_lin;
  } else if (fname == "gau") {
    f = f_gau;
  } else {
    std::ostringstream os;
    os << fname << " is not a known function name. Supported names"
       << "are: exp, lin, gau.";
    std::runtime_error(os.str());
  }

  for (Index i = 0; i < m; ++i) {
    for (Index j = 0; j < n; ++j) {
      Numeric e = f(i, j);
      if (e >= co) {
        row_indices.push_back(i);
        column_indices.push_back(j);
      }
    }
  }

  Vector elements(row_indices.size());
  for (size_t i = 0; i < row_indices.size(); ++i) {
    Index ii = row_indices[i];
    Index jj = column_indices[i];
    elements[i] = sigma_view_1[ii] * sigma_view_2[jj] * f(ii, jj);
  }

  block = MatrixType(m, n);
  insert_elements(block, m, n, row_indices, column_indices, elements);
}

void covmat1D(Matrix& block,
              const Vector& grid1,
              const Vector& grid2,
              const Vector& sigma1,
              const Vector& sigma2,
              const Vector& lc1,
              const Vector& lc2,
              const Numeric& co,
              const String& fname) {
  covmat1D<Matrix>(block, grid1, grid2, sigma1, sigma2, lc1, lc2, co, fname);
}

void covmat1D(Sparse& block,
              const Vector& grid1,
              const Vector& grid2,
              const Vector& sigma1,
              const Vector& sigma2,
              const Vector& lc1,
              const Vector& lc2,
              const Numeric& co,
              const String& fname) {
  covmat1D<Sparse>(block, grid1, grid2, sigma1, sigma2, lc1, lc2, co, fname);
}

////////////////////////////////////////////////////////////////////////////////
// Manipulation of covmat_se and covmat_sx
////////////////////////////////////////////////////////////////////////////////

template <typename MatrixType>
void covmat_seAddBlock(CovarianceMatrix& covmat_se,
                       const MatrixType& block,
                       const Index& i,
                       const Index& j) {
  Index m = block.nrows();
  Index n = block.ncols();

  Index ii(i), jj(j);
  if ((ii < 0) && (jj < 0)) {
    ii = covmat_se.ndiagblocks();
    jj = ii;
  }

  ARTS_USER_ERROR_IF (j < i,
        "The block must be on or above the diagonal, "
        " i.e. *i* <= *j*.");

  Index ndiagblocks = covmat_se.ndiagblocks();

  if ((ii >= ndiagblocks)) {
    if (ii > ndiagblocks) {
      ARTS_USER_ERROR_IF (jj > ii,
            "Off-diagonal block can only be added to rows that already "
            "have a block on the diagonal.");
      ARTS_USER_ERROR (
            "Diagonal block must be added row-by-row starting in the "
            " upper left of the matrix.");
    }
  }

  ARTS_USER_ERROR_IF (covmat_se.has_block(ii, jj),
                      "Block already present in covariance matrix.");

  if (ii == jj) {
    ARTS_USER_ERROR_IF (m != n,
                        "Diagonal blocks must be square.");
    Index start = covmat_se.nrows();
    Range range(start, m);
    std::shared_ptr<MatrixType> mat = std::make_shared<MatrixType>(block);
    covmat_se.add_correlation(Block(range, range, std::make_pair(ii, ii), mat));

  } else {
    const Block* b = covmat_se.get_block(ii, ii);
    ARTS_USER_ERROR_IF (!b,
          "Trying to add an off-diagonal block that"
          " lacks corresponding diagonal block in the "
          " same row.");
    Range row_range = b->get_row_range();

    b = covmat_se.get_block(jj, jj);
    ARTS_USER_ERROR_IF (!b,
          "Trying to add an off-diagonal block that"
          " lacks corresponding diagonal block in the "
          " same column.");
    Range column_range = b->get_column_range();

    ARTS_USER_ERROR_IF ((row_range.extent != m) || (column_range.extent != n),
          "The off-diagonal block is inconsistent "
          "with the corresponding diagonal blocks.");

    std::shared_ptr<MatrixType> mat = std::make_shared<MatrixType>(block);
    covmat_se.add_correlation(
        Block(row_range, column_range, std::make_pair(ii, jj), mat));
  }
}

void covmat_seAddBlock(CovarianceMatrix& covmat_se,
                       const Matrix& block,
                       const Index& i,
                       const Index& j) {
  covmat_seAddBlock<Matrix>(covmat_se, block, i, j);
}

void covmat_seAddBlock(CovarianceMatrix& covmat_se,
                       const Sparse& block,
                       const Index& i,
                       const Index& j) {
  covmat_seAddBlock<Sparse>(covmat_se, block, i, j);
}

template <typename MatrixType>
void covmat_seAddInverseBlock(CovarianceMatrix& covmat_se,
                              const MatrixType& inv_block,
                              const Index& i,
                              const Index& j) {
  ARTS_USER_ERROR_IF (covmat_se.ndiagblocks() == 0,
                      "Need at least one non-inverse block in the matrix"
                      " before an inverse block can be added.");
  Index ii(i), jj(j);
  if ((ii < 0) && (jj < 0)) {
    ii = covmat_se.ndiagblocks() - 1;
    jj = ii;
  }

  Index m = inv_block.nrows();
  Index n = inv_block.ncols();

  const Block* block = covmat_se.get_block(ii, jj);

  ARTS_USER_ERROR_IF (!block,
        "Cannot add inverse  block to the covariance "
        " without corresponding non-inverse block.");

  ARTS_USER_ERROR_IF ((m != inv_block.nrows()) || (n != inv_block.ncols()),
        "Dimensions of block are inconsistent with "
        " non-inverse block.");

  Range row_range = block->get_row_range();
  Range column_range = block->get_column_range();

  std::shared_ptr<MatrixType> mat = std::make_shared<MatrixType>(inv_block);
  covmat_se.add_correlation_inverse(
      Block(row_range, column_range, std::make_pair(ii, jj), mat));
}

void covmat_seAddInverseBlock(CovarianceMatrix& covmat_se,
                              const Matrix& block,
                              const Index& i,
                              const Index& j) {
  covmat_seAddInverseBlock<Matrix>(covmat_se, block, i, j);
}

void covmat_seAddInverseBlock(CovarianceMatrix& covmat_se,
                              const Sparse& block,
                              const Index& i,
                              const Index& j) {
  covmat_seAddInverseBlock<Sparse>(covmat_se, block, i, j);
}

template <typename MatrixType>
void setCovarianceMatrix(CovarianceMatrix& covmat, const MatrixType& block) {
  Index m = block.nrows();
  Index n = block.ncols();

  ARTS_USER_ERROR_IF (n != m,
                      "Covariance matrix must be square!");

  covmat = CovarianceMatrix();
  IndexPair indices = std::make_pair(0, 0);
  Range range = Range(0, n);
  std::shared_ptr<MatrixType> mat_ptr = std::make_shared<MatrixType>(block);
  covmat.add_correlation(Block(range, range, indices, mat_ptr));
}

template <>
void setCovarianceMatrix(CovarianceMatrix& covmat,
                         const CovarianceMatrix& block) {
  covmat = block;
}

template <typename MatrixType>
void covmat_seSet(CovarianceMatrix& covmat,
                  const MatrixType &block) {
  setCovarianceMatrix(covmat, block);
}

void covmat_seSet(CovarianceMatrix& covmat, const CovarianceMatrix& x) {
  covmat_seSet<CovarianceMatrix>(covmat, x);
}

void covmat_seSet(CovarianceMatrix& covmat, const Matrix& x) {
  covmat_seSet<Matrix>(covmat, x);
}

void covmat_seSet(CovarianceMatrix& covmat, const Sparse& x) {
  covmat_seSet<Sparse>(covmat, x);
}

template <typename MatrixType>
void covmat_sxSet(CovarianceMatrix& covmat,
                  const MatrixType &block) {
  setCovarianceMatrix(covmat, block);
}

void covmat_sxSet(CovarianceMatrix& covmat, const CovarianceMatrix& x) {
  covmat_sxSet<CovarianceMatrix>(covmat, x);
}

void covmat_sxSet(CovarianceMatrix& covmat, const Matrix& x) {
  covmat_sxSet<Matrix>(covmat, x);
}

void covmat_sxSet(CovarianceMatrix& covmat, const Sparse& x) {
  covmat_sxSet<Sparse>(covmat, x);
}

template <typename MatrixType>
void covmat_sxAddBlock(CovarianceMatrix& covmat_sx,
                       const ArrayOfRetrievalQuantity& jq,
                       const MatrixType& block,
                       const Index& i,
                       const Index&j ) {
  Index ii{i}, jj{j};
  if ((ii < 0) && (jj < 0)) {
    ii = covmat_sx.ndiagblocks();
    jj = ii;
    ARTS_USER_ERROR_IF ((static_cast<Size>(ii) >= jq.size()) || (static_cast<Size>(jj) >= jq.size()),
          "*covmat_sx* already contains more or as many diagonal"
          " blocks as there are retrieval quantities.");
  } else {
    ARTS_USER_ERROR_IF ((static_cast<Size>(ii) >= jq.size()) || (static_cast<Size>(jj) >= jq.size()),
        "The block indices must either be both -1 (default) or\n"
        "non-negative and smaller than the number of retrieval \n"
        "quantities.");
    ARTS_USER_ERROR_IF (ii > jj,
          "Only blocks above the diagonal can be set, hence"
          "*i* must be less than or equal to *j*.");
  }

  Index m = block.nrows();
  Index n = block.ncols();

  Index jq_m = jq[ii].size();
  if (jq[ii].HasAffine()) {
    jq_m = jq[ii].TransformationMatrix().ncols();
  }
  Index jq_n = jq[jj].size();
  if (jq[jj].HasAffine()) {
    jq_n = jq[jj].TransformationMatrix().ncols();
  }

  ARTS_USER_ERROR_IF ((m != jq_m) || (n != jq_n),
      "The covariance block for retrieval quantities ", ii, " and ",
      jj, " was expected to have dimensions [", jq_m, ", ", jq_n,
      "] but found to have dimensions "
      "[", m, ", ", n, "].")

  ArrayOfArrayOfIndex ji;
  bool any_affine;
  jac_ranges_indices(ji, any_affine, jq);
  Index row_start = ji[ii][0];
  Index row_extent = ji[ii][1] - ji[ii][0] + 1;
  Range row_range(row_start, row_extent);

  Index col_start = ji[jj][0];
  Index col_extent = ji[jj][1] - ji[jj][0] + 1;
  Range col_range(col_start, col_extent);

  std::shared_ptr<MatrixType> mat = std::make_shared<MatrixType>(block);
  covmat_sx.add_correlation(
      Block(row_range, col_range, std::make_pair(ii, jj), mat));
}

void covmat_sxAddBlock(CovarianceMatrix& covmat_sx,
                       const ArrayOfRetrievalQuantity& jq,
                       const Matrix& block,
                       const Index& i,
                       const Index& j) {
  covmat_sxAddBlock<Matrix>(covmat_sx, jq, block, i, j);
}

void covmat_sxAddBlock(CovarianceMatrix& covmat_sx,
                       const ArrayOfRetrievalQuantity& jq,
                       const Sparse& block,
                       const Index& i,
                       const Index& j) {
  covmat_sxAddBlock<Sparse>(covmat_sx, jq, block, i, j);
}

template <typename MatrixType>
void covmat_sxAddInverseBlock(CovarianceMatrix& covmat_sx,
                              const ArrayOfRetrievalQuantity& jq,
                              const MatrixType& block_inv,
                              const Index& i,
                              const Index&j ) {
  ARTS_USER_ERROR_IF (covmat_sx.ndiagblocks() == 0,
                      "Need at least one non-inverse block in the matrix"
                      " before an inverse block can be added.");
  Index ii(i), jj(j);
  if ((ii < 0) && (jj < 0)) {
    ii = covmat_sx.ndiagblocks() - 1;
    jj = ii;
  } else {
    ARTS_USER_ERROR_IF ((static_cast<Size>(ii) >= jq.size()) || (static_cast<Size>(jj) >= jq.size()),
        "The block indices must either be both -1 (default) or\n"
        "non-negative and smaller than the number of retrieval \n"
        "quantities.");
    ARTS_USER_ERROR_IF (ii > jj,
        "Only blocks above the diagonal can be set, hence"
        "*i* must be less than or equal to *j*.");
  }

  Index m = block_inv.nrows();
  Index n = block_inv.ncols();

  Index jq_m = jq[ii].size();
  if (jq[ii].HasAffine()) {
    jq_m = jq[ii].TransformationMatrix().ncols();
  }

  Index jq_n = jq[jj].size();
  if (jq[jj].HasAffine()) {
    jq_n = jq[jj].TransformationMatrix().ncols();
  }

  ARTS_USER_ERROR_IF (!((m == jq_m) && (n == jq_n)),
      "The dimensions of the covariance block ( ", m,
      " x ", n, " )"
      " with the dimensionality of "
      " retrieval quantity ", ii, " and ", jj, ", respectively.")

  ARTS_USER_ERROR_IF (!covmat_sx.has_block(ii, jj),
        "To add the inverse of a block the non-inverse"
        " block must be added first.");

  ArrayOfArrayOfIndex ji;
  bool any_affine;
  jac_ranges_indices(ji, any_affine, jq);
  Index row_start = ji[ii][0];
  Index row_extent = ji[ii][1] - ji[ii][0] + 1;
  Range row_range(row_start, row_extent);

  Index col_start = ji[jj][0];
  Index col_extent = ji[jj][1] - ji[jj][0] + 1;
  Range col_range(col_start, col_extent);

  std::shared_ptr<MatrixType> mat = std::make_shared<MatrixType>(block_inv);
  covmat_sx.add_correlation_inverse(
      Block(row_range, col_range, std::make_pair(ii, jj), mat));
}

void covmat_sxAddInverseBlock(CovarianceMatrix& covmat_sx,
                              const ArrayOfRetrievalQuantity& jq,
                              const Matrix& block_inv,
                              const Index& i,
                              const Index& j) {
  covmat_sxAddInverseBlock<Matrix>(covmat_sx, jq, block_inv, i, j);
}

void covmat_sxAddInverseBlock(CovarianceMatrix& covmat_sx,
                              const ArrayOfRetrievalQuantity& jq,
                              const Sparse& block_inv,
                              const Index& i,
                              const Index& j) {
  covmat_sxAddInverseBlock<Sparse>(covmat_sx, jq, block_inv, i, j);
}

void covmat_sxExtractSqrtDiagonal(Vector& x_norm,
                                  const CovarianceMatrix& covmat_sx) {
  x_norm = covmat_sx.diagonal();
  for (Index i = 0; i < x_norm.size(); ++i) {
    x_norm[i] = sqrt(x_norm[i]);
  }
}

////////////////////////////////////////////////////////////////////////////////
// retrievalAdd Functions
////////////////////////////////////////////////////////////////////////////////

void retrievalAddAbsSpecies(const Workspace& ws,
                            CovarianceMatrix& covmat_sx,
                            ArrayOfRetrievalQuantity& jacobian_quantities,
                            Agenda& jacobian_agenda,
                            const Sparse& covmat_block,
                            const Sparse& covmat_inv_block,
                            const Vector& rq_p_grid,
                            const Vector& rq_lat_grid,
                            const Vector& rq_lon_grid,
                            const String& species,
                            const String& mode,
                            const Index &for_species_tag) {
  jacobianAddAbsSpecies(ws,
                        jacobian_quantities,
                        jacobian_agenda,
                        rq_p_grid,
                        rq_lat_grid,
                        rq_lon_grid,
                        species,
                        mode,
                        for_species_tag);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddFreqShift(const Workspace& ws,
                           CovarianceMatrix& covmat_sx,
                           ArrayOfRetrievalQuantity& jacobian_quantities,
                           Agenda& jacobian_agenda,
                           const Sparse& covmat_block,
                           const Sparse& covmat_inv_block,
                           const Vector& f_grid,
                           const Numeric &df) {
  jacobianAddFreqShift(
      ws, jacobian_quantities, jacobian_agenda, f_grid, df);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      1,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddFreqStretch(const Workspace& ws,
                             CovarianceMatrix& covmat_sx,
                             ArrayOfRetrievalQuantity& jacobian_quantities,
                             Agenda& jacobian_agenda,
                             const Vector& f_grid,
                             const Sparse& covmat_block,
                             const Sparse& covmat_inv_block,
                             const Numeric&df) {
  jacobianAddFreqStretch(
      ws, jacobian_quantities, jacobian_agenda, f_grid, df);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      1,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddCatalogParameter(const Workspace& ws,
                                  CovarianceMatrix& covmat_sx,
                                  ArrayOfRetrievalQuantity& jacobian_quantities,
                                  Agenda& jacobian_agenda,
                                  const QuantumIdentifier& catalog_identity,
                                  const String& catalog_parameter,
                                  const Numeric &var) {
  jacobianAddBasicCatalogParameter(ws,
                                   jacobian_quantities,
                                   jacobian_agenda,
                                   catalog_identity,
                                   catalog_parameter);
  add_scalar_variance(covmat_sx, jacobian_quantities, var);
}

void retrievalAddCatalogParameters(
    const Workspace& ws,
    CovarianceMatrix& covmat_sx,
    ArrayOfRetrievalQuantity& jacobian_quantities,
    Agenda& jacobian_agenda,
    const Sparse& covmat_block,
    const Sparse& covmat_inv_block,
    const ArrayOfQuantumIdentifier& catalog_identities,
    const ArrayOfString &catalog_parameters) {
  jacobianAddBasicCatalogParameters(ws,
                                    jacobian_quantities,
                                    jacobian_agenda,
                                    catalog_identities,
                                    catalog_parameters);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      0,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddMagField(const Workspace& ws,
                          CovarianceMatrix& covmat_sx,
                          ArrayOfRetrievalQuantity& jacobian_quantities,
                          Agenda& jacobian_agenda,
                          const Sparse& covmat_block,
                          const Sparse& covmat_inv_block,
                          const Vector& rq_p_grid,
                          const Vector& rq_lat_grid,
                          const Vector& rq_lon_grid,
                          const String& component,
                          const Numeric &dB) {
  jacobianAddMagField(ws,
                      jacobian_quantities,
                      jacobian_agenda,
                      rq_p_grid,
                      rq_lat_grid,
                      rq_lon_grid,
                      component,
                      dB);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddPointingZa(const Workspace& ws,
                            CovarianceMatrix& covmat_sx,
                            ArrayOfRetrievalQuantity& jacobian_quantities,
                            Agenda& jacobian_agenda,
                            const Sparse& covmat_block,
                            const Sparse& covmat_inv_block,
                            const Matrix& sensor_pos,
                            const ArrayOfTime& sensor_time,
                            const Index& poly_order,
                            const String& calcmode,
                            const Numeric &dza) {
  jacobianAddPointingZa(ws,
                        jacobian_quantities,
                        jacobian_agenda,
                        sensor_pos,
                        sensor_time,
                        poly_order,
                        calcmode,
                        dza);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      1,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddPolyfit(const Workspace& ws,
                         CovarianceMatrix& covmat_sx,
                         ArrayOfRetrievalQuantity& jacobian_quantities,
                         Agenda& jacobian_agenda,
                         const Sparse& covmat_block,
                         const Sparse& covmat_inv_block,
                         const ArrayOfIndex& sensor_response_pol_grid,
                         const Matrix& sensor_response_dlos_grid,
                         const Matrix& sensor_pos,
                         const Index& poly_order,
                         const Index& no_pol_variation,
                         const Index& no_los_variation,
                         const Index &no_mblock_variation) {
  size_t jq_start = jacobian_quantities.size();
  jacobianAddPolyfit(ws,
                     jacobian_quantities,
                     jacobian_agenda,
                     sensor_response_pol_grid,
                     sensor_response_dlos_grid,
                     sensor_pos,
                     poly_order,
                     no_pol_variation,
                     no_los_variation,
                     no_mblock_variation);
  for (Index i = 0; i <= poly_order; ++i) {
    check_and_add_block(covmat_sx,
                        jacobian_quantities[jq_start + i],
                        jq_start + i,
                        4,
                        covmat_block,
                        covmat_inv_block);
  }
}

void retrievalAddScatSpecies(const Workspace& ws,
                             CovarianceMatrix& covmat_sx,
                             ArrayOfRetrievalQuantity& jacobian_quantities,
                             Agenda& jacobian_agenda,
                             const Sparse& covmat_block,
                             const Sparse& covmat_inv_block,
                             const Vector& rq_p_grid,
                             const Vector& rq_lat_grid,
                             const Vector& rq_lon_grid,
                             const String& species,
                             const String& quantity) {
  jacobianAddScatSpecies(ws,
                         jacobian_quantities,
                         jacobian_agenda,
                         rq_p_grid,
                         rq_lat_grid,
                         rq_lon_grid,
                         species,
                         quantity);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddSinefit(const Workspace& ws,
                         CovarianceMatrix& covmat_sx,
                         ArrayOfRetrievalQuantity& jacobian_quantities,
                         Agenda& jacobian_agenda,
                         const Sparse& covmat_block,
                         const Sparse& covmat_inv_block,
                         const ArrayOfIndex& sensor_response_pol_grid,
                         const Matrix& sensor_response_dlos_grid,
                         const Matrix& sensor_pos,
                         const Vector& period_lengths,
                         const Index& no_pol_variation,
                         const Index& no_los_variation,
                         const Index &no_mblock_variation) {
  size_t jq_start = jacobian_quantities.size();
  jacobianAddSinefit(ws,
                     jacobian_quantities,
                     jacobian_agenda,
                     sensor_response_pol_grid,
                     sensor_response_dlos_grid,
                     sensor_pos,
                     period_lengths,
                     no_pol_variation,
                     no_los_variation,
                     no_mblock_variation);
  for (Index i = 0; i < period_lengths.size(); ++i) {
    check_and_add_block(covmat_sx,
                        jacobian_quantities[jq_start + i],
                        jq_start + i,
                        4,
                        covmat_block,
                        covmat_inv_block);
  }
}

void retrievalAddSpecialSpecies(const Workspace& ws,
                                CovarianceMatrix& covmat_sx,
                                ArrayOfRetrievalQuantity& jacobian_quantities,
                                Agenda& jacobian_agenda,
                                const Sparse& covmat_block,
                                const Sparse& covmat_inv_block,
                                const Vector& rq_p_grid,
                                const Vector& rq_lat_grid,
                                const Vector& rq_lon_grid,
                                const String &species) {
  jacobianAddSpecialSpecies(ws,
                            jacobian_quantities,
                            jacobian_agenda,
                            rq_p_grid,
                            rq_lat_grid,
                            rq_lon_grid,
                            species);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddWind(const Workspace& ws,
                      CovarianceMatrix& covmat_sx,
                      ArrayOfRetrievalQuantity& jacobian_quantities,
                      Agenda& jacobian_agenda,
                      const Sparse& covmat_block,
                      const Sparse& covmat_inv_block,
                      const Vector& rq_p_grid,
                      const Vector& rq_lat_grid,
                      const Vector& rq_lon_grid,
                      const String& component,
                      const Numeric &dfrequency) {
  jacobianAddWind(ws,
                  jacobian_quantities,
                  jacobian_agenda,
                  rq_p_grid,
                  rq_lat_grid,
                  rq_lon_grid,
                  component,
                  dfrequency);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddTemperature(const Workspace& ws,
                             CovarianceMatrix& covmat_sx,
                             ArrayOfRetrievalQuantity& jacobian_quantities,
                             Agenda& jacobian_agenda,
                             const Sparse& covmat_block,
                             const Sparse& covmat_inv_block,
                             const Vector& rq_p_grid,
                             const Vector& rq_lat_grid,
                             const Vector& rq_lon_grid,
                             const String& hse) {
  jacobianAddTemperature(ws,
                         jacobian_quantities,
                         jacobian_agenda,
                         rq_p_grid,
                         rq_lat_grid,
                         rq_lon_grid,
                         hse);
  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3,
                      covmat_block,
                      covmat_inv_block);
}

void retrievalAddSurfaceQuantity(const Workspace& ws,
                                 CovarianceMatrix& covmat_sx,
                                 ArrayOfRetrievalQuantity& jacobian_quantities,
                                 Agenda& jacobian_agenda,
                                 const Sparse& covmat_block,
                                 const Sparse& covmat_inv_block,
                                 const Vector& rq_lat_grid,
                                 const Vector& rq_lon_grid,
                                 const String &quantity) {
  jacobianAddSurfaceQuantity(ws,
                             jacobian_quantities,
                             jacobian_agenda,
                             rq_lat_grid,
                             rq_lon_grid,
                             quantity);

  check_and_add_block(covmat_sx,
                      jacobian_quantities.back(),
                      jacobian_quantities.size() - 1,
                      3 - 1,
                      covmat_block,
                      covmat_inv_block);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void retrievalDefClose(Index& jacobian_do,
                       Agenda& jacobian_agenda,
                       Index& retrieval_checked,
                       const CovarianceMatrix& covmat_sx,
                       const ArrayOfRetrievalQuantity &jacobian_quantities) {
  jacobianClose(jacobian_do, jacobian_agenda, jacobian_quantities);

  ArrayOfArrayOfIndex ji_t;
  bool any_affine;
  jac_ranges_indices(ji_t, any_affine, jacobian_quantities);
  ARTS_USER_ERROR_IF (!covmat_sx.has_diagonal_blocks(ji_t),
      "*covmat_sx* does not contain a diagonal block for each retrieval"
      " quantity in the Jacobian.\n"
      " Fails test (!covmat_sx.has_diagonal_blocks(ji_t)) for ji_t ", ji_t,
      "\n")
  ARTS_USER_ERROR_IF (!covmat_sx.is_consistent(ji_t),
      "The blocks in *covmat_sx* are not consistent with the retrieval"
      " quantities in the Jacobian.");
  retrieval_checked = true;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void retrievalDefInit(CovarianceMatrix& covmat_se,
                      CovarianceMatrix& covmat_sx,
                      Sparse& covmat_block,
                      Sparse& covmat_inv_block,
                      ArrayOfRetrievalQuantity& jacobian_quantities,
                      Agenda& jacobian_agenda,
                      const Index &initialize_jacobian) {
  if (initialize_jacobian == 1) {
    jacobianInit(jacobian_quantities, jacobian_agenda);
  }

  covmat_block = Sparse();
  covmat_inv_block = Sparse();
  covmat_sx = CovarianceMatrix();
  covmat_se = CovarianceMatrix();
}

void retrievalErrorsExtract(Vector& retrieval_eo,
                            Vector& retrieval_ss,
                            const Matrix& covmat_so,
                            const Matrix&covmat_ss ) {
  Index n_so = covmat_so.nrows();
  Index n_ss = covmat_ss.nrows();

  retrieval_eo.resize(n_so);
  for (Index i = 0; i < n_so; ++i) {
    retrieval_eo[i] = sqrt(covmat_so(i, i));
  }

  retrieval_ss.resize(n_ss);
  for (Index i = 0; i < n_ss; ++i) {
    retrieval_ss[i] = sqrt(covmat_ss(i, i));
  }
}
