/* Copyright (C) 2017
   Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

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
  \file   test_covariance.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   Thu Jun 29 12:34:56 2017

  \brief  Tests for covariance matrices.
*/

#include <utility>
#include <cstdlib>
#include <random>
#include <limits>

#include "covariance_matrix.h"
#include "lin_alg.h"
#include "jacobian.h"
#include "array.h"
#include "test_utils.h"
#include "xml_io.h"

// Type alias describing a retrieval setup for testing.
using RetrievalData = std::tuple<ArrayOfRetrievalQuantity, ArrayOfArrayOfIndex>;

/**
 * Create a block of a covariance matrix representating a correlation between the
 * retrieval quantities i and j by applying the functional f to the grid positions
 * of each of the retrieval grids.
 *
 * @param i The index of the first retrieval quantity
 * @param j The index of the second retrieval quantity
 * @param jqs The array holding the retrieval quantities.
 * @param f The functional to apply to the correlated grid positions.
 *
 * @return Shared pointer to the Matrix object representing the correlation relation.
 */
template<typename F>
std::shared_ptr<Matrix> create_covariance_matrix_1D(
    Index i,
    Index j,
    ArrayOfRetrievalQuantity jqs,
    F f)
{
    if (i > j) {
        std::swap(i,j);
    }

    const RetrievalQuantity &q1 = jqs[i];
    const RetrievalQuantity &q2 = jqs[j];

    const Vector &gv1 = q1.Grids()[0];
    const Vector &gv2 = q2.Grids()[0];

    std::shared_ptr<Matrix> c = std::make_shared<Matrix>(gv1.nelem(), gv2.nelem());

    for (Index k = 0; k < gv1.nelem(); k++) {
        for (Index l = 0; l < gv2.nelem(); l++) {
            (*c)(k,l) = f(gv1[k], gv2[l]);
        }
    }
    return c;
}

/**
 * Same as create_covariance_matrix_1D but creates a matrix of type sparse.
 */
template<typename F>
std::shared_ptr<Sparse> create_sparse_covariance_matrix_1D(
    Index i,
    Index j,
    ArrayOfRetrievalQuantity jqs,
    F f)
{
    if (i > j) {
        std::swap(i,j);
    }

    const RetrievalQuantity &q1 = jqs[i];
    const RetrievalQuantity &q2 = jqs[j];

    const Vector &gv1 = q1.Grids()[0];
    const Vector &gv2 = q2.Grids()[0];

    std::shared_ptr<Sparse> s = std::make_shared<Sparse>(gv1.nelem(), gv2.nelem());

    Index nelem = gv1.nelem() * gv2.nelem();
    ArrayOfIndex row_indices{}, col_indices{};
    ArrayOfNumeric elements{};
    row_indices.reserve(nelem);
    col_indices.reserve(nelem);
    elements.reserve(nelem);
    for (Index k = 0; k < gv1.nelem(); k++) {
        for (Index l = 0; l < gv2.nelem(); l++) {
            row_indices.push_back(k);
            col_indices.push_back(l);
            elements.push_back(f(gv1[k], gv2[l]));
        }
    }
    s->insert_elements(nelem, row_indices, col_indices, Vector(elements));
    return s;
}

/**
 * Setup a random retrieval case for testing.
 *
 * @return A RetrievalData object holding an array of retrieval quantities and jacobian
 * indices.
 */
RetrievalData setup_retrieval_1D()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> n_rqs_dist(2, 10), n_gs_dist(10,100);

    Index n_rqs = n_rqs_dist(gen);
    Numeric g_begin = 0.0;
    Numeric g_end   = 100.0;
    ArrayOfRetrievalQuantity rqs = ArrayOfRetrievalQuantity();

    for (Index i = 0; i < n_rqs; i++) {
        Index n_gs = n_gs_dist(gen);
        ArrayOfVector gvs(1);
        gvs[0] = Vector(n_gs);
        for (Index j = 0; j < n_gs; j++) {
            Numeric frac = static_cast<Numeric>(j) / static_cast<Numeric>(n_gs);
            gvs[0][j] = g_begin + (g_end - g_begin) * frac;
        }
        rqs.emplace_back("maintag", "subtag", "subsubtag", "mode", 1, 0.0, gvs);
    }

    ArrayOfArrayOfIndex jis(n_rqs);
    Index ji = 0;
    for (Index i = 0; i < n_rqs; i++) {
        Index n_grid_points = 1;
        const ArrayOfVector & gs = rqs[i].Grids();
        for (auto & g : gs) {
            n_grid_points *= g.nelem();
        }
        jis[i] = ArrayOfIndex{ji, ji + n_grid_points - 1};
        ji += n_grid_points;
    }

    return std::make_tuple(rqs, jis);
}

/**
 * Create a covariance matrix with correlations between retrieval quantitites given
 * in rqs and jis.
 *
 * @param rqs The array holding the retrieval quantities.
 * @param jis The array holding the jacobian indices of the retrieval quantities.
 */
CovarianceMatrix random_covariance_matrix(
    const ArrayOfRetrievalQuantity &rqs,
    const ArrayOfArrayOfIndex      &jis
    )
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    CovarianceMatrix covmat{};

    auto g = [](Numeric x, Numeric y) -> Numeric {
        if (x==y) return 1.0; else return 0.0;};
    auto g2 = [](Numeric x, Numeric y) -> Numeric {
        if (x==y) return 0.1; else return 0.0;};

    for (size_t i = 0; i < rqs.size(); i++) {
        Range range(jis[i][0], jis[i][1] - jis[i][0] + 1);
        auto inds = std::make_pair(i, i);

        Numeric p = dis(gen);
        if (p < 0.0) {
            std::shared_ptr<Matrix> m = create_covariance_matrix_1D(i, i, rqs, g);
            covmat.add_correlation(Block(range, range, inds, m));
        } else {
            std::shared_ptr<Sparse> m = create_sparse_covariance_matrix_1D(i, i, rqs, g);
            covmat.add_correlation(Block(range, range, inds, m));
        }
    }

    for (size_t i = 0; i < rqs.size(); i++) {
        for (size_t j = i+1; j < rqs.size(); j++) {
            Range row_range(jis[i][0], jis[i][1] - jis[i][0] + 1);
            Range col_range(jis[j][0], jis[j][1] - jis[j][0] + 1);
            auto inds = std::make_pair(i, j);
            Numeric p = dis(gen);
            if (p > 0.8) {
                p = dis(gen);
                if (p < 0.0) {
                    std::shared_ptr<Matrix> m = create_covariance_matrix_1D(i, j, rqs, g);
                    covmat.add_correlation(Block(row_range, col_range, inds, m));
                } else {
                    std::shared_ptr<Sparse> m = create_sparse_covariance_matrix_1D(i, j, rqs, g2);
                    covmat.add_correlation(Block(row_range, col_range, inds, m));
                }
            }
        }
    }
    return covmat;
}

/**
 * Tests the multiplication of covariance matrices with vectors.
 *
 * @param  n_tests The number of tests to perform.
 * @return The maximum error of the result with respect to the same multiplication
 * performed using an identical matrix of type Matrix
 */
Numeric test_multiplication_by_vector(Index n_tests)
{
    Numeric e = 0.0;
    for (Index i = 0; i < n_tests; i++) {

        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat(random_covariance_matrix(rqs,jis));

        Index n = covmat.ncols();
        Matrix A(covmat);
        Vector v(n), w(n), w_ref(n);
        random_fill_vector(v, 10.0, false);

        mult(w, covmat, v);
        mult(w_ref, A, v);

        e = std::max(e,get_maximum_error(w, w_ref, true));
    }
    return e;
}


/**
 * Tests the multiplication of covariance matrices with matrices.
 *
 * @param  n_tests The number of tests to perform
 * @return The maximum error of the result with respect to the same multiplication
 * performed using an identical matrix of type Matrix
 */
Numeric test_multiplication_by_matrix(Index n_tests)
{
    Numeric e = 0.0;
    for (Index i = 0; i < n_tests; i++) {
        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat(random_covariance_matrix(rqs,jis));

        Index n = covmat.ncols();
        Matrix A(covmat), B(n,n), C(n,n), C_ref(n,n);
        random_fill_matrix(B, 10.0, false);

        mult(C,     covmat, B);
        mult(C_ref, A, B);
        e = std::max(e,get_maximum_error(C, C_ref, true));

        mult(C,B,covmat);
        mult(C_ref, B, A);
        e = std::max(e,get_maximum_error(C, C_ref, true));
    }
    return e;
}

/**
 * Tests the inversion of covariance matrices by computing the products of inverses
 * of covariance matrices with matrices and vectors and comparing the result to
 * what is obtained by performing these operations with an identical matrix of type
 * Matrix
 *
 * @param  n_tests The number of tests to perform
 * @return The maximum error of the result with respect to the same operations
 * performed using an identical matrix of type Matrix
 */
Numeric test_inverse(Index n_tests)
{
    Numeric e = 0.0;
    for (Index i = 0; i < n_tests; i++) {
        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat(random_covariance_matrix(rqs,jis));

        Index n = covmat.ncols();
        Matrix A(covmat), B(n,n), B_ref(n,n), C(n,n);
        random_fill_matrix(B, 10.0, false);

        covmat.compute_inverse();
        mult_inv(B, covmat, A);
        id_mat(B_ref);
        mult(C, covmat, B_ref);
        e = std::max(e,get_maximum_error(B, B_ref, true));
        id_mat(B_ref);
        mult_inv(B, A, covmat);
        e = std::max(e,get_maximum_error(B, B_ref, true));
    }
    return e;
}

/**
 * Test addition of covariance matrices and inverse covariance matrices.
 *
 * @param  n_tests The number of tests to perform
 * @return The maximum error of the result with respect to the same operations
 * performed using an identical matrix of type Matrix
 */
Numeric test_addition(Index n_tests)
{
    Numeric e = 0.0;
    for (Index i = 0; i < n_tests; i++) {
        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat(random_covariance_matrix(rqs,jis));

        Index n = covmat.ncols();
        Matrix A(covmat), B(n,n), B_ref(n,n), C(n,n);
        B = 0.0;
        B_ref = B;

        covmat.compute_inverse();

        B += covmat;
        B_ref += A;
        e = std::max(e, get_maximum_error(B, B_ref, true));

        add_inv(B, covmat);
        inv(C,A);
        B_ref += C;
        e = std::max(e, get_maximum_error(B, B_ref, true));
    }
    return e;
}

/**
 * Test extraction of diagonal.
 *
 *
 * @param  n_tests The number of tests to perform
 * @return The maximum error of the result with respect to the diagonal
 * extracted from an identical matrix of type Matrix
 */
Numeric test_diagonal(Index n_tests)
{
    Numeric e = 0.0;
    for (Index i = 0; i < n_tests; i++) {
        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat(random_covariance_matrix(rqs,jis));

        Vector diag_1 = covmat.diagonal();
        Vector diag_2 = Matrix(covmat).diagonal();
        e = std::max(e, get_maximum_error(diag_1, diag_2, true));
    }
    return e;
}

/**
 * Test input and output of covariance matrices.
 *
 * @param  n_tests The number of tests to perform
 * @return The maximum error of the original covariance matrix and a covariance that has
 * been stored to xml format and reread.
 */
Numeric test_io(Index n_tests)
{
    Numeric e(0.0);
    for (Index i = 0; i < n_tests; i++) {
        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat_1(random_covariance_matrix(rqs,jis)), covmat_2{};
        covmat_1.compute_inverse();

        xml_write_to_file("test.xml", covmat_1, FILE_TYPE_ASCII, 0, Verbosity());
        xml_read_from_file("test.xml", covmat_2, Verbosity());

        Matrix A(covmat_1), B(covmat_2), C(covmat_1);
        e = std::max(e, get_maximum_error(A, B, true));

        id_mat(C);
        mult_inv(A, covmat_1, C);
        mult_inv(B, covmat_2, C);
        e = std::max(e, get_maximum_error(A, B, true));
    }
    return 0.0;
}

template<typename MatrixType>
void covmat_seSet(CovarianceMatrix& covmat,
                  const MatrixType& block,
                  const Verbosity& /*v*/);

template<typename MatrixType>
void covmat_sxSet(CovarianceMatrix& covmat,
                  const MatrixType& block,
                  const Verbosity& /*v*/);

template<typename MatrixType>
void covmat_seAddBlock(CovarianceMatrix& covmat_se,
                       const MatrixType& block,
                       const Index& i,
                       const Index& j,
                       const Verbosity &);

template<typename MatrixType>
void covmat_seAddInverseBlock(CovarianceMatrix& covmat_se,
                       const MatrixType& block,
                       const Index& i,
                       const Index& j,
                       const Verbosity &);


template<typename MatrixType>
void covmat_sxAddBlock(CovarianceMatrix& covmat_se,
                       const ArrayOfRetrievalQuantity& jq,
                       const MatrixType& block,
                       const Index& i,
                       const Index& j,
                       const Verbosity &);

template<typename MatrixType>
void covmat_sxAddInverseBlock(CovarianceMatrix& covmat_se,
                       const ArrayOfRetrievalQuantity& jq,
                       const MatrixType& block,
                       const Index& i,
                       const Index& j,
                       const Verbosity &);
 
/**
 * Test general functionality of CovarianceMatrix class.
 */
void test_workspace_methods()
{
    CovarianceMatrix covmat{};
    Matrix A(10, 10);
    Sparse A_sparse(10, 10);
    Matrix B(20, 20);
    Matrix C(10, 20);

    // covmat_seSet

    covmat_seSet(covmat, A, Verbosity());
    assert(covmat.ncols() == 10);

    // covmat_seAddBlock

    covmat_seAddBlock(covmat, A_sparse, -1, -1, Verbosity());
    assert(covmat.ncols() == 20);

    try {
        covmat_seAddBlock(covmat, A, 3, 3, Verbosity());
        // This should fail.
        assert(false);
    } catch(std::runtime_error) {}

    covmat_seAddBlock(covmat, A, 0, 1, Verbosity());

    covmat_seAddBlock(covmat, A, 2, 2, Verbosity());

    try {
        covmat_seAddBlock(covmat, B, 1, 2, Verbosity());
        // This should fail.
        assert(false);
    } catch(std::runtime_error) {}

    covmat_seAddInverseBlock(covmat, A, 0, 1, Verbosity());

    try {
        covmat_seAddInverseBlock(covmat, B, 3, 3, Verbosity());
        // This should fail.
        assert(false);
    } catch(std::runtime_error) {}


    // covmat_sxSet
    covmat_sxSet(covmat, A, Verbosity());
    assert(covmat.ncols() == 10);

    // covmat_sxAddBlock

    A.resize(1000, 1000);
    A_sparse.resize(8000, 8000);
    C.resize(1000, 8000);

    ArrayOfVector grids_1{Vector(10), Vector(10), Vector(10)};
    ArrayOfVector grids_2{Vector(20), Vector(20), Vector(20)};

    RetrievalQuantity rq_1("mt", "st", "sst", "m", 1, 0.1, grids_1);
    RetrievalQuantity rq_2("mt", "st", "sst", "m", 1, 0.1, grids_2);

    ArrayOfRetrievalQuantity rqs{};
    rqs.push_back(rq_1);
    rqs.push_back(rq_2);

    covmat_sxAddBlock(covmat, rqs, A_sparse, -1, -1, Verbosity());

    try {
        covmat_sxAddBlock(covmat, rqs, A_sparse, 0, 1, Verbosity());
        // This should fail.
        assert(false);
    } catch(std::runtime_error) {}

    covmat_sxAddBlock(covmat, rqs, C, 0, 1, Verbosity());
    covmat_sxAddInverseBlock(covmat, rqs, A, 0, 0, Verbosity());

    try {
        covmat_sxAddInverseBlock(covmat, rqs, C, 1, 1, Verbosity());
        // This should fail.
        assert(false);
    } catch(std::runtime_error) {}
}

/**
 * Tests invlib wrapper for covariance matrices.
 *
 * @param  n_tests The number of tests to perform
 * @return The maximum error of the result with respect to the same multiplication
 * performed using a normal covariance matrix.
 */
namespace {
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/algebra.h"

Numeric test_invlib_wrapper(Index n_tests)
{
    Numeric e = 0.0;
    for (Index i = 0; i < n_tests; i++) {
        ArrayOfArrayOfIndex jis;
        ArrayOfRetrievalQuantity rqs;
        std::tie(rqs, jis) = setup_retrieval_1D();
        CovarianceMatrix covmat(random_covariance_matrix(rqs,jis));
        invlib::Matrix<ArtsCovarianceMatrixWrapper> wrapper(covmat);
        Index n = covmat.ncols();

        // Multiplication by Vector
        invlib::Vector<ArtsVector> v{}, w{}, w_ref{};
        v.resize(n);
        w.resize(n);
        w_ref.resize(n);
        random_fill_vector(v, 10.0, false);

        mult(w, covmat, v);
        w_ref = wrapper * v;
        e = std::max(e,get_maximum_error(w, w_ref, true));

        solve(w, covmat, v);
        w_ref = inv(wrapper) * v;
        e = std::max(e,get_maximum_error(w, w_ref, true));
        w_ref = inv(inv(wrapper)) * v;

        // Multiplication by Matrix
        invlib::Matrix<ArtsMatrix> A(covmat), B{}, C{}, C_ref{};
        B.resize(n,n);
        C.resize(n,n);
        C_ref.resize(n,n);

        random_fill_matrix(B, 10.0, false);
        mult(C, covmat, B);
        C_ref = wrapper * B;
        e = std::max(e,get_maximum_error(C, C_ref, true));

        // Inverse
        mult_inv(C, covmat, B);
        C_ref = invlib::inv(wrapper) * B;
        e = std::max(e,get_maximum_error(C, C_ref, true));

        // Addition to Matrix
        C += covmat;
        C_ref = C_ref + wrapper;
        e = std::max(e,get_maximum_error(C, C_ref, true));

        // Addition of Inverse to Matrix
        add_inv(C, covmat);
        C_ref = C_ref + inv(wrapper);
        e = std::max(e,get_maximum_error(C, C_ref, true));

    }
    return e;
}
}

int main() {
    Numeric e = 0.0, e_max = 0.0;

    std::cout << "Testing covariance matrix: " << std::endl;
    e = test_addition(10);
    std::cout << "\tAddition:                " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    e = test_multiplication_by_vector(10);
    std::cout << "\tVector Multiplication:   " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    e = test_multiplication_by_matrix(10);
    std::cout << "\tMatrix Multiplication:   " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    e = test_inverse(10);
    std::cout << "\tInverse:                 " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    e = test_io(10);
    std::cout << "\tXML IO:                  " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    e = test_invlib_wrapper(10);
    std::cout << "\tinvlib Wrapper:          " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    e = test_diagonal(10);
    std::cout << "\tdiagonal                 " << e << std::endl;
    e_max = std::max(e, e_max);
    if (e_max > 1e-5) {
        return -1;
    }

    std::cout << std::endl << "\tTesting workspace functions ... ";
    test_workspace_methods();
    std::cout << " DONE." << std::endl;

    return 0;
}
