/*!
  \file   test_sparse.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:10:44 2003

  \brief  Tests for sparse matrices.

  Add more tests here as necessary...
*/

#include <iostream>
#include <stdexcept>
#include "lin_alg.h"
#include "matpack_data.h"
#include "matpack_sparse.h"
#include "test_utils.h"
#include "xml_io.h"

void test3() {
  Sparse M(10, 15);

  /*
  std::cout << "M.nrows(), M.ncols() = "
       << M.nrows() << ", " << M.ncols() << "\n";
  */
  for (Index i = 3; i < 10; ++i) M.rw(i, i) = (Numeric)i + 1;
  M.rw(0, 0) = 1;
  M.rw(0, 1) = 2;
  M.rw(0, 2) = 3;

  std::cout << "\nM = \n" << M;

  /*
  // Test Sparse matrix-Matrix multiplication
  Matrix A(10,5);
  Matrix C(15,5,2.0);
  // C = 2;

  mult(A, M, C(joker, joker));
  std::cout << "\nA = \n" << A << "\n";

  */

  // Test Sparse-Sparse multiplication
  //   Sparse A(10,5);
  //   Sparse C(15,5);
  //   for (Index i=0; i<5; i++) {
  //     C.rw(i*3,i) = i*3+1;
  //     C.rw(i*3+1,i) = i*3+2;
  //     C.rw(i*3+2,i) = i*3+3;
  //   }

  //   mult(A,M,C);

  //   std::cout << "\nA = \n" << A;

  /*
  // Test transpose
  Sparse B(15,10);

  transpose(B,M);
  std::cout << "\nM' = \n" << B;
  */

  /*
  // Test rw-operator
  Sparse S(M);
  S.rw(2,0) = 5;

  std::cout << "\nS(2,0) = " << S.ro(2,0) << "\n";

  std::cout << "\nS = \n" << S;
  */

  /*
  // Test vector multiplication
  Vector y(20, 0.0);
  Vector x(1,30,1);

  mult(y[Range(1,10,2)], S, x[Range(1,15,2)]);

  std::cout << "\ny = \n" << y << "\n";
  */
}

// void test38()
// {
//   std::cout << "Test sparse matrix - sparse matrix multiplication\n";

//   Sparse B(757,2271);
//   Sparse C(B.ncols(),B.ncols());
//   Sparse A(B.nrows(),B.ncols());

//   Index i=0;
//   for (Index j=0; j<B.ncols(); j++) {
//     B.rw(i,j) = (j+1.0);
//     if( i<B.nrows()-1 )
//       i++;
//     else
//       i=0;
//   }

//   for (Index i=0; i<C.nrows(); i++)
//     C.rw(i,i) = 1;

//   mult(A,B,C);

//   std::cout << "\n(Sparse) A = \n" << A;
//   std::cout << "\n(Sparse) B = \n" << B;
//   std::cout << "\n(Sparse) C = \n" << C;

//   Matrix a(5,15), b(5,15), c(15,15);

//   i=0;
//   for (Index j=0; j<15; j++) {
//     b(i,j) = j+1;
//     if( i<4 )
//       i++;
//     else
//       i=0;
//   }

//   for (Index i=0; i<15; i++)
//     c(i,i) = 1;

//   mult(a,b,c);

//   //std::cout << "\n(Full) a = \n" << a << "\n";

// //  std::cout << "\n(Sparse) B = \n" << B << "\n";
// //  std::cout << "\n(Full) b = \n" << b << "\n";
// //  std::cout << "\n(Sparse) C = \n" << C << "\n";
// //  std::cout << "\n(Full) c = \n" << c << "\n";

// }

// void test39()
// {
//   //Test sparse transpose function
//   Sparse B(1000,2000);
//   Sparse Bt(B.ncols(), B.nrows());

//   Index i=0;
//   for (Index j=0; j<B.ncols(); j=j+2) {
//     B.rw(i,j) = j+1;
//     if( i<B.nrows()-2 )
//       i += 2;
//     else
//       i=0;
//   }

//   std::cout << "\nB = \n" << B;

//   transpose(Bt,B);

//   std::cout << "\ntranspose(B) = \n" << Bt << "\n";
// }

void test40() {
  std::cout << "Testing the new simplified Sparse matrices:\n";

  Sparse A(3, 3);
  std::cout << "Empty A: " << A << "\n";

  A.rw(0, 0) = 11;
  A.rw(1, 1) = 22;
  A.rw(2, 2) = 33;
  std::cout << "Diagonal A:\n" << A << "\n";

  Vector b=uniform_grid(1, 3, 1), c(3);
  std::cout << "b:\n" << b << "\n";

  mult(c, A, b);
  std::cout << "c = A*b (should be [11,44,99]):\n" << c << "\n";

  Matrix D(3, 2);
  D[joker, 0] = b;
  D[joker, 1] = b;
  D[joker, 1] *= 2;
  std::cout << "D:\n" << D << "\n";

  Matrix E(3, 2);
  mult(E, A, D);
  std::cout << "E = A*D (should be [11,22],[44,88],[99,198]):\n" << E << "\n";
}

void test41() {
  std::cout << "Testing transpose for the new simplified sparse matrices:\n";

  Sparse B(4, 5);
  Index r[] = {0, 1, 1, 2, 2, 2, 3, 1, 3};
  Index c[] = {0, 0, 1, 1, 2, 3, 3, 4, 4};
  for (Index i = 0; i < 9; i++) B.rw(r[i], c[i]) = (Numeric)(i + 1);

  std::cout << "B:\n" << B << "\n";

  Sparse A(5, 4);

  transpose(A, B);

  std::cout << "A:\n" << A << "\n";

  std::cout << "Testing with a fully occupied matrix:\n";

  for (Index ri = 0; ri < 4; ri++)
    for (Index ci = 0; ci < 5; ci++) {
      B.rw(ri, ci) = (Numeric)(ri * 10 + ci);
    }

  std::cout << "B:\n" << B << "\n";
  transpose(A, B);
  std::cout << "A:\n" << A << "\n";
}

void test42() {
  std::cout << "Testing sparse-sparse matrix multiplication:\n";

  Sparse B(4, 5);
  Index r[] = {0, 1, 1, 2, 2, 2, 3, 1, 3};
  Index c[] = {0, 0, 1, 1, 2, 3, 3, 4, 4};
  for (Index i = 0; i < 9; i++) B.rw(r[i], c[i]) = (Numeric)(i + 1);

  Sparse A(4, 4), Bt(5, 4);
  transpose(Bt, B);
  mult(A, B, Bt);

  std::cout << "A:\n" << A << "\n";
}

void test43() {
  std::cout << "Testing sparse copying:\n";

  Sparse B(4, 5);
  Index r[] = {0, 1, 1, 2, 2, 2, 3, 1, 3};
  Index c[] = {0, 0, 1, 1, 2, 3, 3, 4, 4};
  for (Index i = 0; i < 9; i++) B.rw(r[i], c[i]) = (Numeric)(i + 1);

  std::cout << "B:\n" << B << "\n";

  Sparse A;

  A = B;

  std::cout << "A:\n" << A << "\n";

  for (Index i = 0; i < 100; ++i) {
    B.rw(0, 0) += 1;
    A = B;
  }

  std::cout << "A now:\n" << A << "\n";
}

void test44() {
  std::cout << "Test to insert row in sparse:\n";

  Vector v(5, 10);

  Sparse B(4, 5);
  Index r[] = {0, 1, 1, 2, 2, 2, 3, 1, 3};
  Index c[] = {0, 0, 1, 1, 2, 3, 3, 4, 4};
  for (Index i = 0; i < 9; i++) B.rw(r[i], c[i]) = (Numeric)(i + 1);

  std::cout << "B[" << B.nrows() << "," << B.ncols() << "]:\n" << B << "\n";
  std::cout << "v:\n" << v << "\n";

  B.insert_row(3, v);

  std::cout << "B (after insertion):\n" << B << "\n";
}

void test45() {
  std::cout << "Test Sparse-Sparse multiplication reading matrices from xml "
          "files:\n";

  Sparse A, B;
  String a = "antenna.xml";
  String b = "backend.xml";

  try {
    std::cout << "  Reading " << a << "...";
    xml_read_from_file(a, A);
    std::cout << "done.\n  Reading " << b << "...";
    xml_read_from_file(b, B);
    std::cout << "done.\n";
  } catch (const std::runtime_error &e) {
    std::cerr << e.what() << std::endl;
  }

  Sparse C(B.nrows(), A.ncols());
  std::cout << "  Performing multiplication...";
  mult(C, B, A);
  std::cout << "done.\n";

  //std::cout << "C=A*B:\n" << A << "\n";
  try {
    std::cout << "  Writing product to file: test45.xml...";
    xml_write_to_file("test45.xml", C, FileType::ascii, 0);
    std::cout << "done.\n";
  } catch (const std::runtime_error &e) {
    std::cerr << e.what() << std::endl;
  }
}

void test46() {
  std::cout << "Test transpose with large matrix read from xml file:\n";

  Sparse A;
  String a = "backend.xml";

  try {
    std::cout << "  Reading " << a << "...";
    xml_read_from_file(a, A);
    std::cout << "done.\n";
  } catch (const std::runtime_error &e) {
    std::cerr << e.what() << std::endl;
  }

  //std::cout << "A:\n" << A << std::endl;

  Sparse B(A.ncols(), A.nrows());
  transpose(B, A);

  try {
    std::cout << "  Writing transpose(A) to file test46.xml" << std::endl;
    xml_write_to_file("test46.xml", B, FileType::ascii, 0);
  } catch (const std::runtime_error &e) {
    std::cerr << e.what() << std::endl;
  }

  //std::cout << "transpose(A):\n" << B << std::endl;
}

void test47() {
  std::cout << "Test make Identity matrix:\n";

  Sparse A;

  A.resize(5, 5);
  id_mat(A);

  std::cout << "A:\n" << A << std::endl;
}

void test48() {
  std::cout << "Test absolute values of sparse matrix:\n";

  Sparse B(4, 5);
  Index r[] = {0, 1, 1, 2, 2, 2, 3, 1, 3};
  Index c[] = {0, 0, 1, 1, 2, 3, 3, 4, 4};
  for (Index i = 0; i < 9; i++) B.rw(r[i], c[i]) = -(Numeric)i * 0.5;
  std::cout << "B:\n" << B << std::endl;

  Sparse A(B);
  abs(A, B);

  std::cout << "abs(B):\n" << A << std::endl;
}

void test49() {
  std::cout << "Testing sparse adding:\n";

  Sparse B(4, 5);
  Index rb[] = {1, 3};
  Index cb[] = {1, 3};
  for (Index i = 0; i < 2; i++) B.rw(rb[i], cb[i]) = (Numeric)(i + 1);

  Sparse C(4, 5);
  Index rc[] = {0, 1, 2};
  Index cc[] = {0, 1, 2};
  for (Index i = 0; i < 3; i++) {
    C.rw(rc[i], cc[i]) = (Numeric)(i + 1);
    std::cout << C.rw(rc[i], cc[i]) << std::endl;
  }

  std::cout << "B:\n" << B << "\n";
  std::cout << "C:\n" << C << "\n";

  Sparse A;
  add(A, B, C);
  std::cout << "A=B+C:\n" << A << "\n";
  Sparse D;
  sub(D, B, C);
  std::cout << "D=B-C:\n" << D << "\n";
}

Numeric test_xml_io(Index ntests, bool verbose) {
  if (verbose) std::cout << "Testing xml IO." << std::endl << std::endl;

  Index n, m;

  Numeric err_max = 0.0;

  for (int i = 0; i < ntests; i++) {
    m = rand() % 1000 + 1;
    n = rand() % 1000 + 1;

    Sparse A(m, n), B;
    String a("A.xml");
    random_fill_matrix(A, 10, false);
    xml_write_to_file(a, A, FileType::ascii, 0);
    xml_read_from_file(a, B);
    xml_write_to_file("B.xml", B, FileType::ascii, 0);

    Numeric err =
        get_maximum_error(static_cast<Matrix>(B), static_cast<Matrix>(A), true);
    if (err > err_max) err_max = err;

    if (verbose)
      std::cout << "Test " << i << ": "
           << "Max. Rel. Error = " << err << std::endl;
  }

  return err_max;
}

//! Test sparse insert_row function.
/*!

  Performs ntests randomized tests of the Sparse::insert_row(...) function.
  For each test a random vector is inserted as row r into a sparse matrix.
  The sparse matrix is then transformed into a dense matrix and the row r
  compared to the vector. Return the maximum error between the rth row in
  the sparse matrix and the inserted vector, which should be 0.

  \param ntests Number of test to perform.
  \param verbose If verbose == true, the error for each test is printed
  to stdout.

  \return The maximum error between the inserted vector and the corresponding
  row of the matrix.
*/
Numeric test_insert_row(Index ntests, bool verbose) {
  Numeric err_max = 0.0;

  Vector v;
  Matrix A;
  Sparse A_sparse;

  if (verbose) std::cout << std::endl << "Testing insert_row:" << std::endl << std::endl;
  ;

  for (Index i = 0; i < ntests; i++) {
    Index m = (std::rand() % 10) + 1;
    Index n = (std::rand() % 10) + 1;
    Index r = std::rand() % m;

    v.resize(n);
    random_fill_vector(v, 10, false);

    A.resize(m, n);

    // Test resizing as well.
    A_sparse.resize(m, n);
    if ((A_sparse.nrows() != m) || (A_sparse.ncols() != n)) {
      if (verbose) {
        std::cout << "FAILED: Resize failed." << std::endl;
        return 1.0;
      }
    }

    A_sparse.insert_row(r, v);

    A = static_cast<Matrix>(A_sparse);

    Numeric err = get_maximum_error(A[r, joker], v, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::endl;
      std::cout << "Maximum relative error: " << err << std::endl;
    }
  }

  return err_max;
}

//! Test sparse identity matrix.
/*!

  Test the creation of sparse identity matrices using Sparse::make_I(...).
  For each test, a randomly sized sparse matrix is created and set to the
  identity matrix. The result is compared to the result of the dense std::couterpart
  identity(...). ntests sets the number of tests to be performed. Also tests
  if the matrix was correctly resized. If the resize fails the function returns
  1.0.

  \param ntest The number of tests to perform.
  \param verbose If true, test results for each test are printed to stdout.

  \return The maximum relative error between the sparse identity matrix and
  the dense identity matrix taken over all tests. Return 1.0 if the resize
  operation fails.
*/
Numeric test_identity(Index ntests, bool verbose) {
  Numeric err_max = 0.0;

  Matrix A;
  Sparse A_sparse;

  if (verbose) std::cout << std::endl << "Testing sparse identity:" << std::endl << std::endl;

  for (Index i = 0; i < ntests; i++) {
    Index n = (std::rand() % 1000) + 1;

    A.resize(n, n);
    A_sparse.resize(n, n);

    if ((A_sparse.nrows() != n) || (A_sparse.ncols() != n)) {
      if (verbose) {
        std::cout << std::endl << "FAILED: Resize failed." << std::endl;
        return 1.0;
      }
    }

    id_mat(A);
    id_mat(A_sparse);

    Numeric err = get_maximum_error(static_cast<Matrix>(A_sparse), A, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << "\t Maximum relative error: " << err << std::endl;
    }
  }

  return err_max;
}

//! Test sparse matrix construction.
/*!

  Performs ntests of the construction and conversion of sparse matrices to
  dense matrices. In each test a sparse m-times-n matrix is created and
  filled with random values. The matrix is then converted to a dense matrix
  and the error between the two matrices is computed.

  \param m The number of rows of the sparse matrix.
  \param n The number of columns of the sparse matrix.
  \param ntests The number of test to be performed.
  \param verbose If true, the results of each test are printed to stdout.

  \return Returns the maximum relative error between the sparse matrix and the
  converted, dense matrix taken over all ntests tests.
*/
Numeric test_sparse_construction(Index m, Index n, Index ntests, bool verbose) {
  Numeric err_max = 0.0;

  if (verbose) {
    std::cout << std::endl;
    std::cout << "Testing sparse matrix construction and conversion:" << std::endl;
  }

  for (Index i = 0; i < ntests; i++) {
    Matrix A;
    Sparse A_sparse(m, n);

    random_fill_matrix(A_sparse, 10, false);
    A = static_cast<Matrix>(A_sparse);

    Numeric err = get_maximum_error(static_cast<Matrix>(A_sparse), A, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << "Test No.: " << std::setw(5) << i;
      std::cout << " Maximum rel. error: " << err_max << std::endl;
    }
  }
  return err_max;
}

//! Test unary operation on sparse matrices.
/*!

  Perform ntests tests of the unary operation on sparse matrices
  (abs(...), transpose(...)) as well as construction of sparse matrices
  using rw, ro. The test of the construction of sparse matrix is done
  by filling a dense matrix with identical random values at the same positions
  and comparing the dense and the sparse matrices.

  \param m The number of rows of the sparse matrix.
  \param n The number of colums of the sparse matrix.
  \param ntests The number of tests to be performed.
  \param verbose If true, results for each test are printed to stdout.

  \return The maximum relative error between the sparse matrix and the
  dense counterpart taken over all tests.
*/
Numeric test_sparse_unary_operations(Index m,
                                     Index n,
                                     Index ntests,
                                     bool verbose) {
  Numeric err_max = 0.0;

  if (verbose) {
    std::cout << std::endl;
    std::cout << "Testing sparse unary operations:" << std::endl << std::endl;
    std::cout << std::setw(5) << "Test " << std::setw(15) << "Construction";
    std::cout << std::setw(15) << "abs" << std::setw(15) << "transpose" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
  }

  for (Index i = 0; i < ntests; i++) {
    Matrix A(m, n), B(m, n);
    Sparse A_sparse(m, n), B_sparse(m, n), C_sparse(n, m);

    A = 0.0;
    random_fill_matrix(A, A_sparse, 10, false);

    // Construction

    Numeric err = get_maximum_error(static_cast<Matrix>(A_sparse), A, true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(5) << i << std::setw(15) << err;

    // abs

    for (Index j = 0; j < m; j++) {
      for (Index k = 0; k < n; k++) {
        B[j, k] = abs(A[j, k]);
      }
    }

    abs(B_sparse, A_sparse);

    err = get_maximum_error(static_cast<Matrix>(B_sparse), B, true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(15) << err;

    // transpose

    transpose(C_sparse, A_sparse);

    err = get_maximum_error(static_cast<Matrix>(C_sparse), transpose(A), true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(15) << err << std::endl;
  }
  return err_max;
}

//! Test dense-sparse multiplication.
/*!

  Test multiplication B x C of a dense matrix B with a sparse matrix C. Takes a
  random submatrix of a m-times-n matrix B and multiplies it with a
  randomly generated sparse matrix C. The operation is simultaneously performed
  using dense arithmetic and the results are compared. Also tests the
  multiplication of transposed matrices.

  \param m The number of rows of the dense matrix C.
  \param n The number of cols ot the sparse matrix B.
  \param ntests The number of tests to be performed.
  \param verbose If true, the test results of each test are printed to stdout.

  \return Returns the maximum relative error between the sparse and the dense
  operation taken over all tests.
*/
Numeric test_dense_sparse_multiplication(Index m,
                                         Index n,
                                         Index ntests,
                                         bool verbose) {
  Numeric err_max = 0.0;

  if (verbose) {
    std::cout << std::endl;
    std::cout << "Testing sparse multiplication:" << std::endl << std::endl;
    std::cout << std::setw(5) << "Test " << std::setw(5) << "m1" << std::setw(5) << "m2";
    std::cout << std::setw(10) << "c_stride" << std::setw(5) << "n1" << std::setw(5) << "n2";
    std::cout << std::setw(10) << std::setw(10) << "r_stride" << std::setw(15) << "A = B x C";
    std::cout << std::setw(18) << "A^T = C^T x B^T" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
  }

  Matrix A(m, n);
  Matrix A_ref(m, n);
  Matrix B(m, n);
  random_fill_matrix(B, 10.0, true);

  for (Index i = 0; i < ntests; i++) {
    //
    // Generate random sub-matrix.
    //

    Index m1, dm, n1, dn, c_stride, r_stride;

    m1 = rand() % (m - 1);
    n1 = rand() % (n - 1);
    dm = (rand() % (m - m1)) + 1;
    dn = (rand() % (n - n1)) + 1;
    c_stride = (rand() % 10) + 1;
    r_stride = (rand() % 10) + 1;

    c_stride = std::min(dn, c_stride);
    r_stride = std::min(dm, r_stride);
    dn = dn / c_stride;
    dm = dm / r_stride;

    if (verbose) {
      std::cout << std::setw(5) << i << std::setw(5) << m1 << std::setw(5) << dm;
      std::cout << std::setw(10) << r_stride << std::setw(5) << n1 << std::setw(5) << dn;
      std::cout << std::setw(10) << c_stride;
    }

    MatrixView A_view = A[Range(m1, dm, r_stride), joker];
    MatrixView A_ref_view = A_ref[Range(m1, dm, r_stride), joker];
    Matrix C(dn, n);
    MatrixView B_mul = B[Range(m1, dm, r_stride), Range(n1, dn, c_stride)];
    Sparse C_sparse(dn, n), C_sparse_transpose(dn, n);

    //
    // Test standard dense-sparse multiplication.
    //

    random_fill_matrix(C_sparse, 10, false);
    C = static_cast<Matrix>(C_sparse);

    mult(A_ref_view, B_mul, C);
    mult(A_view, B_mul, C_sparse);

    // Compare results
    Numeric err = get_maximum_error(A_view, A_ref_view, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err;
    }

    //
    // Test transposed dense-sparse multiplication.
    //

    A.resize(n, m);
    A_ref.resize(n, m);
    MatrixView A_view_transp = A[joker, Range(m1, dm, r_stride)];
    MatrixView A_ref_view_transp = A_ref[joker, Range(m1, dm, r_stride)];

    C_sparse.resize(n, dn);
    C.resize(n, dn);
    random_fill_matrix(C_sparse, 10, false);
    transpose(C_sparse_transpose, C_sparse);
    C = static_cast<Matrix>(C_sparse);

    MatrixView B_mul_transp =
        B[Range(n1, dn, c_stride), Range(m1, dm, r_stride)];

    mult(transpose(A_ref_view_transp), transpose(B_mul_transp), transpose(C));
    mult(transpose(A_view_transp), transpose(B_mul_transp), C_sparse_transpose);

    // Compare results

    err = get_maximum_error(A_view, A_ref_view, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err << std::endl;
    }
  }
  return err_max;
}

//! Test sparse-dense multiplication.
/*!

  Test multiplication B x C of a sparse matrix B with a dense matrix C. Takes a
  random submatrix of a m-times-n matrix C and multiplies it with a
  randomly generated sparse matrix B. The operation is simultaneously performed
  using dense arithmetic and the results are compared. Also test the
  multiplication of transposed matrices.

  \param m The number of rows of the dense matrix C.
  \param n The number of cols ot the sparse matrix B.
  \param ntests The number of tests to be performed.
  \param verbose If true, the test results of each test are printed to stdout.

  \return Returns the maximum relative error between the sparse and the dense
  operation taken over all tests.
*/
Numeric test_sparse_dense_multiplication(Index m,
                                         Index n,
                                         Index ntests,
                                         bool verbose) {
  Numeric err_max = 0.0;

  if (verbose) {
    std::cout << std::endl;
    std::cout << "Testing sparse multiplication:" << std::endl << std::endl;
    std::cout << std::setw(5) << "Test " << std::setw(5) << "m1" << std::setw(5) << "m2";
    std::cout << std::setw(10) << "c_stride" << std::setw(5) << "n1" << std::setw(5) << "n2";
    std::cout << std::setw(10) << std::setw(10) << "r_stride" << std::setw(15) << "A = B x C";
    std::cout << std::setw(18) << "A^T = C^T x B^T" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
  }

  Matrix C(m, n);
  Matrix A(m, n);
  Matrix A_ref(m, n);
  random_fill_matrix(C, 10.0, true);

  for (Index i = 0; i < ntests; i++) {
    //
    // Generate random sub-matrix.
    //

    Index m1, dm, n1, dn, c_stride, r_stride;

    m1 = rand() % (m - 1);
    n1 = rand() % (n - 1);
    dm = (rand() % (m - m1)) + 1;
    dn = (rand() % (n - n1)) + 1;
    c_stride = (rand() % 10) + 1;
    r_stride = (rand() % 10) + 1;

    if (verbose) {
      std::cout << std::setw(5) << i << std::setw(5) << m1 << std::setw(5) << dm;
      std::cout << std::setw(10) << c_stride << std::setw(5) << n1 << std::setw(5) << dn;
      std::cout << std::setw(10) << r_stride;
    }

    c_stride = std::min(dn, c_stride);
    r_stride = std::min(dm, r_stride);
    dn = dn / c_stride;
    dm = dm / r_stride;

    Matrix B(m, dm);
    MatrixView A_view = A[joker, Range(n1, dn, c_stride)];
    MatrixView A_ref_view = A_ref[joker, Range(n1, dn, c_stride)];
    MatrixView C_mul = C[Range(m1, dm, r_stride), Range(n1, dn, c_stride)];
    Sparse B_sparse(m, dm), B_sparse_transpose(m, dm);

    //
    // Test standard multiplication.
    //

    random_fill_matrix(static_cast<Matrix>(B_sparse), 10, false);
    B = static_cast<Matrix>(B_sparse);

    // Sparse-sparse multiplication
    mult(A_ref_view, B, C_mul);
    mult(A_view, B_sparse, C_mul);

    // Compare results
    Numeric err = get_maximum_error(A_view, A_ref_view, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err;
    }

    //
    // Test transposed matrix views.
    //

    A.resize(n, m);
    A_ref.resize(n, m);
    MatrixView A_view_transp = A[Range(n1, dn, c_stride), joker];
    MatrixView A_ref_view_transp = A_ref[Range(n1, dn, c_stride), joker];

    B_sparse.resize(dm, m);
    B.resize(dm, m);
    random_fill_matrix(B_sparse, 10, false);
    transpose(B_sparse_transpose, B_sparse);
    B = static_cast<Matrix>(B_sparse);

    MatrixView C_mul2 = C[Range(n1, dn, c_stride), Range(m1, dm, r_stride)];

    // Transposed sparse-dense multiplication
    mult(transpose(A_ref_view_transp), transpose(B), transpose(C_mul2));
    mult(transpose(A_view_transp), B_sparse_transpose, transpose(C_mul2));

    // Compare results
    err = get_maximum_error(A_view, A_ref_view, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err << std::endl;
    }
  }
  return err_max;
}

//! Test sparse multiplication.
/*!

  Test multiplication of sparse matrices with sparse and dense matrices as well
  as vectors. Performs ntests test, where in each test the product of two sparse
  matrices, of a sparse and a dense matrix and of a sparse matrix and a vector
  are computed:

      Matrix-matrix product: A = B x C
      Matrix-vector product: y = B * x

  Where A is a m-times-n matrix and B a m-times-k matrix, where k is picked
  randomly for each test from the range [1,999].
  The  results are compared to the results obtained using dense
  arithmetic and the maximum relative error taken over all tests is returned.

  \param m Number of rows of A.
  \param n Number of columns of A.
  \param ntests Number of tests to perform.
  \param verbose If true, the results of each test are printed to stdout.

  \return The maximum relative error taken over all operations and number
  of tests performed.
*/
Numeric test_sparse_multiplication(Index m,
                                   Index n,
                                   Index ntests,
                                   bool verbose) {
  Numeric err_max = 0.0;

  if (verbose) {
    std::cout << std::endl;
    std::cout << "Testing sparse multiplication:" << std::endl << std::endl;
    std::cout << std::setw(5) << "Test " << std::setw(15) << "sparse-sparse";
    std::cout << std::setw(15) << "sparse-dense" << std::setw(15) << "matrix-vector";
    std::cout << "transposed matrix-vector";
    std::cout << std::endl << std::string(80, '-') << std::endl;
  }

  for (Index i = 0; i < ntests; i++) {
    Index k;

    k = (rand() % 1000) + 1;

    Matrix A(m, n), A2(m, n), B, C;
    Vector y(m), y2(m), yt(k), yt2(k), x(k), xt(m);
    Sparse A_sparse(m, n), B_sparse(m, k), C_sparse(k, n);

    random_fill_matrix(B_sparse, 10, false);
    random_fill_matrix(C_sparse, 10, false);
    random_fill_vector(x, 10, false);
    random_fill_vector(xt, 10, false);
    B = static_cast<Matrix>(B_sparse);
    C = static_cast<Matrix>(C_sparse);

    // Sparse-sparse
    mult(A, B, C);
    mult(A_sparse, B_sparse, C_sparse);

    Numeric err = get_maximum_error(static_cast<Matrix>(A_sparse), A, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(5) << i << std::setw(15) << err;
    }

    // Sparse-dense
    mult(A, B, C);
    mult(A2, B_sparse, C);

    err = get_maximum_error(A2, A, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err;
    }

    // Matrix-vector
    mult(y, B, x);
    mult(y2, B_sparse, x);

    err = get_maximum_error(y2, y, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err;
    }

    // transposed Matrix-vector
    mult(yt, transpose(B), xt);
    transpose_mult(yt2, B_sparse, xt);

    err = get_maximum_error(yt2, yt, true);
    if (err > err_max) err_max = err;

    if (verbose) {
      std::cout << std::setw(15) << err << std::endl;
    }
  }
  return err_max;
}

//! Test sparse matrix arithmetic.
/*!

  Test multiplication, addition and subtraction of sparse matrices using
  the corresponding dense operations. The operations performed are

      A x B, C + D, C - D

  for a m-times-n matrix C, a m-times-m matrix A, a m-times-n matrix B
  and a m-times-n matrix D.

  \param m The number of rows of C
  \param n The number of columns of C
  \param ntests The number of test to be performed
  \param verbose If true, the results of each test are printed to stdout.

  \return The maximum relative error taken over all tests and operations.
*/
Numeric test_sparse_arithmetic(Index m, Index n, Index ntests, bool verbose) {
  Numeric err_max = 0.0;

  if (verbose) {
    std::cout << std::endl;
    std::cout << "Testing sparse arithmetic:" << std::endl << std::endl;
    std::cout << std::setw(5) << "Test " << std::setw(15) << "Addition";
    std::cout << std::setw(15) << "Multiplication" << std::setw(15) << "Subtraction";
    std::cout << std::setw(15) << "Scaling" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
  }

  for (Index i = 0; i < ntests; i++) {
    Matrix A, B, C(m, n), D;
    Sparse A_sparse(m, m), B_sparse(m, n), C_sparse(m, n), D_sparse(m, n);

    random_fill_matrix(A_sparse, 10, false);
    random_fill_matrix(B_sparse, 10, false);
    random_fill_matrix(D_sparse, 10, false);
    A = static_cast<Matrix>(A_sparse);
    B = static_cast<Matrix>(B_sparse);
    D = static_cast<Matrix>(D_sparse);

    // Multiplication

    mult(C, A, B);
    mult(C_sparse, A_sparse, B_sparse);

    Numeric err = get_maximum_error(static_cast<Matrix>(C_sparse), C, true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(5) << i << std::setw(15) << err;

    // Addition

    C += D;
    add(B_sparse, C_sparse, D_sparse);

    err = get_maximum_error(static_cast<Matrix>(B_sparse), C, true);
    if (err > err_max) err_max = err;

    C_sparse += D_sparse;

    err = get_maximum_error(static_cast<Matrix>(C_sparse), C, true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(15) << err;

    // Subtraction

    C = static_cast<Matrix>(C_sparse);
    C -= D;
    sub(B_sparse, C_sparse, D_sparse);

    err = get_maximum_error(static_cast<Matrix>(B_sparse), C, true);
    if (err > err_max) err_max = err;

    C_sparse -= D_sparse;

    err = get_maximum_error(static_cast<Matrix>(C_sparse), C, true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(15) << err;

    // Scaling

    C_sparse *= 1.2;
    C *= 1.2;

    err = get_maximum_error(static_cast<Matrix>(C_sparse), C, true);
    if (err > err_max) err_max = err;

    C_sparse /= 1.2;
    C /= 1.2;

    err = get_maximum_error(static_cast<Matrix>(C_sparse), C, true);
    if (err > err_max) err_max = err;

    if (verbose) std::cout << std::setw(15) << err << std::endl;
  }
  return err_max;
}

int main() {
  // test3();
  // test38();
  // test39();
  // test40();
  // test41();
  // test42();
  // test43();
  // test44();
  // test45();
  // test46();
  // test47();
  // test48();

  Numeric err;

  std::cout << "Testing sparse arithmetic: ";
  err = test_sparse_arithmetic(1000, 1000, 100, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing xml IO: ";
  err = test_xml_io(100, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing identity matrix: ";
  err = test_identity(1000, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing inserting of rows: ";
  err = test_insert_row(1000, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing abs(...) and transpose(...): ";
  err = test_sparse_unary_operations(1000, 1000, 1000, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing multiplication: ";
  err = test_sparse_multiplication(1000, 1000, 1000, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing sparse-dense multiplication: ";
  err = test_sparse_dense_multiplication(1000, 1000, 1000, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  std::cout << "Testing dense-sparse multiplication: ";
  err = test_dense_sparse_multiplication(1000, 1000, 1000, false);
  if (err < 1e-11)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED (Error: " << err << ")" << std::endl;

  return 0;
}
