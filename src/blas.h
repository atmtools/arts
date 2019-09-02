/*!
  \file   blas.h
  \author simon <simon@thinks>
  \date   Thu May  7 22:41:17 2015

  \brief  Interface for BLAS library.

  Contains declarations for BLAS functions.

*/

#include <complex>

//! BLAS matrix multiplication.
/*!
  Performs the operation
      C = alpha * (transa(A) * transb(B)) + beta * C
  where transa(A) is either A or A^T and trans(B) either B or B^T.

  Note the difference in memory layout between arts and BLAS. arts matrices
  using row-major order whereas BLAS uses column major order. Therefore, passing
  the raw-data of the arts matrices A, B, C to BLAS will compute

  C = (A^T * B^T)^T = B * A

  It is therefore necessary to pass the factors to dgemm_ in reverse order.

  \param[in] transa If transa points to char 'T', the matrix A is transposed
  for the multiplication. Otherwise the untransposed matrix A is used.
  \param[in] transb If transb points to char 'T', the matrix B is transposed
  for the multiplication. Otherwise the untransposed matrix B is used.
  \param[in] m Pointer to an int variable containing the number of rows of A.
  \param[in] n Pointer to an int variable containing the number of columns of B.
  \param[in] k Pointer to an int variable containing the number of columns of A /
  rows of B the matrix is to be performed over.
  \param[in] alpha Pointer to the factor alpha.
  \param[in] A Pointer to the double array containing A.
  \param[in] lda Pointer to an int variable containing the number of colums of A.
  \param[in] B Pointer to the double array containing B.
  \param[in] ldb Pointer to an int variabl containing the number of rows of B.
  \param[in] beta Pointer to the factor beta
  \param[in,out] C Pointer to the double array containing the matrix C.
  \param[in] ldc Pointer to an int variable containing the number rows of C.
 */
extern "C" void dgemm_(char *transa,
                       char *transb,
                       int *m,
                       int *n,
                       int *k,
                       double *alpha,
                       double *A,
                       int *lda,
                       double *B,
                       int *ldb,
                       double *beta,
                       double *C,
                       int *ldc);

extern "C" void zgemm_(char *transa,
                       char *transb,
                       int *m,
                       int *n,
                       int *k,
                       std::complex<double> *alpha,
                       std::complex<double> *A,
                       int *lda,
                       std::complex<double> *B,
                       int *ldb,
                       std::complex<double> *beta,
                       std::complex<double> *C,
                       int *ldc);

//! Matrix-Vector Multiplication
/*!
  Perform one of the matrix-vector operations

      y = alpha * (A x) + beta * y
      y = alpha * (A^T x) + beta * y

  where * denotes multiplication by a scalar and A x matrix-vector
  multiplication.

  \param trans If trans points to the char 'N' or 'n', the standard
  matrix-vector product (A x) is computed. If trans points to 'C','c','T'
  or 't', the matrix-vector product (A^T x) is computed.
  \param m Number of rows of the matrix A.
  \param n Number of columns of the matrix A.
  \param alpha The scalar alpha.
  \param A Pointer to the matrix A.
  \param LDA The leading dimension of A.
  \param x Pointer to the vector x.
  \param incx The stride of the vector x.
  \param beta The scalar beta.
  \param y Pointer to the vector y.
  \param incy The stride of the vector y.
*/
extern "C" void dgemv_(char *trans,
                       int *m,
                       int *n,
                       double *alpha,
                       double *A,
                       int *LDA,
                       double *x,
                       int *incx,
                       double *beta,
                       double *y,
                       int *incy);

extern "C" void zgemv_(char *trans,
                       int *m,
                       int *n,
                       std::complex<double> *alpha,
                       std::complex<double> *A,
                       int *LDA,
                       std::complex<double> *x,
                       int *incx,
                       std::complex<double> *beta,
                       std::complex<double> *y,
                       int *incy);