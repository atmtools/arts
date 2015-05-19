/*!
  \file   blas.h
  \author simon <simon@thinks>
  \date   Thu May  7 22:41:17 2015

  \brief  Interface for BLAS library.

  Contains declarations for BLAS functions.

*/

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
extern "C" void dgemm_( char *transa,
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
			int *ldc );
