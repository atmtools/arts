/*!
  \file   lapack.h
  \author simon <simon@thinks>
  \date   Thu May  7 22:41:17 2015

  \brief  Interface for the LAPACK library.

  Contains declarations for LAPACK functions.

*/

#include <complex>

namespace lapack {

//! LU decomposition.
/*!
  Performs an LU decomposition of a genereal m-by-n matrix using
  partial pivoting with row interchanges. See LAPACK reference.

  \param[in] m The number of rows of the matrix A.
  \param[in] n The number of column of the matrix A.
  \param[out] A The matrix A.
  \param[in] lda The leading dimension of the matrix A.
  \param[out] ipiv Integer array to hold the pivoting indices.
  \param[out] info Integer indicating if operation was successful: 0 if success,
  otherwise failure.
*/
extern "C" void dgetrf_( int *m,
                         int *n,
                         double *A,
                         int *lda,
                         int *ipiv,
                         int *info );

extern "C" void zgetrf_( int *m,
                         int *n,
                         std::complex<double> *A,
                         int *lda,
                         int *ipiv,
                         int *info );

//! Solve linear system of equations.
/*!
  Solves a linear system  of equations of the form

      trans(A) * x = b

  using the QR decomposition obtained using dgetrf_.

  \param[in] trans Specifies the form of the system of equations:

  trans = 'N' : trans(A) = A
  trans = 'T' : trans(A) = A^T
  trans = 'C' : trans(A) = conj( A^T )

  \param[in] n The size of the system.
  \param[in] nrhs The number of right-hand sides.
  \param[in] A The matrix describing the system.
  \param[in] lda The leading dimension of A.
  \param[in] ipiv The pivot vector as returned by dgetrf_.
  \param[in,out] b The matrix containing the right-hand side vectors.
  \param[in] ldb The leading dimension of b.
  \param[out] info Integer indicating succes of the operation.
*/
extern "C" void dgetrs_( char *trans,
                         int *n,
                         int *nrhs,
                         double *A,
                         int *lda,
                         int *ipiv,
                         double *b,
                         int *ldb,
                         int *info );

//
//! Matrix inversion.
/*!
  Inverts a given matrix using LU decompositions. See LAPACK reference.

  \param[in] n The number of rows and columns of the matrix A.
  \param[in,out] A The matrix to be inverted.
  \param[in] lda The leading dimension of A.
  \param[in] ipiv The vector containing the pivoting indices returned from dgetrf_.
  \param[in] work Work array for improved performance.
  \param[in] lwork Size of the work array.
  \param[out] info Integer indicating if operation was successful: 0 if success,
  otherwise failure.

*/
extern "C" void dgetri_( int* n,
                         double* A,
                         int* lda,
                         int* ipiv,
                         double* work,
                         int* lwork,
                         int* info );

extern "C" void zgetri_( int* n,
                         std::complex<double>* A,
                         int* lda,
                         int* ipiv,
                         std::complex<double>* work,
                         int* lwork,
                         int* info );

//! Optimal parameters for computation.
/*!
  This function returns problem-dependent parameters for the computing
  environment. See LAPACK reference.

  \param[in,out] ispec On input specifies which parameter to be returned and
  on output contains the requested parameter.
  \param[in] name Name of the functions to be called.
  \param[in] opts Character options to the function "Name"
  \param[in] n1 Problem dimension 1 of "Name".
  \param[in] n2 Problem dimension 2 of "Name".
  \param[in] n3 Problem dimension 3 of "Name".
  \param[in] n4 Problem dimension 4 of "Name".
*/
extern "C" void ilaenv_( int *ispec,
                         char *name,
                         char *opts,
                         int *n1,
                         int *n2,
                         int *n3,
                         int *n4 );

//! Solve linear system.
/*!
   Solves the linear system A*X = B using LAPACK's dgesv function. The use of the
   expert function is necessary since we need to solve the system for the inverse
   matrix due to the different memory layouts in arts and BLAS.
   Since most of the parameters are unused they are not fully documented. See
   LAPACK documentation.

  \param[in] fact Specifies whether or not the factored form of the matrix A is
  provided and whether the matrix A should be equilibrated.
  \param[in] trans Specifies the form of the linear system
  = 'N': A * x = B
  = 'N': A^T * x = B
  = 'N': A^H * x = B
  \param[in] n Dimensionality of the system.
  \param[in] nrhs Number of right-hand side columns.
  \param[in] A The double array (column-major order) representing the matrix a A
  that describes the linear system.
  \param[in] lda The leading dimension of A.
  \param[in,out] AF If a is given in factored form AF must contain the factored
  form of A on input. On output contains the factored form of A.
  \param[in] ldaf The leading dimension of AF.
  \param[in,out] ipiv The pivot vector from the LU decomposition.
  \param[in] equed Specifies the form of equilibration that was done.
  \param[in] R The low scale factors for A.
  \param[in] C The column scale factors for A.
  \param[in] B The right hand-side matrix of the system.
  \param[in] ldb The leading dimension of B.
  \param[out] X The solution matrix.
  \param[in] ldx The leading dimension of x.
  \param[out] RCOND The estimate of the reciprocal condition number of the matrix
  A after equilibration (if done).
  \param[out] FERR The estimates error bound for each solution.
  \param[out] BERR The componentwise relative backward error.
  \param[out] WORK
  \param[out] IWORK
  \param[out] info
*/
extern "C" void dgesvx_( char *fact,
                         char * trans,
                         int *n,
                         int *nrhs,
                         double *A,
                         int *lda,
                         double *AF,
                         int *ldaf,
                         int *ipiv,
                         char *equed,
                         double *R,
                         double *C,
                         double *B,
                         int *ldb,
                         double *X,
                         int *ldx,
                         double *RCOND,
                         double *FERR,
                         double *BERR,
                         double *WORK,
                         int *IWORK,
                         int *info );

/* Computes eigenvalues and eigenvectors for the Complex n-by-n Matrix A
 * Use-case example diag(VR*A*inv(VR)) = W.

    \param[in] jobvl calculate left eigenvectors if 'V', otherwise 'N'
    \param[in] jobvl calculate right eigen vectors if 'V', otherwise 'N'
    \param[in] n Dimensionality of the system.
    \param[in,out] A matrix to find eigenvalues and eigenvectors.  Changed in execution.
    \param[in] lda The leading dimension of A.
    \param[out] W/WR/WI The eigenvalues in a vector (double matrix can have complex eigenvalues)
    \param[out] VL The left eigenvectors of A.
    \param[in] ldvl The leading dimension of VL.
    \param[out] VR The right eigenvectors of A.
    \param[in] ldvr The leading dimension of VR.
    \param[out] WORK
    \param[out] IWORK
    \param[out] RWORK
    \param[out] info
 */

extern "C" void dgeev_(  char *jobvl,
                         char *jobvr,
                         int  *n,
                         double *A,
                         int *lda,
                         double *WR,
                         double *WI,
                         double *VL,
                         int *ldvl,
                         double *VR,
                         int *ldvr,
                         double *work,
                         int *lwork,
                         double *rwork,
                         int *info );

extern "C" void zgeev_(  char *jobvl,
                         char *jobvr,
                         int  *n,
                         std::complex<double> *A,
                         int *lda,
                         std::complex<double> *W,
                         std::complex<double> *VL,
                         int *ldvl,
                         std::complex<double> *VR,
                         int *ldvr,
                         std::complex<double> *work,
                         int *lwork,
                         double *rwork,
                         int *info );
}
