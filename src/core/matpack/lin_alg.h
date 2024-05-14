/*!
     \file   lin_alg.h
     \author Claudia Emde <claudia.emde@dlr.de>
     \date   Thu May  2 14:34:05 2002

     \brief  Linear algebra functions.

   */

#ifndef linalg_h
#define linalg_h

#include "array.h"
#include "matpack_data.h"

// LU decomposition
void ludcmp(Matrix& LU, ArrayOfIndex& indx, ConstMatrixView A);

// LU backsubstitution
void lubacksub(VectorView x,
               ConstMatrixView LU,
               ConstVectorView b,
               const ArrayOfIndex& indx);

// Solve linear system
void solve(VectorView x, ConstMatrixView A, ConstVectorView b);

struct solve_workdata {
  std::size_t N{};
  std::vector<int> ipiv{};

  constexpr solve_workdata() = default;
  constexpr solve_workdata(const solve_workdata&) = default;
  constexpr solve_workdata(solve_workdata&&) = default;
  constexpr solve_workdata& operator=(const solve_workdata&) = default;
  constexpr solve_workdata& operator=(solve_workdata&&) = default;

  constexpr solve_workdata(std::size_t N_) : N(N_), ipiv(N) {}
  constexpr void resize(std::size_t N_) {
    N = N_;
    ipiv.resize(N);
  }
};

/*! Solves A X = B inplace using dgesv.
  * 
  * Returns the Lapack ipiv array.
  *
  * @param[in,out] X   As equation, on input it is B on output is is X
  * @param[in]     A   As equation, it is destroyed on output (LU decomposition)
  * @throws If the system cannot be solved according to Lapack info
  */
void solve_inplace(ExhaustiveVectorView X,
                   ExhaustiveMatrixView A,
                   solve_workdata& wo);

//! As above but allocates WO
void solve_inplace(ExhaustiveVectorView X, ExhaustiveMatrixView A);

struct inv_workdata {
  std::size_t N{};
  std::vector<int> ipiv{};
  std::vector<Numeric> work{};

  constexpr inv_workdata() = default;
  constexpr inv_workdata(const inv_workdata&) = default;
  constexpr inv_workdata(inv_workdata&&) = default;
  constexpr inv_workdata& operator=(const inv_workdata&) = default;
  constexpr inv_workdata& operator=(inv_workdata&&) = default;

  constexpr inv_workdata(std::size_t N_) : N(N_), ipiv(N), work(N) {}
  constexpr void resize(std::size_t N_) {
    N = N_;
    ipiv.resize(N);
    work.resize(N);
  }
};

// Matrix inverse
void inv(MatrixView Ainv, ConstMatrixView A);

// Matrix inverse in place with destructive consequences
void inv_inplace(ExhaustiveMatrixView A);

// Matrix inverse in place with destructive consequences
void inv_inplace(ExhaustiveMatrixView A, inv_workdata& wo);

// Matrix inverse
void inv(ComplexMatrixView Ainv, const ConstComplexMatrixView A);

struct diagonalize_workdata {
  std::size_t N{};
  std::vector<Numeric> w{};

  constexpr diagonalize_workdata() = default;
  constexpr diagonalize_workdata(const diagonalize_workdata&) = default;
  constexpr diagonalize_workdata(diagonalize_workdata&&) = default;
  constexpr diagonalize_workdata& operator=(const diagonalize_workdata&) =
      default;
  constexpr diagonalize_workdata& operator=(diagonalize_workdata&&) = default;

  constexpr diagonalize_workdata(std::size_t N_) : N(N_), w(4 * N + N * N) {}
  constexpr Numeric* work() { return w.data(); }
  constexpr Numeric* rwork() { return w.data() + 2 * N; }
};

// Matrix diagonalization with lapack
void diagonalize(MatrixView P, VectorView WR, VectorView WI, ConstMatrixView A);

// Matrix diagonalization with lapack
void diagonalize(MatrixView P,
                 VectorView WR,
                 VectorView WI,
                 ConstMatrixView A,
                 diagonalize_workdata& wo);

// Same as diagonalize but inplace manilpulation of input with destructive consqeuences
void diagonalize_inplace(ExhaustiveMatrixView P,
                         ExhaustiveVectorView WR,
                         ExhaustiveVectorView WI,
                         ExhaustiveMatrixView A);

// Same as diagonalize but inplace manilpulation of input with destructive consqeuences
void diagonalize_inplace(ExhaustiveMatrixView P,
                         ExhaustiveVectorView WR,
                         ExhaustiveVectorView WI,
                         ExhaustiveMatrixView A,
                         diagonalize_workdata& wo);

// Matrix diagonalization with lapack
void diagonalize(ComplexMatrixView P,
                 ComplexVectorView W,
                 const ConstComplexMatrixView A);

// Exponential of a Matrix
void matrix_exp(MatrixView F, ConstMatrixView A, const Index& q = 10);

// 2-norm of a vector
Numeric norm2(ConstVectorView v);

// Maximum absolute row sum norm
Numeric norm_inf(ConstMatrixView A);

// Identity Matrix
void id_mat(MatrixView I);

Numeric det(ConstMatrixView A);

void linreg(Vector& p, ConstVectorView x, ConstVectorView y);

/** Least squares fitting by solving x for known A and y
 * 
 * (A^T A)x = A^T y
 * 
 * Returns the squared residual, i.e., <(A^T A)x-A^T y|(A^T A)x-A^T y>.
 * 
 * @param[in]  x   As equation
 * @param[in]  A   As equation
 * @param[in]  y   As equation
 * @param[in]  residual (optional) Returns the residual if true
 * @return Squared residual or 0
 */
Numeric lsf(VectorView x,
            ConstMatrixView A,
            ConstVectorView y,
            bool residual = true) noexcept;

#endif  // linalg_h
