#ifndef ALGEBRA_SOLVERS
#define ALGEBRA_SOLVERS

#include "invlib/algebra.h"
#include <iostream>

namespace invlib
{

/** file solvers.h
 * \brief Solver for linear systems.
 *
 * This file contains class that provide solvers for linear systems of the form
 *
 * \f[
 *  A x = b
 * \f]
 *
 * where \f$ A\f$ an full rank, square matrix \f$ A \in \mathbb{R}^{n \times n}\f$
 * and \f$ b \in \mathbb{R}^n \f$ a vector.
 *
 */

// ------------------ //
//  Standard Solver   //
// ------------------ //

/** \brief Standard solver forwarding to underlying member functions.
 *
 * The Standard solver simply forwards the call to solve() to the underlying
 * member function MatrixType::solve(const VectorType&).
 */
class Standard
{
public:

    /*! Solve the linear system using the solve method of the fundamental Matrix
    *  type.
    */
    template < typename VectorType, typename MatrixType>
    VectorType solve(const MatrixType&A, const VectorType& v);

};

// -------------------------  //
//  Conjugate Gradient Solver //
// -------------------------  //

/*! \brief Conjugate gradient solver.
 *
 * The conjugate gradient (CG) solver uses an iterative method to solve the given
 * linear system. Its advantage is that the system does not have to be
 * explicitly computed to be solved.The convergence criterion used is the
 * Euclidean norm of the residual.
 *
 */
class ConjugateGradient
{

public:

    /*! Create CG solver object with given convergence tolerance and
     * given coordinate transform, to be applied to the system.
     *
     * \param tol The convergence tolerance
     * \param trans The coordinate transformation. Defaults to the identity
     * transformation.
     */
    ConjugateGradient(double tol);

    /*! Solve linear system using the conjugate gradient method.
     *
     * Takes an arbitrary algebraic expression representing a matrix
     * \f$A\f$ that supports multiplication from the right by a vector
     * and solves the corresponding linear system \f$Ax = v\f$. The
     * iteration is stopped when the norm of the residual
     * \f$r = Ax - v \f$ falls below the given convergence tolerance.
     *
     * \tparam MatrixType The algebraic expression type of representing the
     * linear system.
     * \tparam The fundamental vector type. Note: must be the fundamental type.
     * \param A The algebraic expression representing the linear system.
     * \param v The RHS vector \f$v\f$ of the linear system.
     */
    template <typename VectorType, typename MatrixType>
    VectorType solve(const MatrixType&A, const VectorType& v);

private:

    double tolerance;

};

template
<
typename Solver,
typename Transformation
>
class PreconditionedSolver
{

public:

    PreconditionedSolver(Solver s, Transformation t);

    template <typename VectorType, typename MatrixType>
    VectorType solve(const MatrixType &A, const VectorType &v);

private:
    Solver solver;
    Transformation transformation;
};

#include "solvers.cpp"

}      // namespace invlib

#endif // ALGEBRA_SOLVERS

