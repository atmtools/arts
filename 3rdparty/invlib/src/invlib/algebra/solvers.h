#ifndef ALGEBRA_SOLVERS
#define ALGEBRA_SOLVERS

#include "invlib/algebra.h"
#include <iostream>
#include "invlib/log.h"

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
    VectorType solve(const MatrixType&A,
                     const VectorType& v);

};

// -------------------------------- //
//  Setting Functors for CG Solvers //
// -------------------------------- //

/*! \brief Default Settings for Conjugate Gradient Solvers
 *
 * The default settings for the CG solver initialize the start vector
 * as the 0 vector.
 *
 * As convergence criterion the ratio of the residual norm and the
 * norm of the RHS vector are compared to the tolerance setting
 * of the ConjugateGradient object.
 */
class CGDefaultSettings
{

public:

    CGDefaultSettings(double);

    CGDefaultSettings(const CGDefaultSettings & )             = default;
    CGDefaultSettings(      CGDefaultSettings &&)             = default;
    CGDefaultSettings & operator=(const CGDefaultSettings & ) = default;
    CGDefaultSettings & operator=(      CGDefaultSettings &&) = default;

    /*! Returns the 0 vector as start vector for the CG method. */
    template <typename VectorType>
    VectorType start_vector(const VectorType &) const;

    /*! Check whether the norm of the residual normalized by the RHS norm
     * is smaller than tolerance */
    template <typename VectorType>
    bool converged(const VectorType &, const VectorType &) const;

private:

    double tolerance;

};

/*! \brief Step Limit for CG Solver.
 *
 * Performs a fixed number of steps of the conjugate gradient solver.
 *
 * \tparam maximum_steps The number of steps to be performed.
 */
template<size_t maximum_steps>
class CGStepLimit
{

public:

    CGStepLimit(double);

    CGStepLimit(const CGStepLimit & )             = default;
    CGStepLimit(      CGStepLimit &&)             = default;
    CGStepLimit & operator=(const CGStepLimit & ) = default;
    CGStepLimit & operator=(      CGStepLimit &&) = default;

    /*! Returns the 0 vector as start vector for the CG method. */
    template <typename VectorType>
    VectorType start_vector(const VectorType &);

    /*! Check whether the norm of the residual normalized by the RHS norm
     * is smaller than tolerance */
    template <typename VectorType>
    bool converged(const VectorType &, const VectorType &);

private:

    size_t steps;

};

/*! \brief Step Limit for CG Solver.
 *
 * Performs a fixed number of steps of the conjugate gradient solver.
 *
 * \tparam maximum_steps The number of steps to be performed.
 */
template<typename VectorType, size_t maximum_steps>
class CGContinued
{

public:

    CGContinued(double);

    CGContinued(const CGContinued & )             = default;
    CGContinued(      CGContinued &&)             = default;
    CGContinued & operator=(const CGContinued & ) = default;
    CGContinued & operator=(      CGContinued &&) = default;

    /*! Returns the 0 vector as start vector for the CG method. */
    VectorType & start_vector(const VectorType &);

    /*! Check whether the norm of the residual normalized by the RHS norm
     * is smaller than tolerance */
    bool converged(const VectorType &, const VectorType &);

private:

    size_t     steps;
    VectorType v;

};




// -------------------------  //
//  Conjugate Gradient Solver //
// -------------------------  //

/*! \brief Conjugate gradient solver.
 *
 * The conjugate gradient (CG) solver uses an iterative method to solve the given
 * linear system. Its advantage is that the system does not have to be
 * explicitly computed to be solved.
 *
 * The start vector and stopping criterion are defined using a settings
 * functor. The settings functor  must provide the member functions:
 *
 *  - VectorType start_vector(const VectorType &)
 *  - bool converged(const VectorType &, const VectorType &)
 *
 * that return a start vector and the evaluated convergence criterion,
 * respectively.
 *
 */
template <typename CGSettings = CGDefaultSettings>
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
    ConjugateGradient(double tol, int verbosity = 0);

    /*! Solve linear system using the conjugate gradient method.
     *
     * Takes an arbitrary algebraic expression representing a matrix
     * \f$A\f$ that supports multiplication from the right by a vector
     * and solves the corresponding linear system \f$Ax = v\f$.
     *
     * \tparam MatrixType The algebraic expression type of representing the
     * linear system.
     * \tparam The fundamental vector type. Note: must be the fundamental type.
     * \param A The algebraic expression representing the linear system.
     * \param v The RHS vector \f$v\f$ of the linear system.
     */
    template
    <
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log = StandardLog
    >
    VectorType solve(const MatrixType&A, const VectorType& v);

private:

    int        verbosity;
    double     tolerance;
    CGSettings settings;

};

// ----------------------------------------  //
//  Preconditioned Conjugate Gradient Solver //
// ----------------------------------------  //
/*! \brief Preconditioned Conjugate gradient solver.
 *
 * Solves the linear system
 * \f[
 *  (C^{-1})^T} A C^{-1} x = (C^{-1})^{-T} b
 * \f]
 * whose convergence properties thus depend on the eigenvalues of the matrix
 * ((C^{-1})^TAC^{-1}. For a suitably chosen transformation $C$ this may be
 * used to accelerate the convergence of the conjugate gradient method.
 * The transformation must be provided by the user in the functor that performs
 * the action \f$\mathbf{M}\mathbf{x}\f$ of the matrix
 *
 * \f[
 *  \mathbf{M} = C^{T} C
 * \f]
 * on a given vector \f$\mathbf{x}\f$.
 *
 * \tparam The type of the transformation.
 */
template<typename F, bool Cached = true>
class PreconditionedConjugateGradient;

template<typename F>
class PreconditionedConjugateGradient<F, true>
{

public:

    /*! Create a preconditioned CG solver object with given transformation,
     *  convergence tolerance and verbosity.
     *
     * \param f The functor implementing the action of the preconditioner \f$M = C^{T} C\f$
     * on an arbitrary vector.
     * \param tol The convergence criterion with respect the relative residual
     * \f$\frac{|\mathbf{r}_k|}{|\mathbf{b}|}\f$.
     * \param verbosity If verbosity > 0, iteration progress is printed to standard out.
     */
    PreconditionedConjugateGradient(const F &f, double tol, int verbosity = 0);

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
    template
    <
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log = StandardLog
    >
    VectorType solve(const MatrixType&A, const VectorType& v);

private:

    const F & f;
    int    verbosity;
    double tolerance;
};

template<typename F>
class PreconditionedConjugateGradient<F, false>
{

public:

    /*! Create a non-cached preconditioned CG solver.
     *
     * For each call to the solve(...) member function, this solver
     * will create a new preconditioner. This should be used if the solver
     * is used for an iterative method with varying arguments and the
     * the preconditioner cannot be reused.
     *
     * \param tol The tolerance on the relative residual up to which the
     * iteration is continued.
     * \param verbosity If verbosity > 0, log output is printed to standard out.
     */
    PreconditionedConjugateGradient(double tol, int verbosity = 0);

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
    template
    <
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log = StandardLog
    >
    VectorType solve(const MatrixType&A, const VectorType& v);

private:

    int    verbosity;
    double tolerance;

};

#include "solvers.cpp"

}      // namespace invlib

#endif // ALGEBRA_SOLVERS

