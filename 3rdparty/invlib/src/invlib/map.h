#ifndef MAP_H
#define MAP_H

#include <iostream>
#include <chrono>

using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include "invlib/algebra.h"
#include "invlib/algebra/solvers.h"
#include "invlib/log.h"
#include "invlib/traits.h"

/** file map.h
 * \brief Maximum A Posteriori Estimators
 *
 * This file provides class templates for the inversion of a given forward model
 * using maximum a posteriori estimators. Under the assumption of a Gaussian
 * prior and measurement error, that means that a solution of the inverse problem
 * is obtained by minimizing the cost function
 *
 * \f[
 *   J(\mathbf{x}) &= \frac{1}{2}
 *   \left((\mathbf{F}(\mathbf{x}) - \mathbf{y}) \mathbf{S}_e
 *         (\mathbf{F}(\mathbf{x}) - \mathbf{y})
 *        +(\mathbf{x} - \mathbf{x}_a)\mathbf{S}_a
 *         (\mathbf{x} - \mathbf{x}_a) \right)
 * \f]
 *
 * To this end, this file provides the class template MAP, and three partial
 * specializations, one for each of the three different possible formulations of
 * the problem, here called the standard, n-form and m-form as given in formulas
 * (5.8), (5.9), (5.10) in \cite rodgers, respectively.
 *
 * Methods common to all three methods are aggregated in the base class MAPBase.
 *
 */

namespace invlib
{

/** \brief Formulation enum specifiying which formulation of the MAP estimator
 * to use.
 *
 * For details on the form see template specializations.
 */
enum class Formulation {STANDARD = 0, NFORM = 1, MFORM = 2};

/*!
 * Exception class representing falure of the forward model evaluation.
 */
struct ForwardModelEvaluationException
{
    std::string message()
    {
        return "Could not evaluate forward model";
    }
};

/**
 * \brief MAP base class
 *
 * Implements methods common to all MAP estimators independent of
 * formulation.  Holds references to the forward model, the a priori
 * state vector as well as the a priori and measurement space
 * covariance matrices. Provides a member to hold a pointer to the
 * measurement vector. This is necessary for the class to be able to
 * provide the cost_function(const VectorType &x) function taking only
 * the current state space vector as an argument, which is used by the
 * optimizer to minimize the cost function for a given measurement
 * vector.
 *
 * To allow for maximum flexibility the type of the ForwardModel used is a
 * class template parameter. To be used with the MAP class, the functions
 * ForwardModel type must provide the following member functions:
 *
 * - VectorType evaluate(const VectorType& x): Evaluate the forward model at the given
 * state space vector.
 * - MatrixType Jacobian(const VectorType& x, VectorType&y): Compute the Jacobian
 * of the forward model and set the forward model prediciton y.
 *
 * \tparam ForwardModel The forward model type to be used.
 * \tparam RealType   The floating point type used for scalars.
 * \tparam VectorType The vector type used for vectors.
 * \tparam MatrixType The matrix type used.
 * \tparam SaType The type of the a priori matrix used for the computations.
 * \tparam SeType The type of the measurement space covariance matrix to be
 * used.
 */
template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType = typename MatrixType::VectorType
>
class MAPBase
{

public:

    /*! The basic scalar type. */
    using RealType   = typename MatrixType::RealType;

private:

    // Helper functions to determine type of gradient and Jacobian
    // as returned by the forward model.
    auto evaluate_helper(ForwardModel & f, const VectorType & x)
        -> decltype(f.evaluate(x));
    auto Jacobian_helper(ForwardModel & f, const VectorType & x, VectorType & y)
        -> decltype(f.Jacobian(x,y));

public:

    /*! The type of the gradient vector as returned by the forward model. */
    using FMVectorType =
        return_type<decltype(&MAPBase::evaluate_helper)(MAPBase,
                                                        ForwardModel&,
                                                        const VectorType &)>;

    /*! The type of the Jacobian matrix as returned by the forward model. */
    using FMJacobianType =
        return_type<decltype(&MAPBase::Jacobian_helper)(MAPBase,
                                                        ForwardModel &,
                                                        const VectorType &,
                                                        VectorType &)>;

    /*! Jacobian type that can be assigned to. */
    using JacobianType = CopyWrapper<FMJacobianType>;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Construct MAP estimator for the given forward model F, a priori vector
     *  xa and a priori covariance matrix Sa and measurement space covariance
     * matrix Se.
     *
     * \param F Reference to the forward model object that provides
     * evaluate(...) and Jacobian(...) member function.
     * \param xa The a priori state vector.
     * \param Sa the a priori covariance matrix.
     * \param The measurement covariance matrix.
     */
    MAPBase(ForwardModel &F_,
            const VectorType   &xa_,
            const SaType &Sa_,
            const SeType &Se_ );

    // ----------------------- //
    //  MAP Utility Functions  //
    // ----------------------- //

    /*! Evaluate MAP cost function.
     *
     * \param x A given state vector of dimension m.
     * \param y The measurement vector y.
     * \param yi The forward model prediction of the measurement vector
     * \$f\mathbf{F}(\vec{x}) = \vec{y}_i \$f.
     */
    RealType cost_function(const VectorType &x,
                           const VectorType &y,
                           const VectorType &yi);
    /*! Evaluate MAP cost function.
     *
     * Evaluates forward model at given vector x to obtain the
     * prediction \f$\vec{y}_i\f$. Evaluates the cost function using x
     * and \f$\vec{y}_i\f$ and the measurement vector currently pointed to
     * by the y_ptr private memeber, which is set in the compute(...)
     * member function.
     *
     * \param x A given state vector of dimension m.
     */
    RealType cost_function(const VectorType &x);

    /*! Compute the the contribution of the estimated state vector
     * to the total cost function, which is given by
     *
     * \f[
     *     J_\vec{x}(\vec{x}) &= (\vec{x} - \vec{x}_a)^T
     *                           S_a^{-1}(\vec{x} - \vec{x}_a)
     * \f]
     *
     * \param x The vector \f$x\f$ at which to evaluate the cost function.
     * \return The state space related part of the cost function as given
     * by the equation above.
     */
    RealType cost_x(const VectorType &x);

    /*! Compute the the contribution of the deviation in measurement space
     * to the total cost function, which is given by
     *
     * \f[
     *     J_\vec{y}(\vec{y},\vec{x}) &= (\mathbf{F}(\vec{x}) - \vec{y})^T
     *                           S_\epsilon^{-1}(\mathbf{F}(\vec{x}) - \vec{y})
     * \f]
     *
     * \param y The measurement vector \f$\vec{y}\f$
     * \param yi The predicted measurement vector \f$\vec{y}_i =
     * \mathbf{F}(\vec{x})\f$ as given by the forward model.
     * \return The measurement space related part of the cost function as given by the
     * formula above.
     */
    RealType cost_y(const VectorType &y,
                    const VectorType &yi);

    /*! Exception safe wrapper for the evaluate function of the forward
     * model.
     */
    FMVectorType evaluate(const VectorType &x);

    /*! Exception safe wrapper for the Jaobian computation function of the
     * forward model.
     */
    JacobianType Jacobian(const VectorType &x, VectorType &y);

    /*! Compute the gain matrix at the given state vector x.
     *
     * Computes the gain matrix
     * \f[
     *    G = (K^T {S}_\epsilon K + S_a)^{-1}S_\epsilon K \vec{x}
     * \f]
     *
     * \param x The state vector at which to evaluate the gain matrix.
     */
    MatrixType gain_matrix(const VectorType &x);

protected:

    unsigned int m, n;

    ForwardModel &F;
    const VectorType   &xa;
    const VectorType   *y_ptr;
    const SaType &Sa;
    const SeType &Se;

    duration<double> evaluate_time;
    duration<double> Jacobian_time;

};

// -------------- //
//   MAP Class    //
// -------------- //

template
<
typename ForwardModel,
typename MatrixType,
typename SaType     = MatrixType,
typename SeType     = MatrixType,
typename VectorType = typename MatrixType::VectorType,
Formulation Form    = Formulation::STANDARD
>
class MAP;

// ------------- //
//   Standard    //
// ------------- //

/*! MAP Estimator using the standard form as given by formula (5.8) in
 * \cite Rodgers.
 *
 * The standard form is the most flexible form in terms that is can be used
 * with any optimization method, whereas the two other methods only support
 * the Gauss-Newton minimizer. The formulation uses the inverses of the
 * covariance matrices so it is advisable to give the precision matrices
 * directly in order to avoid repeated computation of their inverses. Each step
 * requires the solution of an n-times-n system of linear equations.
 */
template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
class MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::STANDARD>
    : public MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
{

public:

    /*! The basic scalar type. */
    using RealType   = typename MatrixType::RealType;
    /*! The base class. */
    using Base = MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>;
    /*! The Jacobian Type */
    using typename Base::JacobianType;
    /*! The vector type as returned by the forward model. */
    using typename Base::FMVectorType;

    /*! Make Base memeber directly available. */
    using Base::m; using Base::n;
    using Base::y_ptr; using Base::xa;
    using Base::F; using Base::Sa; using Base::Se;
    using Base::cost_function;
    using Base::evaluate; using Base::Jacobian;
    using Base::evaluate_time; using Base::Jacobian_time;

    MAP( ForwardModel &F_,
         const VectorType &xa_,
         const SaType &Sa_,
         const SeType &Se_ );

    /*! Compute the maximum a posteriori estimator for the inverse problem
     * represented by this MAP object and the given measurement vector
     * \f$\vec{y}\f$ using the given minimizer.
     *
     * Since this specialization uses the standard form, in each
     * iteration, the next state vector \f$\vec{x}_{i+1}\f$ is
     * computed using the formula
     *
     * \f[
     *      \vec{x}_{i+1} = x_i  - (K^T S_\epsilon^{-1} K + S_a^{-1})^{-1}
     *                       K^TS_\epsilon^{-1} (F(\vec{x_i}) - \vec{y})
     *                       +S_a^{-1}(\vec{x} - \vec{x}_i)
     * \f]
     *
     * which requires the solution of a n-times-n linear system of equations.
     *
     * \param x The maximum a posteriori state vector.
     * \param y The measured measurement vector.
     * \param M A minimizer object of representing the minimization method
     * gradient before solving the subproblem.
     * that should be used to minimize the likelihood.
     */
    template<typename Minimizer, template <LogType> class Log = StandardLog>
    int compute(VectorType       &x,
                const VectorType &y,
                Minimizer M,
                int verbosity = 0);

    RealType cost, cost_x, cost_y;
    unsigned int iterations;

};

// --------------- //
//     N-form      //
// --------------- //

/*! MAP Estimator using the n-form form as given by formula (5.9) in
 * \cite Rodgers.
 *
 * Note that the n-form can only be used with the Gauss-Newton optimizer.
 * As with the standard form the formulation uses the inverses of the
 * covariance matrices, so it is advisable to give the precision matrices
 * when defining the inverse problem.
 *
 * Each step requires the solution of an n-times-n linear system of
 * equations.
 *
 */
template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
class MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::NFORM>
    : public MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
{

public:

    /*! The basic scalar type. */
    using RealType   = typename MatrixType::RealType;
    /*! The base class. */
    using Base = MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>;
    /*! The Jacobian Type */
    using typename Base::JacobianType;
    /*! The vector type as returned by the forward model. */
    using typename Base::FMVectorType;

    /*! Make Base memeber directly available. */
    using Base::m; using Base::n;
    using Base::y_ptr; using Base::xa;
    using Base::F; using Base::Sa; using Base::Se;
    using Base::cost_function;
    using Base::evaluate; using Base::Jacobian;
    using Base::evaluate_time; using Base::Jacobian_time;

    MAP( ForwardModel &F_,
         const VectorType   &xa_,
         const SaType &Sa_,
         const SeType &Se_ );

    /*! Compute the maximum a posteriori estimator for the inverse problem
     * represented by this MAP object and the given measurement vector
     * \f$\vec{y}\f$ using the given minimizer.
     *
     * This specialization uses the n-form, there for in each
     * iteration, the next state vector \f$\vec{x}_{i+1}\f$ is
     * computed using the formula
     *
     * \f[
     *      \vec{x}_{i+1} =  x_i - (K^T S_\epsilon^{-1} K + S_a^{-1})^{-1}
     *                       K^TS_\epsilon^{-1} (F(\vec{x_i}) - \vec{y}
     *                                         + K(\vec{x} - \vec{x}_i))
     * \f]
     *
     * which requires the solution of an n-times-n system of linear equations.
     *
     * \param x The maximum a posteriori state vector.
     * \param y The measured measurement vector.
     * \param M A minimizer object of representing the minimization method
     * gradient before solving the subproblem.
     * that should be used to minimize the likelihood.
     */
    template<typename Minimizer, template <LogType> class Log = StandardLog>
    int compute(VectorType       &x,
                const VectorType &y,
                Minimizer M,
                int verbosity = 0);

    RealType cost, cost_x, cost_y;
    unsigned int iterations;

};

// --------------- //
//     M-form      //
// --------------- //

/*! MAP Estimator using the m-form form as given by formula (5.10) in
 * \cite Rodgers.
 *
 * Note that the m-form can only be used with the Gauss-Newton optimizer.
 * The formulation uses the covariance matrices directly, so there is no
 * need to provide precision matrices when defining the inverse problem.
 *
 * Each step requires the solution of an m-times-m linear system of
 * equations.
 *
 */
template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
class MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::MFORM>
    : public MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
{

public:

    using RealType   = typename MatrixType::RealType;
    /*! The base class. */
    using Base = MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>;
    /*! The assignable Jacobian type. */
    using typename Base::JacobianType;
    /*! The vector type as returned by the forward model. */
    using typename Base::FMVectorType;

    /*! Make Base memeber directly available. */
    using Base::m; using Base::n;
    using Base::y_ptr; using Base::xa;
    using Base::F; using Base::Sa; using Base::Se;
    using Base::cost_function;
    using Base::evaluate; using Base::Jacobian;
    using Base::evaluate_time; using Base::Jacobian_time;

    MAP( ForwardModel &F_,
         const VectorType   &xa_,
         const SaType &Sa_,
         const SeType &Se_ );

    MatrixType gain_matrix(const VectorType &x);

    RealType cost_function(const VectorType &x,
                           const VectorType &y,
                           const VectorType &yi);
    RealType cost_function(const VectorType &x);

    /*! Compute the maximum a posteriori estimator for the inverse problem
     * represented by this MAP object and the given measurement vector
     * \f$\vec{y}\f$ using the given minimizer.
     *
     * This specialization uses the m-form, there for in each
     * iteration, the next state vector \f$\vec{x}_{i+1}\f$ is
     * computed using the formula
     *
     * \f[
     *      \vec{x}_{i+1} =  x_i - S_aK^T(K S_a K^T + S_e)^{-1}
     *                       (F(\vec{x_i}) - \vec{y}
     *                                         + K(\vec{x} - \vec{x}_i))
     * \f]
     *
     * which requires the solution of an m-times-m system of linear equations.
     *
     * \param x The maximum a posteriori state vector.
     * \param y The measured measurement vector.
     * \param M A minimizer object of representing the minimization method
     * gradient before solving the subproblem.
     * that should be used to minimize the likelihood.
     */
    template<typename Minimizer, template <LogType> class Log = StandardLog>
    int compute(VectorType       &x,
                const VectorType &y,
                Minimizer M,
                int verbosity = 0);

    RealType cost, cost_x, cost_y;
    unsigned int iterations;

};

#include "map.cpp"

}      // namespace invlib

#endif // MAP_H
