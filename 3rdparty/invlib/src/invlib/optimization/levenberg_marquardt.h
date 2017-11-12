/** \file optimization/levenberg_marquardt.h
 *
 * \brief Contains the LevenbergMarquardt class template implementing the
 * Gauss-Newton optimization scheme.
 *
 */

#ifndef OPTIMIZATION_LEVENBERG_MARQUARDT_H
#define OPTIMIZATION_LEVENBERG_MARQUARDT_H

#include "invlib/algebra/solvers.h"
#include <iostream>
#include <limits>

namespace invlib
{

/**
 * \brief Levenberg-Marquardt method.
 *
 * Class template for a generic Levenberg-Marquardt optimizer. Provides the
 * function step, which computes one step \f$ d\vec{x} \f$ of the
 * Levenberg-Marquardt method:
 *
 * \f[
 *    d\vec{x} &= -(\mathbf{H} + \lambda \mathbf{D})^{-1} \vec{g}
 * \f]
 *
 * where \f$ \mathbf{H} \f$ is the (approximate) Hessian of the cost function
 * and \f$ \vec{g} \f$ its gradient and \f$ \mathbf{D} \f$ is a user provided
 * positive definite matrix.
 *
 * The value of \f$ \lambda \f$ is adapted depending on how well the cost
 * function can be approximated by a quadratic function.
 *
 * The next iteration vector can be computed from \f$ d\vec{x} \f$ using
 * \f$ \vec{x}_{i+1} = \vec{x}_{i} + d\vec{x} \f$.
 *
 * The method used for the solution of the subproblem
 *
 * \f[
 *    \mathbf{H}} d\vec{x} = -\vec{g}
 * \f]
 *
 * can be defined using the @Solver type. The default is to use the
 * the solve() member function of the given Hessian.
 *
 * \tparam RealType The floating point type used to represent scalars.
 * \tparam Solver The Solver type to be used for the subproblem.
 */
template
<
typename RealType,
typename DampingMatrix,
typename Solver = Standard
>
class LevenbergMarquardt
{

public:

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Construct Levenberg-Marquardt optimizer with given damping
     * matrix and solver. If no solver is provided to default constructed
     * Solver object created from Solver() is used.
     */
    LevenbergMarquardt(const DampingMatrix &D_, Solver solver = Solver());

    // ------------------------- //
    //    Getters and Setters    //
    // ------------------------- //

    unsigned int get_maximum_iterations() const;
    void set_maximum_iterations(unsigned int);

    RealType get_tolerance() const;
    void set_tolerance(RealType);

    RealType get_lambda() const;
    void set_lambda(RealType);

    RealType get_lambda_maximum() const;
    void set_lambda_maximum(RealType);

    RealType get_lambda_decrease() const;
    void set_lambda_decrease(RealType);

    RealType get_lambda_increase() const;
    void set_lambda_increase(RealType);

    RealType get_lambda_threshold() const;
    void set_lambda_threshold(RealType);

    RealType get_lambda_constraint() const;
    void set_lambda_constraint(RealType);

    // --------------------------- //
    //  Perform Minimization Step  //
    // --------------------------- //

    /*! Perform Levenberg-Marquardt step. Computes a tentaive step \f$d\vec{x}\f$
     * by solving the system
     *
     * \f[
     *    d\vec{x} &= -(\mathbf{H} + \lambda \mathbf{D})^{-1} \vec{g}
     * \f]
     *
     * If the step reduces the cost function at
     * \f$\vec{x}_{i+1} = \vec{x}_i + d\vec{x}\f$, the step is accepted and the
     * current lambda values reduced by a factor of lambda_decrease. If the value
     * of the cost function is not decreased the values for lambda is increased
     * by a factor of lambda_increase and \f$d\vec{x}\f$ is recomputed. This is
     * repeated until a suitable \f$d\vec{x}\f$ is found or lambda reaches the
     * the value of lambda_maximum. If lambda falls below lambda_threshold, lambda
     * is set to zero and the Levenberg-Marquardt step effectively becomes a
     * Gauss-Newton step.
     */
    template
    <
    typename VectorType,
    typename MatrixType,
    typename CostFunction
    >
    VectorType step(const VectorType &x,
                    const VectorType &g,
                    const MatrixType &B,
                    CostFunction     &J);

    bool stop_iteration() const {return stop;}

private:

    RealType current_cost, tolerance, lambda, lambda_maximum, lambda_increase,
    lambda_decrease, lambda_threshold, lambda_constraint;
    unsigned int maximum_iterations, step_count;
    bool stop;

    // Positive definite matrix defining the trust region sphere r < ||Mx||.
    const DampingMatrix &D;

    Solver s;

};

#include "levenberg_marquardt.cpp"

}      // namespace invlib

#endif // OPTIMIZATION_LEVENBERG_MARQUARDT_H

