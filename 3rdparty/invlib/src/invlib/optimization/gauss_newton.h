/** \file optimization/gauss_newton.h
 *
 * \brief Contains GaussNewton class template implementing the Gauss-Newton
 * optimization scheme.
 *
 */

#ifndef OPTIMIZATION_GAUSS_NEWTON_H
#define OPTIMIZATION_GAUSS_NEWTON_H

#include <stdio.h>
#include "invlib/algebra/solvers.h"

namespace invlib
{

/**
 * \brief Gauss-Newton method.
 *
 * Class template for a generic Gauss-Newton optimizer. Provides the function
 * step, which computes one step \f$ d\vec{x} \f$ of the Gauss-Newton method:
 *
 * \f[
 *    d\vec{x} &= -\mathbf{H}^{-1} \vec{g}
 * \f]
 *
 * where \f$ \mathbf{H} \f$ is the (approximate) Hessian of the cost function
 * and \f$ \vec{g} \f$ its gradient. The next iteration vector can then be
 * computed using \f$ \vec{x}_{i+1} = \vec{x}_{i} + d\vec{x} \f$.
 *
 * The method used for the solution of the subproblem
 *
 * \f[
 *    \mathbf{H}} d\vec{x} = -\vec{g}
 * \f]
 *
 * can be defined using the Solver type. The default is to use the
 * the solve() member function of the given Hessian.
 *
 * \tparam Real The floating point type used to represent scalars.
 * \tparam Solver The Solver type to be used for the subproblem.
 */
template
<
typename RealType,
typename Solver = Standard
>
class GaussNewton
{

public:

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Construct GaussNewton optimizer with default tolerance of
     * 1e-5 and iteration maximum 1000. If no solver argument is given,
     * the Solver type must provide a default constructor that takes no
     * arguments to create the solver object used for the minimization
     * step.
     *
     * \param solver_ The solver to be used to solve the linear subproblem. If not
     * provided the default constructed solver Solver() is used.
     */
    GaussNewton(Solver solver_ = Solver());

    /*! Construct GaussNewton optimizer with given tolerance, maximum
    iterations and * solver used to solve the linear subproblem
    subproblem. If no solver arguments is provided, * the default
    constructed solver Solver() is used as default.
     *
     * \param tolerance The tolerance used to determine the
     * convergence of the optimization method.  \param
     * maximum_iterations_ The maximum of iterations to perform before
     * the minimization process should be aborted.
     * \param solver A solver object used to solve the subproblem of the
     * minimization.
     */
    GaussNewton( RealType tolerance,
                 unsigned int maximum_iterations_,
                 Solver solver_ = Solver() );

    // -------------------------- //
    //    Getters and Setters     //
    // -------------------------- //

    unsigned int get_maximum_iterations() const;
    void set_maximum_iterations(unsigned int n);

    RealType get_tolerance() const;
    void set_tolerance(RealType tolerance_);

    // --------------------------- //
    //  Perform Minimization Step  //
    // --------------------------- //

    /*! Compute minimization step \f$d\vec{x}\f$ for the Gauss-Newton
     *  method.  This method should be called in an iterative
     *  optimization loop. The next iteration vector \f$\vec{x}_{i+1}\f$ can
     *  then be computed using  \f$\vec{x}_{i+1} = \vec{x}_{i} + d\vec{x}\f$.
     *  Only uses the provided gradient vector \f$g\f$ and the (approximate)
     *  Hessian \f$B\f$, the other arguments are only provided for compatibilty
     *  with other optimization methods.
     *
     * \tparam VectorType A fundamental vector. Must not be an algebraic
     * expression.
     * \tparam MatrixType The type of the (approximate) Hessian matrix. May
     * be any algebraic expression that represents a matrix.
     * \param g The gradient vector at the current iteration step.
     * \param B An expression for the (approximate) Hessian.
     */
    template
    <
    typename VectorType,
    typename MatrixType,
    typename CostFunction
    >
    VectorType step(const VectorType &,
                    const VectorType &g,
                    const MatrixType &B,
                    const CostFunction &);

private:

    RealType tolerance;
    unsigned int maximum_iterations;
    Solver solver;

};

#include "gauss_newton.cpp"

}      // namespace invlib

#endif //OPTIMIZATION_GAUSS_NEWTON
