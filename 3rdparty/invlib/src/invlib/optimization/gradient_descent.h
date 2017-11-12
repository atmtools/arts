/** \file optimization/gradient_descent.h
 *
 * \brief Contains GradientDescent class template implementing the gradient
 * descent optimization.
 *
 */

#ifndef OPTIMIZATION_GRADIENT_DESCENT_H
#define OPTIMIZATION_GRADIENT_DESCENT_H

#include <stdio.h>

namespace invlib
{

/*!
 * Implements gradient descent optimization, i.e. the iteration step \f$d\vec{x}\f$
 * computed by the step(...) member function yields
 *
 * \f[
 *    d\vec{x} &= -\alpha \vec{g}
 * \f]
 *
 * where \alpha is a suitable step length so that the objective function is
 * is reduced.
 *
 * \tparam Real The floating point type used to represent scalars.
 */
template
<
typename RealType
>
class GradientDescent
{

public:

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Construct GradientDescent optimizer with default tolerance of
     * 1e-5 and iteration maximum 1000.
     */
    GradientDescent();

    /*! Construct GradientDescent optimizer with given initial step length.
     *
     * \param initial_step The step length to start the line search in the first
     * iteration.
     * \param maximum_iterations The maximum number of operations to be performed
     * by the optimizer.
     */
    GradientDescent(RealType tolerance,
                    unsigned int maximum_iterations,
                    RealType initial_step = 0.1,
                    RealType scale = 10);

    // -------------------------- //
    //    Getters and Setters     //
    // -------------------------- //

    unsigned int get_maximum_iterations() const;
    void set_maximum_iterations(unsigned int n);

    RealType get_tolerance() const;
    void set_tolerance(RealType tolerance_);

    RealType get_initial_step() const;
    void set_initial_step(RealType initial_step);

    RealType get_scale() const;
    void set_scale(RealType );

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
                    CostFunction & J);

    bool stop_iteration() {return false;}

private:

    RealType tolerance, step_length, scale, minimum_step_length, current_cost;
    unsigned int maximum_iterations, step_count;

};

#include "gradient_descent.cpp"

}      // namespace invlib

#endif //OPTIMIZATION_GRADIENT_DESCENT
