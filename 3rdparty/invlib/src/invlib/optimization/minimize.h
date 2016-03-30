/** \file optimization/minimize.h
 *
 * \brief Provides the function template minimize() which implements a genereic
 * minimization algorithm using a provided optimizer object.
 *
 */

#ifndef OPTIMIZATION_MINIMIZE_H
#define OPTIMIZATION_MINIMIZE_H

namespace invlib
{

/**
 * \brief Generic minimization
 *
 * Minimizes the given cost function @J using the given minimizer @M. Generic
 * implementation of an iterative minimization procedure, where in each
 * iteration a new step is using the the given minimizer @M.
 *
 * The cost function object @J is assumed to provide the following functions:
 *
 * - \code J.gradient(const Vector& x)</tt> \endcode: Computes the gradient of
 * the cost function at <tt>x</tt>
 *
 * - \code J.Hessian(const Vector& x)\endcode: Computes the Hessian of the
 *  cost function at \code x \endcode. Return type is deducted using
 *  \code auto \endcode.
 *
 * The Minimizer object @J is assumed to provide the following functions:
 *
 * - \code M.step(Vector& dx, Vector& x, Vector& g, Matrix& H, CostFunction &J)
 *  \endcode: Computes the next iteration step * \code dx \endcode from the
 *  current position \code x \endcode from the gradient \code g \endcode and
 *  the Hessian \code H \endcode.
 *
 * \param CostFunction The cost function
 * \param Minimizer The minimizer object
 * \param x0 The start vector.
 * \param The solution vector. Used to store the current vector in the iteration
 * process.
 * \param max_iter The maximum number of iterations before the iterations aborts
 * \param tol The value of the convergence criterion of the cost function
 * required for convergence.
 *
 * \return Error code indicating success or failure of the minimization process.
 */
template
<
typename Real,
typename Vector,
typename CostFunction,
typename Minimizer
>
int minimize( CostFunction &J,
              Minimizer &M,
              const Vector &x0,
              Vector       &xi,
              unsigned int max_iter,
              Real         tol )
{
    bool converged     = false;
    unsigned int iter  = 0;

    Vector dx;
    xi = x0;

    // Minimization loop.
    while (!converged && (iter < max_iter))
    {
        auto g =  J.gradient(xi);
        auto H =  J.Hessian(xi);
        dx = M.step(xi, g, H, J );

        // Check for convergence.
        if (J.criterion(xi, dx) < tol)
            converged = true;

        xi += dx;
        iter++;
    }
    return 0;
}

}

#endif // OPTIMIZATION_MINIMIZE
