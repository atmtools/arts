#ifndef TEST_OPTIMIZATION_TEST_FUNCTIONS_H
#define TEST_OPTIMIZATION_TEST_FUNCTIONS_H

#include "invlib/algebra.h"
#include <iostream>

/**
 * \brief Sphere function
 *
 * Quadratic cost function to test optimization methods given by
 * \f[
 *    J(\vec{x}) = \sum_i x_i^2
 * \f]
 * Implementes the interface required by the generic minimization function
 * minimize().
 */
template
<
typename R,
typename V,
typename M
>
class SphereFunction
{
public:

    using Real   = R;
    using Vector = V;
    using Matrix = M;

    /*
     * \brief Create sphere function in @n dimensions.
     */
    SphereFunction( unsigned int n_ )
        : n(n_) {}

    /*
     * \brief Evaluate sphere function a @x.
     */
    double cost_function(const Vector &x) const
    {
        return dot(x,x);
    }

    /*
     * \brief Convergence criterion is step size.
     */
    Real criterion( const Vector &x,
                    const Vector &dx ) const
    {
        return dx.norm();
    }

    /*
     * \brief Gradient of the sphere function.
     */
    Vector gradient( const Vector &x ) const
    {
        Vector g; g.resize(n);

        for (unsigned int i = 0; i < n; i++)
        {
            g(i) = 2.0 * x(i);
        }

        return g;
    }

    /*
     * \brief Hessian of the sphere function.
     */
    Matrix Hessian(const Vector &x) const
    {
        Matrix H{}; H.resize(n,n);
        for (unsigned int i = 0; i < n; i++)
        {
            for (unsigned int j = 0; j < n; j++)
            {
                H(i, j) = 0.0;
            }
            H(i, i) = 2.0;
        }
        return H;
    }

    const unsigned int n;
};

/**
 * \brief Random powers.
 *
 * Cost function of the form
 * \f[
 *    J(\vec{x}) = \sum_i x_i^{e_i}
 * \f]
 * where \f$e_i \in [2,4,6] \f$. Has a minimum at \f$\vec{0}\f$.
 * Implementes the interface required by the generic minimization function
 * minimize().
 */
template
<
typename R,
typename V,
typename M
>
class RandomPowerFunction
{
public:

    using Real   = R;
    using Vector = V;
    using Matrix = M;

    /*
     * \brief Create sphere function in @n dimensions.
     */
    RandomPowerFunction( unsigned int n_ )
        : n(n_), e(n_) {

        for (unsigned int i = 0; i < n; i++)
        {
            e[i] = 2 * (rand() % 4 + 1);
        }
    }

    /*
     * \brief Evaluate sphere function a @x.
     */
    Real cost_function(const Vector &x) const
    {
        Real res = 0.0;

        for (unsigned int i = 0; i < n; i++)
        {
            res += std::pow(x(i), e[i]);
        }
        return res;
    }

    /*
     * \brief Convergence criterion is step size.
     */
    Real criterion( const Vector &x,
                    const Vector &dx ) const
    {
        return cost_function(x);
    }

    /*
     * \brief Gradient of the sphere function.
     */
    Vector gradient( const Vector &x ) const
    {
        Vector g; g.resize(n);

        for (unsigned int i = 0; i < n; i++)
        {
            g(i) = e[i] * pow(x(i), e[i]-1);
        }

        return g;
    }

    /*
     * \brief Hessian of the sphere function.
     */
    Matrix Hessian(const Vector &x) const
    {
        Matrix H{}; H.resize(n,n);
        for (unsigned int i = 0; i < n; i++)
        {
            for (unsigned int j = 0; j < n; j++)
            {
                H(i, j) = 0.0;
            }
            H(i, i) = e[i] * (e[i] - 1.0) * pow(x(i) , e[i] - 2);
        }
        return H;
    }

    const unsigned int n;
    std::vector<int> e;
};

template
<
typename Real,
typename Matrix,
typename Vector
>
class SumOfPowers
{

public:

    SumOfPowers( unsigned int n )
        : n_(n) {}

    inline Real m() {return n_;}
    inline Real n() {return n_;}

    Real criterion( const Vector &x,
                    const Vector &dx ) const
    {
        return dx.norm();
    }

    Vector step( const Vector &x ) const
    {
        Vector J(); J.resize(n_);
        Matrix H(); H.resize(n_, n_);


        for (unsigned int i = 0; i < n_; i++)
        {
            J(i) = ((double) (i + 2)) * pow(x(i), i+1);
            H(i, i) = ((double) (i + 1) * (i + 2)) * pow(x(i), i);
        }

        Vector dx = inv(H) * J;
        return dx;
    }

private:

    unsigned int n_;

};

#endif // TEST_OPTIMIZATION_TEST_FUNCTIONS_H
