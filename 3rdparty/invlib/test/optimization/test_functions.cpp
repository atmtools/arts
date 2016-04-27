#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Optimization, Test Functions"
#endif

#include <boost/test/included/unit_test.hpp>
#include <iostream>

#include "invlib/algebra.h"
#include "invlib/optimization.h"

#include "optimization/test_functions.h"
#include "utility.h"
#include "test_types.h"

using namespace invlib;

// Test non-linear optimization using the RandomPowerFunction test function
// which has a minimum at zero and uses the current function value as abortion
// criterion. After the minimization the value of the cost function should be
// lower than the chosen convergence tolerance.
template
<
typename T
>
void random_powers_test(unsigned int n)
{
    using RealType   = typename T::RealType;
    using VectorType = typename T::VectorType;
    using MatrixType = typename T::MatrixType;
    using Identity   = MatrixIdentity<MatrixType>;

    VectorType x0 = random<VectorType>(n);
    VectorType dx; dx.resize(n);

    RandomPowerFunction<RealType, VectorType, MatrixType> J(n);
    Identity I{};
    LevenbergMarquardt<RealType, Identity> LM(I);
    GaussNewton<RealType> GN{};
    GradientDescent<RealType> GD{};

    VectorType x;
    minimize(J, LM, x0, x, 1000, EPS);
    BOOST_TEST((J.cost_function(x) < EPS), "J(x) = " << J.cost_function(x));

    minimize(J, GN, x0, x, 1000, EPS);
    BOOST_TEST((J.cost_function(x) < EPS), "J(x) = " << J.cost_function(x));

    minimize(J, GD, x0, x, 1000, EPS);
    BOOST_TEST((J.cost_function(x) < EPS), "J(x) = " << J.cost_function(x));
}

// Test non-linear optimization using the Sphere test function
// which has a minimum at zero and uses the current function value as abortion
// criterion. After the minimization the value of the cost function should be
// lower than the chosen convergence tolerance.
template
<
typename T
>
void sphere_function_test(unsigned int n)
{
    using RealType   = typename T::RealType;
    using VectorType = typename T::VectorType;
    using MatrixType = typename T::MatrixType;
    using Identity   = MatrixIdentity<MatrixType>;

    VectorType x0 = random<VectorType>(n);
    VectorType dx; dx.resize(n);

    SphereFunction<RealType, VectorType, MatrixType> J(n);

    Identity I{};
    LevenbergMarquardt<RealType, Identity> LM(I);
    GaussNewton<RealType> GN{};
    GradientDescent<RealType> GD{};

    VectorType x;
    minimize(J, LM, x0, x, 1000, EPS);
    BOOST_TEST((J.cost_function(x) < EPS), "J(x) = " << J.cost_function(x));

    minimize(J, GN, x0, x, 1000, EPS);
    BOOST_TEST((J.cost_function(x) < EPS), "J(x) = " << J.cost_function(x));

    minimize(J, GD, x0, x, 1000, EPS);
    BOOST_TEST((J.cost_function(x) < EPS), "J(x) = " << J.cost_function(x));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(random_powers, T, matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 1 + rand() % 100;
        random_powers_test<T>(n);
        sphere_function_test<T>(n);
    }
}
