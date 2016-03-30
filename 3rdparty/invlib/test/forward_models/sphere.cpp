#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Forward Models, Sphere"
#endif

#include <boost/test/included/unit_test.hpp>
#include <iostream>

#include "invlib/algebra.h"
#include "invlib/map.h"
#include "invlib/optimization.h"

#include "forward_models/sphere.h"
#include "utility.h"
#include "test_types.h"

using namespace invlib;

// Use the sphere function forward model to test the equivalence of the
// standard, n-form and m-form when using the Gauss-Newton optimizer.
template
<
typename T
>
void sphere_test(unsigned int n)
{
    using RealType   = typename T::RealType;
    using VectorType = typename T::VectorType;
    using MatrixType = typename T::MatrixType;
    using Model      = Sphere<MatrixType>;

    MatrixType Se = random_positive_definite<MatrixType>(1);
    MatrixType Sa = random_positive_definite<MatrixType>(n);
    VectorType xa = random<VectorType>(n);
    VectorType y  = random<VectorType>(1);

    Model F(n);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::STANDARD>
        std(F, xa, Sa, Se);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::NFORM>
        nform(F, xa, Sa, Se);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::MFORM>
        mform(F, xa, Sa, Se);

    GaussNewton<RealType> GN{};
    GN.set_tolerance(1e-15); GN.set_maximum_iterations(1000);

    VectorType x_std, x_n, x_m;
    std.compute(x_std, y, GN);
    nform.compute(x_n, y, GN);
    mform.compute(x_m, y, GN);

    RealType e1, e2;
    e1 = maximum_error(x_std, x_m);
    e2 = maximum_error(x_std, x_n);

    BOOST_TEST((e1 < EPS), "Error STD - NFORM = " << e1);
    BOOST_TEST((e2 < EPS), "Error STD - MFORM =" << e2);

    // Test inversion using CG solver.

    ConjugateGradient cg(1e-15);
    GaussNewton<RealType, ConjugateGradient> GN_CG(cg);
    GN_CG.set_tolerance(1e-15); GN_CG.set_maximum_iterations(1000);

    std.compute(x_std, y, GN_CG);
    nform.compute(x_n, y, GN_CG);
    mform.compute(x_m, y, GN_CG);

    e1 = maximum_error(x_std, x_m);
    e2 = maximum_error(x_std, x_n);

    BOOST_TEST((e1 < EPS), "Error STD - NFORM CG = " << e1);
    BOOST_TEST((e2 < EPS), "Error STD - MFORM CG = " << e2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sphere, T, matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 10;
        sphere_test<T>(n);
    }
}
