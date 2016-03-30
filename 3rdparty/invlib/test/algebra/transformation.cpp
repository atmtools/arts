#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Algebra, Transformation"
#endif

#include <boost/test/included/unit_test.hpp>
#include "invlib/algebra.h"
#include "utility.h"
#include "test_types.h"

namespace invlib
{

template
<
typename T
>
void transformation_test(unsigned int n)
{

    using RealType   = typename T::RealType;
    using VectorType = typename T::VectorType;
    using MatrixType = typename T::MatrixType;

    auto A = random_positive_definite<MatrixType>(n);
    auto B = random_positive_definite<MatrixType>(n);
    auto v = random<VectorType>(n);

    RealType error;

    // Identity Transformation
    Identity I{};
    VectorType w1 = I.apply(A * B) * I.apply(v);
    VectorType w2 = A * B * v;
    error = maximum_error(w1, w2);
    BOOST_TEST((error < EPS),"maximum_error(w1, w2) = " << error);

    // NormalizeDiagonal Transform
    NormalizeDiagonal<MatrixType> t(A);
    w1 = t.apply(inv(t.apply(A)) * t.apply(v));
    w2 = inv(A) * v;
    error = maximum_error(w1, w2);
    BOOST_TEST((error < EPS),"maximum_error(w1, w2) = " << error);

}

BOOST_AUTO_TEST_CASE_TEMPLATE(transformation, T, matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 1 + rand() % 100;
        transformation_test<T>(n);
    }
}

}
