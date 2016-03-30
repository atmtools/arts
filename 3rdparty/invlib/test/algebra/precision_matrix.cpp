#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Algebra, Precision Matrix"
#endif

#include <boost/test/included/unit_test.hpp>
#include "invlib/algebra.h"
#include "invlib/algebra/precision_matrix.h"
#include "utility.h"
#include "test_types.h"
#include <iostream>

using namespace invlib;

template <typename T>
void foo(T t);

// Test behaviour of precision matrix. The PrecisionMatrix wrapper should
// make a matrix act like its inverse. This is tested below by comparing
// the precision matrix constructed from A to its inverse and vice versa.
template
<
typename MatrixType
>
void precision_test(unsigned int n)
{
    using VectorType = typename MatrixType::VectorType;

    MatrixType A  = random_positive_definite<MatrixType>(n);
    PrecisionMatrix<MatrixType> P(A);
    VectorType v = random<VectorType>(n);

    //foo(inv(P));
    //foo(P);
    MatrixType B = P * A;
    MatrixType C = inv(P) * inv(A);
    MatrixType I; I.resize(n, n); set_identity(I);
    double error = maximum_error(B, I);
    BOOST_TEST((error < EPS), "Deviation from identity:" << error);
    error = maximum_error(C, I);
    BOOST_TEST((error < EPS), "Deviations from identity:" << error);

    MatrixType D = P;
    MatrixType E = inv(A);
    error = maximum_error(D, E);
    BOOST_TEST((error < EPS), "Deviation from inv(A): " << error);

    VectorType w1 = inv(A) * v;
    VectorType w2 = P * v;
    error = maximum_error(w1, w2);
    BOOST_TEST((error < EPS), "Vector mult error: " << error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(precision,
                              T,
                              matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 1 + rand() % 100;
        precision_test<T>(4);
    }
}
