#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Algebra, Identities"
#endif

#include <boost/test/included/unit_test.hpp>
#include "invlib/algebra.h"
#include "utility.h"
#include "test_types.h"

using namespace invlib;

template
<
typename Matrix
>
void identities_test( unsigned int m,
                      unsigned int k,
                      unsigned int n )
{

    auto A = random<Matrix>(m, k);
    auto B = random<Matrix>(m, k);
    auto C = random<Matrix>(k, n);
    auto D = random<Matrix>(k, n);

    // Addition.

    Matrix R1 = A + B;
    Matrix R2 = B + A;
    double error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    R1 = (A + B) + B;
    R2 = B + (A + B);
    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    // Subtraction.

    R1 = A - B;
    R2 = -1.0 * (B - A);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    R1 = (B - A) + B;
    R2 = B - (A - B);

    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    // Multiplication.

    R1 = transp(A * C);
    R2 = transp(C) * transp(A);
    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    R1 = (A + B) * (C + D);
    R2 = A*C + A*D + B*C + B*D;
    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    R1 = (A - B) * (C - D);
    R2 = A*C - A*D - B*C + B*D;
    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    R1 = transp((A + B) * (C + D));
    R2 = transp(A*C) + transp(A*D) + transp(B*C) + transp(B*D);
    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    R2 = transp(C)*transp(A) + transp(D)*transp(A)
        + transp(C)*transp(B) + transp(D)*transp(B);
    error = maximum_error(R1, R2);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);

    // Inversion.

    Matrix E = random_positive_definite<Matrix>(m);
    R1 = inv(E) * E;
    R2 = E * inv(E);
    set_identity(E);
    error = maximum_error(R1, E);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);
    error = maximum_error(R2, E);
    BOOST_TEST((error < EPS), "k = " << k << ", m = " << m << ", n = " << n << ", error = " << error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(identities,
                              T,
                              matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int k = 1 + rand() % 100;
        unsigned int m = 1 + rand() % 100;
        unsigned int n = 1 + rand() % 100;
        identities_test<T>(m, k, n);
    }
}
