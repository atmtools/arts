#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Algebra, Solvers"
#endif

#include <boost/test/included/unit_test.hpp>
#include "invlib/algebra.h"
#include "invlib/algebra/solvers.h"
#include "utility.h"
#include "test_types.h"

using namespace invlib;

// Test solvers by computing A * inv(A) * v for a random vector v and a
// random positive definite matrix A. The resulting vector should be equal
// to v up to the precision of the underlying solver.
template
<
typename MatrixType
>
void solver_test(unsigned int n)
{

    using VectorType = typename MatrixType::VectorType;
    using Preconditioner = NormalizeDiagonal<MatrixType>;

    auto A  = random_positive_definite<MatrixType>(n);
    auto v = random<VectorType>(n);
    VectorType w; w.resize(n);

    Standard std{};
    ConjugateGradient cg(1e-20);
    Preconditioner pre(A);
    PreconditionedSolver<Standard, Preconditioner>          std_n(std, pre);
    PreconditionedSolver<ConjugateGradient, Preconditioner> cg_n(cg, pre);

    w = A * std.solve(A, v);
    double error = maximum_error(v, w);
    BOOST_TEST((error < EPS), "Standard solver error: " << error);

    w = A * cg.solve(A, v);
    error = maximum_error(v, w);
    BOOST_TEST((error < EPS), "CG solver error: " << error);

    w = A * std_n.solve(A, v);
    error = maximum_error(v, w);
    BOOST_TEST((error < EPS), "Standard preconditioned solver error: " << error);

    w = A * cg_n.solve(A, v);
    error = maximum_error(v, w);
    BOOST_TEST((error < EPS), "CG preconditioned solver error: " << error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(solver,
                              T,
                              matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 1 + rand() % 100;
        solver_test<T>(n);
    }
}
