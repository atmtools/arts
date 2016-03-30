#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "Forward Models, Linear"
#endif

#include <boost/test/included/unit_test.hpp>
#include <iostream>

#include "invlib/algebra.h"
#include "invlib/algebra/solvers.h"
#include "invlib/map.h"
#include "invlib/optimization.h"

#include "forward_models/linear.h"
#include "utility.h"
#include "test_types.h"

using namespace invlib;

// Use a random linear forward model to test the equivalence of the
// standard, n-form and m-form when using the Gauss-Newton optimizer
// and the standard form using Levenberg-Marquardt and Gauss-Newton
// optimization.
template
<
typename T
>
void linear_test(unsigned int n)
{
    using RealType   = typename T::RealType;
    using VectorType = typename T::VectorType;
    using MatrixType = typename T::MatrixType;
    using Id     = MatrixIdentity<MatrixType>;
    using Model  = Linear<MatrixType>;

    MatrixType Se = random_positive_definite<MatrixType>(n);
    MatrixType Sa = random_positive_definite<MatrixType>(n);
    VectorType xa = random<VectorType>(n);
    VectorType y  = random<VectorType>(n);

    Model F(n,n);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::STANDARD>
        std(F, xa, Sa, Se);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::NFORM>
        nform(F, xa, Sa, Se);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::MFORM>
        mform(F, xa, Sa, Se);

    // Test inversion using standard solver.
    Id I{};
    GaussNewton<RealType> gn{};
    gn.set_tolerance(1e-9); gn.set_maximum_iterations(1000);
    LevenbergMarquardt<RealType, Id> lm(I);
    lm.set_tolerance(1e-9); lm.set_maximum_iterations(1000);

    VectorType x_std_lm, x_std_gn, x_n_gn, x_m_gn;
    std.compute(x_std_lm, y, lm, 1);
    std.compute(x_std_gn, y, gn);
    nform.compute(x_n_gn, y, gn);
    mform.compute(x_m_gn, y, gn);

    RealType e1, e2, e3;
    e1 = maximum_error(x_std_lm, x_std_gn);
    e2 = maximum_error(x_std_gn, x_n_gn);
    e3 = maximum_error(x_std_gn, x_m_gn);

    BOOST_TEST((e1 < EPS), "Error STD - NFORM = " << e1);
    BOOST_TEST((e2 < EPS), "Error STD - MFORM = " << e2);
    BOOST_TEST((e3 < EPS), "Error STD - MFORM = " << e3);

    // Test inversion using CG solver.
    ConjugateGradient cg(1e-9);
    GaussNewton<RealType, ConjugateGradient> gn_cg(cg);
    gn_cg.set_tolerance(1e-9); gn_cg.set_maximum_iterations(1000);
    LevenbergMarquardt<RealType, Id, ConjugateGradient> lm_cg(I, cg);
    lm_cg.set_tolerance(1e-9); lm_cg.set_maximum_iterations(1000);

    std.compute(x_std_lm, y, lm_cg);
    std.compute(x_std_gn, y, gn_cg);
    nform.compute(x_n_gn, y, gn_cg);
    mform.compute(x_m_gn, y, gn_cg);

    e1 = maximum_error(x_std_lm, x_std_gn);
    e2 = maximum_error(x_std_gn, x_n_gn);
    e3 = maximum_error(x_std_gn, x_m_gn);

    BOOST_TEST((e1 < EPS), "Error STD - NFORM CG = " << e1);
    BOOST_TEST((e2 < EPS), "Error STD - MFORM CG = " << e2);
    BOOST_TEST((e3 < EPS), "Error STD - MFORM CG = " << e3);
}

// Same test as above but applying the NormalizingDiagonal preconditioner.
template
<
typename T
>
void linear_test_transformed(unsigned int n)
{
    using RealType       = typename T::RealType;
    using VectorType     = typename T::VectorType;
    using MatrixType     = typename T::MatrixType;
    using Id             = MatrixIdentity<MatrixType>;
    using Model          = Linear<MatrixType>;
    using Preconditioner = NormalizeDiagonal<MatrixType>;
    using StdPre         = PreconditionedSolver<Standard, Preconditioner>;
    using CGPre          = PreconditionedSolver<ConjugateGradient, Preconditioner>;

    MatrixType Se = random_positive_definite<MatrixType>(n);
    MatrixType Sa = random_positive_definite<MatrixType>(n);
    VectorType xa = random<VectorType>(n);
    VectorType y  = random<VectorType>(n);

    Model F(n,n);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::STANDARD>
        std(F, xa, Sa, Se);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::NFORM>
        nform(F, xa, Sa, Se);
    MAP<Model, MatrixType, MatrixType, MatrixType, Formulation::MFORM>
        mform(F, xa, Sa, Se);

    // Test inversion using standard, preconditioned solver.
    Preconditioner pre(Sa);
    StdPre std_pre(Standard(), pre);
    Id I{};
    GaussNewton<RealType> gn{};
    GaussNewton<RealType, StdPre> gn_pre(std_pre);
    gn_pre.set_tolerance(1e-9); gn_pre.set_maximum_iterations(1000);
    LevenbergMarquardt<RealType, Id, StdPre> lm_pre(I, std_pre);
    lm_pre.set_tolerance(1e-9); lm_pre.set_maximum_iterations(1000);

    VectorType x_ref, x_std_lm, x_std_gn, x_n_gn, x_m_gn;
    std.compute(x_ref, y, gn);
    std.compute(x_std_lm, y, lm_pre);
    std.compute(x_std_gn, y, gn_pre);
    nform.compute(x_n_gn, y, gn_pre);
    mform.compute(x_m_gn, y, gn_pre);

    RealType e1, e2, e3, e4;
    e1 = maximum_error(x_std_lm, x_ref);
    e2 = maximum_error(x_std_gn, x_ref);
    e3 = maximum_error(x_n_gn, x_ref);
    e4 = maximum_error(x_m_gn, x_ref);

    BOOST_TEST((e1 < EPS), "Error STD - NFORM = " << e1);
    BOOST_TEST((e2 < EPS), "Error STD - MFORM = " << e2);
    BOOST_TEST((e3 < EPS), "Error STD - MFORM = " << e3);
    BOOST_TEST((e4 < EPS), "Error STD - MFORM = " << e4);

    // Test inversion using CG solver.
    ConjugateGradient cg(1e-9);
    CGPre cg_pre(cg, pre);
    GaussNewton<RealType, CGPre> gn_cg_pre(cg_pre);
    gn_cg_pre.set_tolerance(1e-9); gn_cg_pre.set_maximum_iterations(1000);
    LevenbergMarquardt<RealType, Id, CGPre> lm_cg_pre(I, cg_pre);
    lm_cg_pre.set_tolerance(1e-9); lm_cg_pre.set_maximum_iterations(1000);

    std.compute(x_std_lm, y, lm_cg_pre);
    std.compute(x_std_gn, y, gn_cg_pre);
    nform.compute(x_n_gn, y, gn_cg_pre);
    mform.compute(x_m_gn, y, gn_cg_pre);

    e1 = maximum_error(x_std_lm, x_ref);
    e2 = maximum_error(x_std_gn, x_ref);
    e3 = maximum_error(x_n_gn, x_ref);
    e4 = maximum_error(x_m_gn, x_ref);

    BOOST_TEST((e1 < EPS), "Error STD - NFORM = " << e1);
    BOOST_TEST((e2 < EPS), "Error STD - MFORM = " << e2);
    BOOST_TEST((e3 < EPS), "Error STD - MFORM = " << e3);
    BOOST_TEST((e4 < EPS), "Error STD - MFORM = " << e4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(linear, T, matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 1 + rand() % 50;
        linear_test<T>(n);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(linear_transformed, T, matrix_types)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < ntests; i++)
    {
        unsigned int n = 1 + rand() % 50;
        linear_test_transformed<T>(n);
    }
}
