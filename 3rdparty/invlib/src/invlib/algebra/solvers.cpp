// ------------------ //
//  Standard Solver   //
// ------------------ //

template<typename VectorType, typename MatrixType>
auto Standard::solve(const MatrixType &A,const VectorType &v)
    -> VectorType
{
    return A.solve(v);
}

// -------------------------  //
//  Conjugate Gradient Solver //
// -------------------------  //

ConjugateGradient::ConjugateGradient(double tol)
    : tolerance(tol)
{
    // Nothing to do here.
}

template<typename VectorType, typename MatrixType>
auto ConjugateGradient::solve(const MatrixType &A,
                               const VectorType &v)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    unsigned int n = v.rows();
    RealType tol, alpha, beta, rnorm;
    VectorType x, r, p, xnew, rnew, pnew;

    x = v;
    r = A * x - v;
    p = -1.0 * r;

    int i = 0;
    while (r.norm() > tolerance)
    {
        alpha = dot(r, r) / dot(p, A * p);
        xnew  = x + alpha *     p;
        rnew  = r + alpha * A * p;
        beta  = dot(rnew, rnew) / dot(r, r);
        pnew  = beta * p - rnew;

        x = xnew;
        r = rnew;
        p = pnew;
    }

    return x;
}

// -----------------------  //
//  Preconditioned Solver   //
// -----------------------  //

template<typename Solver, typename Transformation>
PreconditionedSolver<Solver, Transformation>::
PreconditionedSolver(Solver s, Transformation t)
    : solver(s), transformation(t)
{
    // Nothin to do here.
}

template<typename Solver, typename Transformation>
    template<typename VectorType, typename MatrixType>
auto PreconditionedSolver<Solver, Transformation>::solve(const MatrixType &A,
                                                         const VectorType &v)
    -> VectorType
{
    VectorType vt = transformation.apply(v);
    auto At       = transformation.apply(A);
    VectorType x  = transformation.apply(solver.solve(At, vt));
    return x;
}
