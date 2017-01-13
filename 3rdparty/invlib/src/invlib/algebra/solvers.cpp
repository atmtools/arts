// ------------------ //
//  Standard Solver   //
// ------------------ //

    template<typename VectorType, typename MatrixType>
auto Standard::solve(const MatrixType &A,const VectorType &v)
    -> VectorType
{
    VectorType x = A.solve(v);
    return x;
}

// -------------------------  //
//  Conjugate Gradient Solver //
// -------------------------  //

ConjugateGradient::ConjugateGradient(double tol, int verbosity_)
: verbosity(verbosity_), tolerance(tol)
{
    // Nothing to do here.
}

template
<
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
auto ConjugateGradient::solve(const MatrixType &A,
        const VectorType &v)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    Log<LogType::SOL_CG> log(verbosity);

    RealType tol, alpha, beta, rnorm, vnorm;
    VectorType x, r, p, xnew, rnew, pnew;

    x = 0.0 * v;
    r = A * x - v;
    p = -1.0 * r;
    vnorm = v.norm();
    rnorm = r.norm();

    log.init(tolerance, rnorm, vnorm);
    int i = 0;
    while (rnorm / vnorm > tolerance)
    {
        alpha = dot(r, r) / dot(p, A * p);
        xnew  = x + alpha *     p;
        rnew  = r + alpha * A * p;
        beta  = dot(rnew, rnew) / dot(r, r);
        pnew  = beta * p - rnew;

        x = xnew;
        r = rnew; rnorm = r.norm();
        p = pnew;

        i++;
        if (i % 10 == 0) {
            log.step(i, rnorm / vnorm);
        }
    }

    log.finalize(i);
    return x;
}

// ----------------------------------------- //
//  Preconditioned Conjugate Gradient Solver //
// ----------------------------------------- //

template<typename F>
PreconditionedConjugateGradient<F, true>::PreconditionedConjugateGradient(
    const F &f_,
    double tolerance_,
    int verbosity_)
    : f(f_), verbosity(verbosity_), tolerance(tolerance_)
{
    // Nothing to do here.
}

template <typename F>
template
<
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
auto PreconditionedConjugateGradient<F, true>::solve(const MatrixType &A,
                                                     const VectorType &v)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    Log<LogType::SOL_CG> log(verbosity);

    RealType tol, alpha, beta, rnorm, r0;
    VectorType x, y, r, p, xnew, ynew, rnew, pnew;

    x = 0.0 * v;
    r = A * x - v;
    y = f(r);
    p = -1.0 * y;
    rnorm = r.norm();
    r0    = rnorm;

    log.init(tolerance, rnorm, v.norm());
    int i = 0;
    while (rnorm / r0 > tolerance)
    {
        alpha = dot(r, y) / dot(p, A * p);
        xnew  = x + alpha *     p;
        rnew  = r + alpha * A * p;
        ynew  = f(rnew);
        beta  = dot(rnew, ynew) / dot(r, y);
        pnew  = beta * p - ynew;

        x = xnew;
        r = rnew; rnorm = r.norm();
        p = pnew;
        y = ynew;

        i++;
        if (i % 10 == 0) {
            log.step(i, rnorm / r0);
        }
    }

    log.finalize(i);
    return x;
}

template<typename F>
PreconditionedConjugateGradient<F, false>::PreconditionedConjugateGradient(
    double tolerance_,
    int verbosity_)
    : verbosity(verbosity_), tolerance(tolerance_)
{
    // Nothing to do here.
}

template <typename F>
template
<
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
auto PreconditionedConjugateGradient<F, false>::solve(const MatrixType &A,
                                                      const VectorType &v)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    F f(A);
    Log<LogType::SOL_CG> log(verbosity);

    RealType tol, alpha, beta, rnorm, r0;
    VectorType x, y, r, p, xnew, ynew, rnew, pnew;

    x = 0.0 * v;
    r = A * x - v;
    y = f(r);
    p = -1.0 * y;
    rnorm = r.norm();
    r0    = rnorm;

    log.init(tolerance, rnorm, v.norm());
    int i = 0;
    while (rnorm / r0 > tolerance)
    {
        alpha = dot(r, y) / dot(p, A * p);
        xnew  = x + alpha *     p;
        rnew  = r + alpha * A * p;
        ynew  = f(rnew);
        beta  = dot(rnew, ynew) / dot(r, y);
        pnew  = beta * p - ynew;

        x = xnew;
        r = rnew; rnorm = r.norm();
        p = pnew;
        y = ynew;

        i++;
        if (i % 10 == 0) {
            log.step(i, rnorm / r0);
        }
    }

    log.finalize(i);
    return x;
}
