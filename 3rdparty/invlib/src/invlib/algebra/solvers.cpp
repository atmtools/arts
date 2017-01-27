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

// -------------------------------- //
//  Setting Functors for CG Solvers //
// -------------------------------- //

CGDefaultSettings::CGDefaultSettings(double tolerance_)
    : tolerance(tolerance_)
{
    // Nothing to do here.
}

template<typename VectorType>
VectorType CGDefaultSettings::start_vector(const VectorType & v) const
{
    VectorType w = 0.0 * v;
    return w;
}

template<typename VectorType>
bool CGDefaultSettings::converged(const VectorType & r,
                                  const VectorType & v) const
{
    return ((r.norm() / v.norm()) < tolerance);
}

template<size_t maximum_steps>
CGStepLimit<maximum_steps>::CGStepLimit(double /* unused */)
    : steps(0)
{
    // Nothing to do here.
}

template<size_t maximum_steps>
template<typename VectorType>
VectorType CGStepLimit<maximum_steps>::start_vector(const VectorType & v)
{
    steps = 0;
    VectorType w = 0.0 * v;
    return w;
}

template<size_t maximum_steps>
template<typename VectorType>
bool CGStepLimit<maximum_steps>::converged(const VectorType & /*r*/,
                                           const VectorType & /*v*/)
{
    steps++;
    return (steps > maximum_steps);
}

template<typename VectorType, size_t maximum_steps>
CGContinued<VectorType, maximum_steps>::CGContinued(double /* unused */)
    : steps(0)
{
    // Nothing to do here.
}

template<typename VectorType, size_t maximum_steps>
VectorType & CGContinued<VectorType, maximum_steps>::start_vector(const VectorType & w)
{
    if (steps == 0)
    {
        v = w;
        v.scale(0.0);
    }

    steps = 0;
    return v;
}

template<typename VectorType, size_t maximum_steps>
bool CGContinued<VectorType, maximum_steps>::converged(const VectorType & /*r*/,
                                           const VectorType & /*v*/)
{
    steps++;
    return (steps > maximum_steps);
}

// -------------------------  //
//  Conjugate Gradient Solver //
// -------------------------  //

template<typename CGSettings>
ConjugateGradient<CGSettings>::ConjugateGradient(double tol, int verbosity_)
    : verbosity(verbosity_), tolerance(tol), settings(tol)
{
    // Nothing to do here.
}

template
<
    typename CGSettings
>
template
<
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
auto ConjugateGradient<CGSettings>::solve(const MatrixType &A,
                                          const VectorType &v)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    Log<LogType::SOL_CG> log(verbosity);

    RealType tol, alpha, beta, rnorm, vnorm;
    VectorType r, p, xnew, rnew, pnew;

    auto x = settings.start_vector(v);
    r = A * x - v;
    p = -1.0 * r;
    vnorm = v.norm();
    rnorm = r.norm();

    log.init(tolerance, rnorm, vnorm);
    int i = 0;
    while (!settings.converged(r, v))
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
