// ------------------ //
//  Standard Solver   //
// ------------------ //

    template<typename VectorType, typename MatrixType>
inline auto Standard::solve(const MatrixType &A,const VectorType &v)
    -> VectorType
{
    VectorType x = A.solve(v);
    return x;
}

// -------------------------------- //
//  Setting Functors for CG Solvers //
// -------------------------------- //

inline CGDefaultSettings::CGDefaultSettings(double tolerance_)
    : tolerance(tolerance_)
{
    // Nothing to do here.
}

template<typename VectorType>
inline VectorType CGDefaultSettings::start_vector(const VectorType & v) const
{
    VectorType w = 0.0 * v;
    return w;
}

template<typename VectorType>
inline bool CGDefaultSettings::converged(const VectorType & r,
                                  const VectorType & v) const
{
    const auto rnorm = r.norm();
    const auto vnorm = v.norm();
    return rnorm == 0.0 || (vnorm > 0.0 && rnorm / vnorm < tolerance);
}

template<size_t maximum_steps>
inline CGStepLimit<maximum_steps>::CGStepLimit(double /* unused */)
    : steps(0)
{
    // Nothing to do here.
}

template<size_t maximum_steps>
template<typename VectorType>
inline VectorType CGStepLimit<maximum_steps>::start_vector(const VectorType & v)
{
    steps = 0;
    VectorType w = 0.0 * v;
    return w;
}

template<size_t maximum_steps>
template<typename VectorType>
inline bool CGStepLimit<maximum_steps>::converged(const VectorType & /*r*/,
                                           const VectorType & /*v*/)
{
    steps++;
    return (steps > maximum_steps);
}

template<typename VectorType, size_t maximum_steps>
inline CGContinued<VectorType, maximum_steps>::CGContinued(double /* unused */)
    : steps(0)
{
    // Nothing to do here.
}

template<typename VectorType, size_t maximum_steps>
inline VectorType & CGContinued<VectorType, maximum_steps>::start_vector(const VectorType & w)
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
inline bool CGContinued<VectorType, maximum_steps>::converged(const VectorType & /*r*/,
                                                              const VectorType & /*v*/)
{
    steps++;
    return (steps > maximum_steps);
}

// -------------------------  //
//  Conjugate Gradient Solver //
// -------------------------  //

template<typename CGSettings>
inline ConjugateGradient<CGSettings>::ConjugateGradient(double tol, int verbosity_,
                                                       int max_iterations_)
    : verbosity(verbosity_), tolerance(tol), max_iterations(max_iterations_), settings(tol)
{
    assert(std::isfinite(tolerance) && tolerance > 0.0);
    assert(max_iterations > 0);
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
inline auto ConjugateGradient<CGSettings>::solve(const MatrixType &A,
                                                 const VectorType &v)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    Log<LogType::SOL_CG> log(verbosity);

    RealType alpha, beta, rnorm, vnorm;
    VectorType r, p, ap, xnew, rnew, pnew;

    vnorm = v.norm();
    assert(std::isfinite(vnorm) && "right-hand side norm");
    auto x = settings.start_vector(v);
    assert(std::isfinite(x.norm()) && "initial solution norm");
    r = A * x - v;
    p = -1.0 * r;
    rnorm = r.norm();
    assert(std::isfinite(rnorm) && "residual norm");

    log.init(tolerance, rnorm, vnorm);
    int i = 0;
    // The safety limit belongs to the solver, independently of the custom
    // convergence policy. An exact solution also stops fixed-step policies
    // before a subsequent iteration could divide zero by zero.
    while (rnorm != 0.0 && !settings.converged(r, v))
    {
        if (i >= max_iterations) {
            if (iteration_limit_warning) iteration_limit_warning();
            break;
        }
        const RealType rr = invlib::dot(r, r);
        assert(std::isfinite(rr) && rr > 0.0 && "squared residual norm");
        ap = A * p;
        const RealType curvature = invlib::dot(p, ap);
        assert(std::isfinite(curvature) && curvature > 0.0 && "curvature (p^T A p)");
        alpha = rr / curvature;
        assert(std::isfinite(alpha) && alpha > 0.0 && "step length");
        xnew  = x + alpha *     p;
        rnew  = r + alpha * ap;
        assert(std::isfinite(xnew.norm()) && "solution norm");
        rnorm = rnew.norm();
        assert(std::isfinite(rnorm) && "residual norm");
        beta  = invlib::dot(rnew, rnew) / rr;
        assert(std::isfinite(beta) && "direction coefficient");
        pnew  = beta * p - rnew;
        assert(std::isfinite(pnew.norm()) && "search direction norm");

        x = xnew;
        r = rnew;
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

namespace solver_detail
{

template
<
    typename F,
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
inline auto solve_preconditioned_cg(const F &f, const MatrixType &A,
                                   const VectorType &v, double tolerance,
                                   int verbosity, int max_iterations,
                                   const std::function<void()>& iteration_limit_warning)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    Log<LogType::SOL_CG> log(verbosity);

    RealType alpha, beta, rnorm, r0, vnorm;
    VectorType x, y, r, p, ap, xnew, ynew, rnew, pnew;

    vnorm = v.norm();
    assert(std::isfinite(vnorm) && "right-hand side norm");
    x = 0.0 * v;
    r = A * x - v;
    rnorm = r.norm();
    assert(std::isfinite(rnorm) && "residual norm");
    r0    = rnorm;

    log.init(tolerance, rnorm, vnorm);
    int i = 0;
    if (rnorm == 0.0 || rnorm / r0 <= tolerance)
    {
        log.finalize(i);
        return x;
    }
    y = f(r);
    assert(std::isfinite(y.norm()) && "preconditioned residual norm");
    p = -1.0 * y;

    while (rnorm != 0.0 && rnorm / r0 > tolerance)
    {
        if (i >= max_iterations) {
            if (iteration_limit_warning) iteration_limit_warning();
            break;
        }
        const RealType ry = invlib::dot(r, y);
        assert(std::isfinite(ry) && ry > 0.0 && "preconditioned residual product (r^T M r)");
        ap = A * p;
        const RealType curvature = invlib::dot(p, ap);
        assert(std::isfinite(curvature) && curvature > 0.0 && "curvature (p^T A p)");
        alpha = ry / curvature;
        assert(std::isfinite(alpha) && alpha > 0.0 && "step length");
        xnew  = x + alpha *     p;
        rnew  = r + alpha * ap;
        assert(std::isfinite(xnew.norm()) && "solution norm");
        rnorm = rnew.norm();
        assert(std::isfinite(rnorm) && "residual norm");
        x = xnew;
        i++;
        if (i % 10 == 0) {
            log.step(i, rnorm / r0);
        }
        if (rnorm == 0.0 || rnorm / r0 <= tolerance)
        {
            break;
        }
        ynew  = f(rnew);
        assert(std::isfinite(ynew.norm()) && "preconditioned residual norm");
        const RealType rynew = invlib::dot(rnew, ynew);
        assert(std::isfinite(rynew) && rynew > 0.0 && "preconditioned residual product (r^T M r)");
        beta  = rynew / ry;
        assert(std::isfinite(beta) && "direction coefficient");
        pnew  = beta * p - ynew;
        assert(std::isfinite(pnew.norm()) && "search direction norm");

        r = rnew;
        p = pnew;
        y = ynew;
    }

    log.finalize(i);
    return x;
}

} // namespace solver_detail

template<typename F>
inline PreconditionedConjugateGradient<F, true>::PreconditionedConjugateGradient(
    const F &f_,
    double tolerance_,
    int verbosity_,
    int max_iterations_)
    : f(f_), verbosity(verbosity_), tolerance(tolerance_), max_iterations(max_iterations_)
{
    assert(std::isfinite(tolerance) && tolerance > 0.0);
    assert(max_iterations > 0);
}

template <typename F>
template
<
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
inline auto PreconditionedConjugateGradient<F, true>::solve(const MatrixType &A,
                                                            const VectorType &v)
    -> VectorType
{
    return solver_detail::solve_preconditioned_cg<F, VectorType, MatrixType, Log>(
        f, A, v, tolerance, verbosity, max_iterations, iteration_limit_warning);
}

template<typename F>
inline PreconditionedConjugateGradient<F, false>::PreconditionedConjugateGradient(
    double tolerance_,
    int verbosity_,
    int max_iterations_)
    : verbosity(verbosity_), tolerance(tolerance_), max_iterations(max_iterations_)
{
    assert(std::isfinite(tolerance) && tolerance > 0.0);
    assert(max_iterations > 0);
}

template <typename F>
template
<
    typename VectorType,
    typename MatrixType,
    template <LogType> class Log
>
inline auto PreconditionedConjugateGradient<F, false>::solve(const MatrixType &A,
                                                             const VectorType &v)
    -> VectorType
{
    F f(A);
    return solver_detail::solve_preconditioned_cg<F, VectorType, MatrixType, Log>(
        f, A, v, tolerance, verbosity, max_iterations, iteration_limit_warning);
}
