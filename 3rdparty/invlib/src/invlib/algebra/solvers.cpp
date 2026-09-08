namespace solver_detail
{

inline void validate_cg_settings(double tolerance, int max_iterations)
{
    if (!std::isfinite(tolerance) || tolerance <= 0.0)
    {
        throw std::invalid_argument(
            "Conjugate gradient tolerance must be finite and positive.");
    }
    if (max_iterations <= 0)
    {
        throw std::invalid_argument(
            "Conjugate gradient max_iterations must be positive.");
    }
}

template<typename RealType>
inline void require_finite(RealType value, const char *quantity)
{
    if (!std::isfinite(value))
    {
        throw std::runtime_error(std::string("Conjugate gradient requires a finite ")
                                 + quantity + ". Check the linear system and scaling.");
    }
}

template<typename RealType>
inline void require_positive(RealType value, const char *quantity)
{
    require_finite(value, quantity);
    if (value <= 0.0)
    {
        throw std::runtime_error(std::string("Conjugate gradient encountered nonpositive ")
                                 + quantity
                                 + ". Check positive definiteness and scaling.");
    }
}

inline void check_iteration_limit(int iteration, int max_iterations)
{
    if (iteration >= max_iterations)
    {
        throw std::runtime_error("Conjugate gradient iteration limit ("
                                 + std::to_string(max_iterations)
                                 + ") reached before convergence. Check conditioning, "
                                   "tolerance, or increase max_iterations.");
    }
}

} // namespace solver_detail

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
    solver_detail::validate_cg_settings(tolerance, max_iterations);
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
    solver_detail::require_finite(vnorm, "right-hand side norm");
    auto x = settings.start_vector(v);
    solver_detail::require_finite(x.norm(), "initial solution norm");
    r = A * x - v;
    p = -1.0 * r;
    rnorm = r.norm();
    solver_detail::require_finite(rnorm, "residual norm");

    log.init(tolerance, rnorm, vnorm);
    int i = 0;
    // The safety limit belongs to the solver, independently of the custom
    // convergence policy. An exact solution also stops fixed-step policies
    // before a subsequent iteration could divide zero by zero.
    while (rnorm != 0.0 && !settings.converged(r, v))
    {
        solver_detail::check_iteration_limit(i, max_iterations);
        const RealType rr = invlib::dot(r, r);
        solver_detail::require_positive(rr, "squared residual norm");
        ap = A * p;
        const RealType curvature = invlib::dot(p, ap);
        solver_detail::require_positive(curvature, "curvature (p^T A p)");
        alpha = rr / curvature;
        solver_detail::require_positive(alpha, "step length");
        xnew  = x + alpha *     p;
        rnew  = r + alpha * ap;
        solver_detail::require_finite(xnew.norm(), "solution norm");
        rnorm = rnew.norm();
        solver_detail::require_finite(rnorm, "residual norm");
        beta  = invlib::dot(rnew, rnew) / rr;
        solver_detail::require_finite(beta, "direction coefficient");
        pnew  = beta * p - rnew;
        solver_detail::require_finite(pnew.norm(), "search direction norm");

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
                                   int verbosity, int max_iterations)
    -> VectorType
{
    using RealType = typename VectorType::RealType;

    Log<LogType::SOL_CG> log(verbosity);

    RealType alpha, beta, rnorm, r0, vnorm;
    VectorType x, y, r, p, ap, xnew, ynew, rnew, pnew;

    vnorm = v.norm();
    require_finite(vnorm, "right-hand side norm");
    x = 0.0 * v;
    r = A * x - v;
    rnorm = r.norm();
    require_finite(rnorm, "residual norm");
    r0    = rnorm;

    log.init(tolerance, rnorm, vnorm);
    int i = 0;
    if (rnorm == 0.0 || rnorm / r0 <= tolerance)
    {
        log.finalize(i);
        return x;
    }
    y = f(r);
    require_finite(y.norm(), "preconditioned residual norm");
    p = -1.0 * y;

    while (rnorm != 0.0 && rnorm / r0 > tolerance)
    {
        check_iteration_limit(i, max_iterations);
        const RealType ry = invlib::dot(r, y);
        require_positive(ry, "preconditioned residual product (r^T M r)");
        ap = A * p;
        const RealType curvature = invlib::dot(p, ap);
        require_positive(curvature, "curvature (p^T A p)");
        alpha = ry / curvature;
        require_positive(alpha, "step length");
        xnew  = x + alpha *     p;
        rnew  = r + alpha * ap;
        require_finite(xnew.norm(), "solution norm");
        rnorm = rnew.norm();
        require_finite(rnorm, "residual norm");
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
        require_finite(ynew.norm(), "preconditioned residual norm");
        const RealType rynew = invlib::dot(rnew, ynew);
        require_positive(rynew, "preconditioned residual product (r^T M r)");
        beta  = rynew / ry;
        require_finite(beta, "direction coefficient");
        pnew  = beta * p - ynew;
        require_finite(pnew.norm(), "search direction norm");

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
    solver_detail::validate_cg_settings(tolerance, max_iterations);
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
        f, A, v, tolerance, verbosity, max_iterations);
}

template<typename F>
inline PreconditionedConjugateGradient<F, false>::PreconditionedConjugateGradient(
    double tolerance_,
    int verbosity_,
    int max_iterations_)
    : verbosity(verbosity_), tolerance(tolerance_), max_iterations(max_iterations_)
{
    solver_detail::validate_cg_settings(tolerance, max_iterations);
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
        f, A, v, tolerance, verbosity, max_iterations);
}
