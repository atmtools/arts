// ----------------- //
//   MAP Base Class  //
// ----------------- //

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
MAPBase<ForwardModel, MatrixType, SaType, SeType>
::MAPBase(ForwardModel     &F_,
          const VectorType &xa_,
          const SaType     &Sa_,
          const SeType     &Se_)
    : m(F_.n), n(F_.m), F(F_), xa(xa_), yi_cached(), y_ptr(nullptr), Sa(Sa_), Se(Se_)
{
    // Nothing to do here.
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::cost_function(const VectorType &x,
                const VectorType &y,
                const VectorType &yi)
    -> RealType
{
    VectorType dy = y - yi;
    VectorType dx = xa - x;
    return dot(dy, inv(Se) * dy) + dot(dx, inv(Sa) * dx);
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::cost_function(const VectorType &x)
    -> RealType
{
    yi_cached = evaluate(x);
    cache_valid = true;

    VectorType dy = yi_cached - *y_ptr;
    VectorType dx = xa - x;
    return dot(dy, inv(Se) * dy) + dot(dx, inv(Sa) * dx);
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::cost_x(const VectorType &x)
    -> RealType
{
    VectorType dx = (xa - x);
    return dot(dx, inv(Sa) * dx);
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::cost_y(const VectorType &y,
         const VectorType &yi)
    -> RealType
{
    VectorType dy = y - yi;
    return dot(dy, inv(Se) * dy);
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::evaluate(const VectorType& x)
    -> GradientType
{
    try
    {
        return F.evaluate(x);
    }
    catch(...)
    {
        throw ForwardModelEvaluationException{};
    }
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::evaluate_cached(const VectorType& x)
    -> GradientType
{
    if (cache_valid)
    {
        return yi_cached;
    }
    else
    {
        try
        {
            return F.evaluate(x);
        }
        catch(...)
        {
            throw ForwardModelEvaluationException{};
        }
    }
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::Jacobian(const VectorType& x)
    -> JacobianType
{
    try
    {
        return F.Jacobian(x);
    }
    catch(...)
    {
        throw ForwardModelEvaluationException{};
    }
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType>
::gain_matrix(const VectorType &x)
    -> MatrixType
{
    auto &&K = Jacobian(x);
    MatrixType tmp = transp(K) * inv(Se);
    MatrixType G = inv(tmp * K + inv(Sa)) * tmp;
    return G;
}

// ----------------- //
//   Standard Form   //
// ----------------- //

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::STANDARD>
::MAP( ForwardModel &F_,
       const VectorType   &xa_,
       const SaType &Sa_,
       const SeType &Se_ )
    : Base(F_, xa_, Sa_, Se_), cost(-1.0), cost_x(-1.0), cost_y(-1.0)
{
    // Nothing to do here.
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
template<typename Minimizer, template <LogType> class Log>
auto MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::STANDARD>
::compute(VectorType       &x,
          const VectorType &y,
          Minimizer M,
          int verbosity)
    -> int
{

    Log<LogType::MAP> log(verbosity);
    log.init(Formulation::STANDARD, M);

    y_ptr = &y;
    x = xa;
    VectorType yi = evaluate(x);
    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    bool converged = false;
    iterations     = 0;

    while (iterations < M.get_maximum_iterations())
    {
        auto &&K = Jacobian(x);
        auto tmp = transp(K) * inv(Se);

        // Compute Hessian and transform.
        auto H  = tmp * K + inv(Sa);

        // Compute gradient and transform.
        VectorType g  = tmp * (yi - y) + inv(Sa) * (x - xa);

        if ((g.norm() / n) < M.get_tolerance())
        {
            converged = true;
            break;
        }

        cache_valid = false;
        dx = M.step(x, g, H, (*this));

        x += dx;
        yi = evaluate_cached(x);
        iterations++;

        cost_x = this->Base::cost_x(x);
        cost_y = this->Base::cost_y(y, yi);
        cost   = cost_x + cost_y;

        log.step(iterations, cost, cost_x, cost_y, M);
    }

    log.finalize(converged, iterations, cost, cost_x, cost_y);

    return 0;
}

// --------------- //
//     N-form      //
// --------------- //

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::NFORM>
::MAP( ForwardModel &F_,
       const VectorType   &xa_,
       const SaType &Sa_,
       const SeType &Se_ )
    : Base(F_, xa_, Sa_, Se_), cost(-1.0), cost_x(-1.0), cost_y(-1.0)
{
    // Nothing to do here.
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
template<typename Minimizer, template <LogType> class Log>
auto MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::NFORM>
::compute(VectorType       &x,
          const VectorType &y,
          Minimizer M,
          int verbosity)
    -> int
{

    Log<LogType::MAP> log(verbosity);

    y_ptr = &y;
    x = xa;
    auto &&yi = evaluate(x);
    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    bool converged = false;
    iterations = 0;

    while (iterations < M.get_maximum_iterations())
    {
        auto &&K  = Jacobian(x);
        auto tmp = transp(K) * inv(Se);

        // Compute true gradient for convergence test.
        VectorType g  = tmp * (yi - y) + inv(Sa) * (x - xa);

        if ((g.norm() / n) < M.get_tolerance())
        {
            converged = true;
            break;
        }

        // Compute Hessian and transform.
        auto H  = tmp * K + inv(Sa);

        // Compute gradient and transform.
        g = tmp * (y - yi + (K * (x - xa)));

        cache_valid = false;
        dx = M.step(xa, g, H, (*this));

        x = xa - dx;
        yi = evaluate_cached(x);
        iterations++;

        cost_x = this->Base::cost_x(x);
        cost_y = this->Base::cost_y(y, yi);
        cost   = cost_x + cost_y;

        log.step(iterations, cost, cost_x, cost_y, M);
    }

    log.finalize(converged, iterations, cost, cost_x, cost_y);

    return 0;
}

// --------------- //
//     M-form      //
// --------------- //

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::MFORM>
::MAP( ForwardModel &F_,
       const VectorType   &xa_,
       const SaType &Sa_,
       const SeType &Se_ )
    : Base(F_, xa_, Sa_, Se_), cost(-1.0), cost_x(-1.0), cost_y(-1.0)
{
    // Nothing to do here.
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
auto MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::MFORM>
::gain_matrix(const VectorType &x)
    -> MatrixType
{
    auto &&K = Jacobian(x);
    MatrixType SaKT = Sa * transp(K);
    MatrixType G = SaKT * inv(K * SaKT + Se);
    return G;
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType
>
template<typename Minimizer, template <LogType> class Log>
auto MAP<ForwardModel, MatrixType, SaType, SeType, Formulation::MFORM>
::compute(VectorType       &x,
          const VectorType &y,
          Minimizer M,
          int verbosity)
    -> int
{
    Log<LogType::MAP> log(verbosity);

    y_ptr = &y;
    x = xa;
    auto &&yi = evaluate(x);
    VectorType dx, yold;

    bool converged = false;
    iterations = 0;

    while (iterations < M.get_maximum_iterations())
    {
        auto &&K = Jacobian(x);
        auto tmp = Sa * transp(K);

        // Compute Hessian.
        auto H   = Se + K * tmp;

        // Compute gradient.
        VectorType g  = y - yi + K * (x - xa);

        cache_valid = false;
        dx = M.step(xa, g, H, (*this));

        x = xa - tmp * dx;

        yold = yi;
        yi = evaluate_cached(x);
        VectorType dy = yi - yold;
        VectorType r = Se * H * Se * dy;

        if ((dot(dy, r) / m) < M.get_tolerance())
        {
            converged = true;
            break;
        }
        iterations++;
    }
    return 0;
}
