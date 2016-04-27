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
    : m(F_.m), n(F_.n), F(F_), xa(xa_), y_ptr(nullptr), Sa(Sa_), Se(Se_)
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
    VectorType y = evaluate(x);
    VectorType dy = y - *y_ptr;
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
    -> FMVectorType
{
    try
    {
        auto t1 = std::chrono::steady_clock::now();
        auto y = F.evaluate(x);
        auto t2 = std::chrono::steady_clock::now();
        evaluate_time += duration_cast<duration<double>>(t2 - t1);

        return y;
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
::Jacobian(const VectorType& x, VectorType &y)
    -> JacobianType
{
    try
    {
        auto t1 = std::chrono::steady_clock::now();
        JacobianType J = F.Jacobian(x, y);
        auto t2 = std::chrono::steady_clock::now();
        Jacobian_time += duration_cast<duration<double>>(t2 - t1);

        return J;
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
    VectorType y; y.resize(m);
    auto &&K = Jacobian(x, y);
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

    auto t1 = std::chrono::steady_clock::now();

    y_ptr = &y;
    x = xa;
    FMVectorType yi; yi.resize(m);
    JacobianType K = Jacobian(x, yi);
    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    bool converged = false;
    iterations     = 0;

    while ((iterations < M.get_maximum_iterations()) && !converged)
    {

        // Compute next step.
        auto tmp = transp(K) * inv(Se);
        auto H  = tmp * K + inv(Sa);
        VectorType g  = tmp * (yi - y) + inv(Sa) * (x - xa);
        dx = M.step(x, g, H, (*this));
        x += dx;

        // Check for convergence.
        RealType conv = dot(dx, H * dx) / n;
        if (conv < M.get_tolerance())
        {
            converged = true;
            yi = evaluate(x);
        }
        else
        {
            K = Jacobian(x, yi);
        }

        // Log output.
        iterations++;
        cost_x = this->Base::cost_x(x);
        cost_y = this->Base::cost_y(y, yi);
        cost   = cost_x + cost_y;
        log.step(iterations, cost, cost_x, cost_y, conv, M);
    }

    log.finalize(converged, iterations, cost, cost_x, cost_y);

    // Timing output.
    auto t2 = std::chrono::steady_clock::now();
    auto compute_time = duration_cast<duration<double>>(t2 - t1);
    log.time(compute_time.count(), evaluate_time.count(), Jacobian_time.count());

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
    auto t1 = std::chrono::steady_clock::now();

    y_ptr = &y;
    x = xa;
    FMVectorType yi; yi.resize(m);
    JacobianType K = Jacobian(x, yi);

    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    bool converged = false;
    iterations = 0;

    while (iterations < M.get_maximum_iterations() && !converged)
    {
        auto tmp = transp(K) * inv(Se);

        // Compute step.
        VectorType g = tmp * (y - yi + (K * (x - xa)));
        auto H  = tmp * K + inv(Sa);
        dx = M.step(xa, g, H, (*this));
        x = xa - dx;

        // Test for convergence.
        g  = tmp * (yi - y) + inv(Sa) * (x - xa);
        RealType conv = dot(dx, H * dx) / n;
        if (conv < M.get_tolerance())
        {
            converged = true;
            yi = evaluate(x);
        }
        else
        {
            K = Jacobian(x, yi);
        }

        // Log output.
        iterations++;
        cost_x = this->Base::cost_x(x);
        cost_y = this->Base::cost_y(y, yi);
        cost   = cost_x + cost_y;
        log.step(iterations, cost, cost_x, cost_y, conv, M);
    }

    log.finalize(converged, iterations, cost, cost_x, cost_y);

    // Timing output.
    auto t2 = std::chrono::steady_clock::now();
    auto compute_time = duration_cast<duration<double>>(t2 - t1);
    log.time(compute_time.count(), evaluate_time.count(), Jacobian_time.count());

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
    VectorType y; y.resize(m);
    auto &&K = Jacobian(x, y);
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
    auto t1 = std::chrono::steady_clock::now();

    y_ptr = &y;
    x = xa;
    FMVectorType yi; yi.resize(m);
    JacobianType K = Jacobian(x, yi);
    VectorType dx, yold;

    bool converged = false;
    iterations = 0;

    while (iterations < M.get_maximum_iterations() && !converged)
    {
        // Compute step.
        auto tmp = Sa * transp(K);
        auto H   = Se + K * tmp;
        VectorType g  = y - yi + K * (x - xa);
        dx = M.step(xa, g, H, (*this));
        x = xa - tmp * dx;

        // Check convergence.
        yold = yi;
        yi = evaluate(x);
        VectorType dy = yi - yold;
        VectorType r = inv(Se) * H * inv(Se) * dy;
        RealType conv = dot(dy, r) / m;
        if (conv < M.get_tolerance())
        {
            converged = true;
            break;
        }

        // Book keeping.
        iterations++;

        if (!converged)
            K = Jacobian(x , yi);
    }
    return 0;

    // Timing output.
    auto t2 = std::chrono::steady_clock::now();
    auto compute_time = duration_cast<duration<double>>(t2 - t1);
    log.time(compute_time.count(), evaluate_time.count(), Jacobian_time.count());
}
