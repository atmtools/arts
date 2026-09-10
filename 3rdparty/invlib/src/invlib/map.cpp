// A minimizer's explicit stop outcome takes precedence over a small returned
// step. Rejected trials commonly return zero and must not imply convergence.
template<typename Minimizer, typename RealType>
bool minimizer_converged(Minimizer &M, RealType criterion)
{
    if (M.stop_iteration()) {
        if constexpr (requires { M.converged(); }) {
            return M.converged();
        }
        return false;
    }
    return std::isfinite(criterion) && criterion < M.get_tolerance();
}

// Only the built-in state-step criteria can skip the new simulated
// measurement. Custom criteria, including derived overrides, stay conservative.
template<typename Criterion>
constexpr bool criterion_needs_measurement = true;

template<typename VectorType>
constexpr bool criterion_needs_measurement<Rodgers530<VectorType>> = false;

template<typename VectorType>
constexpr bool criterion_needs_measurement<Rodgers531<VectorType>> = false;

// ----------------- //
//   MAP Base Class  //
// ----------------- //

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
::MAPBase(ForwardModel     &F_,
          const VectorType &xa_,
          const SaType     &Sa_,
          const SeType     &Se_)
    : m(F_.m), n(F_.n), F(F_), xa(xa_), y_ptr(nullptr), Sa(Sa_), Se(Se_),
      evaluate_time(duration<double>::zero()), Jacobian_time(duration<double>::zero())
{
    // Nothing to do here.
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
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
typename SeType,
typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
::cost_function(const VectorType &x,
                bool robust)
    -> RealType
{
    try {
        // evaluate() returns owned storage: reuse it for the residual instead
        // of copying the full measurement vector a second time per LM trial.
        VectorType dy = evaluate(x);
        dy.subtract(*y_ptr);
        VectorType dx = xa - x;
        return dot(dy, inv(Se) * dy) + dot(dx, inv(Sa) * dx);
    } catch (...) {
        if (robust) {
            return std::numeric_limits<RealType>::max();
        }
        throw;
    }
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
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
typename SeType,
typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
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
typename SeType,
typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
::evaluate(const VectorType& x)
    -> MeasurementVectorType
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
        std::throw_with_nested(std::runtime_error("Forward Model Evaluation Error"));
    }
}

template
<
typename ForwardModel,
typename MatrixType,
typename SaType,
typename SeType,
typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
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
        std::throw_with_nested(std::runtime_error("Forward Model Evaluation Error"));
    }
}

template
<
    typename ForwardModel,
    typename MatrixType,
    typename SaType,
    typename SeType,
    typename VectorType
>
auto MAPBase<ForwardModel, MatrixType, SaType, SeType, VectorType>
::gain_matrix(const VectorType &x)
    -> MatrixType
{
    VectorType y; y.resize(m);
    auto && K = Jacobian(x, y);
    MatrixType tmp = transp(K) * inv(Se);
    MatrixType G = inv(tmp * K + inv(Sa)) * tmp;
    return G;
}

template<typename ReferenceType, typename OtherType>
auto inline operator *(std::reference_wrapper<ReferenceType> & A, const OtherType & B)
    -> decltype(remove_reference_wrapper(A) * B)
{
    return remove_reference_wrapper(A) * B;
}

// ----------------- //
//   Standard Form   //
// ----------------- //

template
<
    typename ForwardModel,
    typename MatrixType,
    typename SaType,
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::STANDARD, ConvergenceCriterion>
::MAP( ForwardModel &F_,
       const VectorType   &xa_,
       const SaType &Sa_,
       const SeType &Se_ )
    : Base(F_, xa_, Sa_, Se_), cost(-1.0), cost_x(-1.0), cost_y(-1.0),
      iterations(0)
{
    // Nothing to do here.
}

template
<
    typename ForwardModel,
    typename MatrixType,
    typename SaType,
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
template<typename Minimizer, template <LogType> class Log, typename ... LogParams>
auto MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::STANDARD, ConvergenceCriterion>
::compute(VectorType       &x,
          const VectorType &y,
          Minimizer M,
          LogParams && ... log_params)
    -> int
{

    Log<LogType::MAP> log(log_params ...);
    Formulation f = Formulation::STANDARD;
    log.init(F, xa, Sa, Se, y, M, f);

    auto t1 = std::chrono::steady_clock::now();

    y_ptr = &y;

    if (x.rows() != n) {
        x = xa;
    }

    MeasurementVectorType yi; yi.resize(m);
    JacobianType K = Jacobian(x, yi);
    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    bool converged = false;
    iterations     = 0;

    RealType conv = NAN;
    log.step(iterations, cost, cost_x, cost_y, conv, M);

    ConvergenceCriterion<VectorType> criterion{};
    criterion(x, yi, y, x, K, Sa, Se);

    while ((iterations < M.get_maximum_iterations())
           && !M.stop_iteration()
           && !converged)
    {

        // Compute next step.
        auto tmp = transp(K) * inv(Se);
        auto H  = tmp * K + inv(Sa);
        VectorType g  = tmp * (yi - y) + inv(Sa) * (x - xa);
        dx = M.step(x, g, H, (*this));
        x += dx;

        // State-step criteria need no new forward value. A continuing
        // iteration obtains both value and derivative from one Jacobian call.
        constexpr bool needs_measurement = criterion_needs_measurement<decltype(criterion)>;
        if constexpr (needs_measurement) yi = evaluate(x);
        conv = criterion(x, yi, y, g, K, Sa, Se);

        converged = minimizer_converged(M, conv);
        if (!converged && !M.stop_iteration()
            && iterations + 1 < M.get_maximum_iterations()) {
            K = Jacobian(x, yi);
        } else if constexpr (!needs_measurement) {
            yi = evaluate(x);
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

    if (converged) {
        return 0;
    } else {
        return 1;
    }
}

// --------------- //
//     N-form      //
// --------------- //

template
<
    typename ForwardModel,
    typename MatrixType,
    typename SaType,
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::NFORM, ConvergenceCriterion>
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
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
template<typename Minimizer, template <LogType> class Log, typename ... LogParams>
auto MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::NFORM, ConvergenceCriterion>
::compute(VectorType       &x,
          const VectorType &y,
          Minimizer M,
          LogParams && ... log_params)
    -> int
{

    Log<LogType::MAP> log(log_params ...);
    Formulation f = Formulation::NFORM;
    log.init(F, xa, Sa, Se, y, M, f);

    auto t1 = std::chrono::steady_clock::now();

    y_ptr = &y;
    if (x.rows() != n) {
        x = xa;
    }
    MeasurementVectorType yi; yi.resize(m);
    JacobianType K = Jacobian(x, yi);

    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    RealType conv = NAN;
    log.step(iterations, cost, cost_x, cost_y, conv, M);

    ConvergenceCriterion<VectorType> criterion{};
    criterion(x, yi, y, x, K, Sa, Se);

    bool converged = false;
    iterations = 0;

    while ((iterations < M.get_maximum_iterations())
           && !M.stop_iteration()
           && !converged)
    {
        auto tmp = transp(K) * inv(Se);

        // Compute step.
        VectorType g = tmp * (y - yi + (K * (x - xa)));
        auto H  = tmp * K + inv(Sa);
        dx = M.step(xa, g, H, (*this));
        x = xa - dx;

        // State-step criteria need no new forward value. A continuing
        // iteration obtains both value and derivative from one Jacobian call.
        constexpr bool needs_measurement = criterion_needs_measurement<decltype(criterion)>;
        if constexpr (needs_measurement) yi = evaluate(x);
        conv = criterion(x, yi, y, g, K, Sa, Se);

        converged = minimizer_converged(M, conv);
        if (!converged && !M.stop_iteration()
            && iterations + 1 < M.get_maximum_iterations()) {
            K = Jacobian(x, yi);
        } else if constexpr (!needs_measurement) {
            yi = evaluate(x);
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

    if (converged) {
        return 0;
    } else {
        return 1;
    }
}

// --------------- //
//     M-form      //
// --------------- //

template
<
    typename ForwardModel,
    typename MatrixType,
    typename SaType,
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::MFORM, ConvergenceCriterion>
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
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
auto MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::MFORM, ConvergenceCriterion>
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
    typename SeType,
    typename VectorType,
    template <class> class ConvergenceCriterion
>
template<typename Minimizer, template <LogType> class Log, typename ... LogParams>
auto MAP<ForwardModel, MatrixType, SaType, SeType, VectorType, Formulation::MFORM, ConvergenceCriterion>
::compute(VectorType       &x,
          const VectorType &y,
          Minimizer M,
          LogParams && ... log_params)
    -> int
{
    Log<LogType::MAP> log(log_params ...);
    auto t1 = std::chrono::steady_clock::now();
    Formulation f = Formulation::MFORM;
    log.init(F, xa, Sa, Se, y, M, f);

    y_ptr = &y;
    if (x.rows() != n) {
        x = xa;
    }
    MeasurementVectorType yi; yi.resize(m);
    JacobianType K = Jacobian(x, yi);
    VectorType dx;

    cost_x = this->Base::cost_x(x);
    cost_y = this->Base::cost_y(y, yi);
    cost   = cost_x + cost_y;

    bool converged = false;
    iterations = 0;

    RealType conv = NAN;
    log.step(iterations, cost, cost_x, cost_y, conv, M);

    ConvergenceCriterion<VectorType> criterion{};
    criterion(x, yi, y, x, K, Sa, Se);

    while ((iterations < M.get_maximum_iterations())
           && !M.stop_iteration()
           && !converged)
    {
        // Compute step.
        constexpr bool dense_system = [] {
            if constexpr (requires { std::remove_cvref_t<Minimizer>::dense_measurement_system; })
                return std::remove_cvref_t<Minimizer>::dense_measurement_system;
            else return false;
        }();
        VectorType g = y - yi + K * (x - xa);
        if constexpr (dense_system) {
            // Materialize the covariance/Jacobian product once and reuse it
            // for both assembly and mapping the solution into state space.
            MatrixType KT = transp(K);
            MatrixType tmp = Sa * KT;
            MatrixType H = [&]() -> MatrixType {
                if constexpr (requires { K.multiply_add(tmp, Se); }) {
                    return K.multiply_add(tmp, Se);
                } else {
                    MatrixType result = K * tmp;
                    result += Se;
                    return result;
                }
            }();
            dx = M.step(xa, g, H, (*this));
            x = xa - tmp * dx;
        } else {
            auto tmp = Sa * transp(K);
            auto H = Se + K * tmp;
            dx = M.step(xa, g, H, (*this));
            x = xa - tmp * dx;
        }

        // State-step criteria need no new forward value. A continuing
        // iteration obtains both value and derivative from one Jacobian call.
        constexpr bool needs_measurement = criterion_needs_measurement<decltype(criterion)>;
        if constexpr (needs_measurement) yi = evaluate(x);
        conv = criterion(x, yi, y, g, K, Sa, Se);

        converged = minimizer_converged(M, conv);
        if (!converged && !M.stop_iteration()
            && iterations + 1 < M.get_maximum_iterations()) {
            K = Jacobian(x, yi);
        } else if constexpr (!needs_measurement) {
            yi = evaluate(x);
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

    if (converged) {
        return 0;
    } else {
        return 1;
    }
}
