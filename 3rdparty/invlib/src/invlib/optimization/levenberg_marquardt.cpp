
    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
LevenbergMarquardt<RealType, DampingMatrix, Solver>
::LevenbergMarquardt(const DampingMatrix &D_,
                     Solver solver)
    : current_cost(0.0), tolerance(1e-5), lambda(4.0), lambda_maximum(100.0),
      lambda_increase(2.0), lambda_decrease(3.0), lambda_threshold(1.0),
      lambda_constraint(std::numeric_limits<RealType>::min()),
      maximum_iterations(100), maximum_trials(100), step_count(0),
      stop_reason(LMStopReason::None), D(D_), s(solver)
{
    // Nothing to do here.
}

// ------------------------- //
//    Getters and Setters    //
// ------------------------- //

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_maximum_iterations() const
    -> unsigned int
{
    return maximum_iterations;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_maximum_iterations(unsigned int maximum_iterations_)
{
    maximum_iterations = maximum_iterations_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_maximum_trials() const
    -> unsigned int
{
    return maximum_trials;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_maximum_trials(unsigned int maximum_trials_)
{
    if (maximum_trials_ == 0) {
        throw std::invalid_argument("Levenberg-Marquardt maximum trials must be positive.");
    }
    maximum_trials = maximum_trials_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_tolerance() const
    -> RealType
{
    if (lambda > lambda_constraint) {
        // The caller uses a strict comparison. Even a rounded zero state
        // change must not pass while the damping constraint is unsatisfied.
        return 0.0;
    } else {
        return tolerance;
    }
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_tolerance(RealType tolerance_)
{
    tolerance = tolerance_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_lambda() const
    -> RealType
{
    return lambda;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_lambda(RealType lambda_)
{
    lambda = lambda_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_lambda_maximum() const
    -> RealType
{
    return lambda_maximum;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_lambda_maximum(RealType lambda_maximum_)
{
    lambda_maximum = lambda_maximum_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_lambda_decrease() const
    -> RealType
{
    return lambda_decrease;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_lambda_decrease(RealType lambda_decrease_)
{
    lambda_decrease = lambda_decrease_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_lambda_increase() const
    -> RealType
{
    return lambda_increase;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_lambda_increase(RealType lambda_increase_)
{
    lambda_increase = lambda_increase_;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_lambda_threshold() const
    -> RealType
{
    return lambda_threshold;
}

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_lambda_threshold(RealType lambda_threshold_)
{
    lambda_threshold = lambda_threshold_;
}

template
<
    typename RealType,
    typename DampingMatrix,
    typename Solver
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::get_lambda_constraint() const
    -> RealType
{
    return lambda_constraint;
}

template
<
    typename RealType,
    typename DampingMatrix,
    typename Solver
>
void LevenbergMarquardt<RealType, DampingMatrix, Solver>
::set_lambda_constraint(RealType lambda_constraint_)
{
    lambda_constraint = lambda_constraint_;
}

// --------------------------- //
//  Perform Minimization Step  //
// --------------------------- //

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
template
<
typename VectorType,
typename MatrixType,
typename CostFunction
>
auto LevenbergMarquardt<RealType, DampingMatrix, Solver>
::step(const VectorType &x,
       const VectorType &g,
       const MatrixType &B,
       CostFunction     &J)
    -> VectorType
{
    if (stop_iteration()) {
        throw std::logic_error("Levenberg-Marquardt step requested after termination.");
    }
    if (step_count == 0) {
        current_cost = J.cost_function(x);
    }
    if (!std::isfinite(current_cost)) {
        stop_reason = LMStopReason::NumericalFailure;
        throw std::runtime_error("Levenberg-Marquardt current cost is not finite.");
    }

    // Generic objectives supply matching cost, gradient and Hessian. MAP
    // supplies half-cost normal equations but reports full squared costs.
    const RealType cost_scale = [&] {
        if constexpr (requires { J.model_cost_scale(); }) {
            return J.model_cost_scale();
        } else {
            return RealType{1};
        }
    }();
    const auto finite = [](const VectorType &v) {
        for (decltype(v.rows()) i = 0; i < v.rows(); ++i) {
            if (!std::isfinite(v(i))) return false;
        }
        return true;
    };
    if (!finite(x) || !finite(g)) {
        stop_reason = LMStopReason::NumericalFailure;
        throw std::runtime_error("Levenberg-Marquardt state or gradient is not finite.");
    }

    VectorType zero(x);
    zero.scale(0.0);
    bool zero_gradient = true;
    for (decltype(g.rows()) i = 0; i < g.rows(); ++i) {
        zero_gradient = zero_gradient && g(i) == 0.0;
    }
    if (zero_gradient) {
        lambda = 0.0;
        stop_reason = LMStopReason::Stationary;
        ++step_count;
        return zero;
    }

    unsigned int trials = 0;
    const auto solve = [&](const auto &matrix) -> VectorType {
        if (trials == maximum_trials) {
            stop_reason = LMStopReason::TrialLimit;
            throw std::runtime_error(
                "Levenberg-Marquardt trial limit reached after "
                + std::to_string(trials) + " solves in one step.");
        }
        ++trials;
        VectorType result;
        try {
            result = -1.0 * s.solve(matrix, g);
        } catch (...) {
            stop_reason = LMStopReason::LinearSolverFailure;
            std::throw_with_nested(std::runtime_error(
                "Linear System Solution Error in Levenberg-Marquardt Method."));
        }
        if (!finite(result)) {
            stop_reason = LMStopReason::NumericalFailure;
            throw std::runtime_error("Levenberg-Marquardt linear solution is not finite.");
        }
        return result;
    };
    const auto predicted_reduction = [&](const VectorType &step) {
        return cost_scale * (-invlib::dot(g, step) - 0.5 * invlib::dot(step, B * step));
    };
    const auto roundoff = [](RealType old_cost, RealType new_cost) {
        return (32.0 * std::numeric_limits<RealType>::epsilon())
               * std::max(std::abs(old_cost), std::abs(new_cost));
    };
    bool first_step = true;
    bool stationarity_checked = false;

    while (true) {
        VectorType dx = solve(B + lambda * D);
        VectorType xnew = x + dx;
        const RealType new_cost = finite(xnew)
            ? J.cost_function(xnew, lambda < lambda_maximum)
            : std::numeric_limits<RealType>::infinity();
        const RealType predicted = predicted_reduction(dx);
        const RealType actual = current_cost - new_cost;
        const RealType noise = roundoff(current_cost, new_cost);

        // Do this BEFORE forming a ratio: a rounded zero (or a tiny random
        // decrease) divided by a vanishing prediction cannot prove progress.
        const bool ambiguous = std::isfinite(new_cost) && std::isfinite(predicted)
            && (std::abs(predicted) <= noise || std::abs(actual) <= noise);
        if (ambiguous && !stationarity_checked) {
            stationarity_checked = true;
            // A strongly damped step can be tiny far from the solution. Only
            // the undamped normal equations measure remaining state error.
            VectorType gn = lambda == 0.0 ? dx : solve(B);
            const RealType decrement = -invlib::dot(g, gn);
            const RealType gn_prediction = predicted_reduction(gn);
            if (std::isfinite(decrement) && decrement >= 0.0
                && decrement / static_cast<RealType>(x.rows()) < tolerance
                && std::isfinite(gn_prediction) && gn_prediction >= 0.0
                && gn_prediction <= roundoff(current_cost, current_cost)) {
                VectorType candidate = x + gn;
                const RealType candidate_cost = lambda == 0.0 ? new_cost
                    : (finite(candidate) ? J.cost_function(candidate, true)
                                         : std::numeric_limits<RealType>::infinity());
                if (std::isfinite(candidate_cost)
                    && std::abs(candidate_cost - current_cost)
                           <= roundoff(current_cost, candidate_cost)) {
                    current_cost = candidate_cost;
                    lambda = 0.0;
                    stop_reason = LMStopReason::Stationary;
                    ++step_count;
                    return gn;
                }
            }
        }

        // Require resolved, finite descent. In particular, NaN must never
        // escape the retry loop as an implicitly accepted trial.
        if (!ambiguous && std::isfinite(new_cost) && std::isfinite(predicted)
            && std::isfinite(actual) && predicted > 0.0 && actual > 0.0) {
            const RealType ratio = actual / predicted;
            if (std::isfinite(ratio) && ratio >= 0.5) {
                if (ratio > 0.75 && first_step) {
                    const RealType decreased = lambda / lambda_decrease;
                    lambda = decreased >= lambda_threshold ? decreased : 0.0;
                }
                current_cost = new_cost;
                ++step_count;
                return dx;
            }
        }

        // A failed trial is never applied, regardless of the reduction sign.
        // Keep lambda physical: maximum+1 is not representable at large maxima.
        if (lambda >= lambda_maximum) {
            stop_reason = LMStopReason::DampingLimit;
            ++step_count;
            return zero;
        }
        const RealType previous_lambda = lambda;
        if (lambda < lambda_threshold) {
            lambda = lambda_threshold;
        } else if (lambda >= lambda_maximum / lambda_increase) {
            lambda = lambda_maximum;
        } else {
            lambda *= lambda_increase;
        }
        if (!std::isfinite(lambda) || lambda <= previous_lambda) {
            stop_reason = LMStopReason::DampingStalled;
            throw std::runtime_error(
                "Levenberg-Marquardt damping did not increase to a finite value "
                "after a rejected trial; check the damping threshold and increase factor.");
        }
        first_step = false;
    }
}
