
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
      maximum_iterations(100), step_count(0),D(D_), s(solver)
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
::get_tolerance() const
    -> RealType
{
    return tolerance;
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
    if (step_count == 0)
        current_cost = J.cost_function(x);

    VectorType dx;

    RealType new_cost = 0.0;
    RealType c = -1.0;
    while (c < 0.0)
    {
        // Compute step.
        auto C = B + lambda * D;
        dx = -1.0 * s.solve(C, g);
        VectorType xnew = x + dx;

        // Compute model accuracy.
        new_cost = J.cost_function(xnew);
        c = (new_cost - current_cost) / (dot(g,dx) + 0.5 * dot(dx, B * dx));

        if (c > 0.75)
        {
            if (lambda >= (lambda_threshold * lambda_decrease))
                lambda /= lambda_decrease;
            else
                lambda = 0.0;

            current_cost = new_cost;
        }
        if (c < 0.2)
        {
            if (lambda < lambda_threshold)
                lambda = lambda_threshold;
            else
            {
                if (lambda < lambda_maximum)
                {
                    lambda *= lambda_increase;
                    if (lambda > lambda_maximum)
                        lambda = lambda_maximum;
                }
                else
                {
                    lambda = lambda_maximum + 1.0;
                    break;

                }
            }
        }
    }
    current_cost = new_cost;
    step_count++;

    return dx;
}
