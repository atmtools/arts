
// ------------------------------- //
//  Constructors and Destructors   //
// ------------------------------- //

template
<
typename RealType
>
GradientDescent<RealType>::GradientDescent()
    : tolerance(1e-5), step_length(1e-2), scale(2.0), minimum_step_length(1e-12),
    maximum_iterations(1000), step_count(0)
{
    // Nothing to do here.
}

template
<
typename RealType
>
GradientDescent<RealType>::GradientDescent(RealType tolerance_,
                                           unsigned int maximum_iterations_,
                                           RealType initial_step,
                                           RealType scale_)
    : tolerance(tolerance_), step_length(initial_step), scale(scale_),
      minimum_step_length(1e-12), maximum_iterations(maximum_iterations_),
      step_count(0)
{
    // Nothing to do here.
}

// -------------------------- //
//    Getters and Setters     //
// -------------------------- //

template
<
typename RealType
>
auto GradientDescent<RealType>::get_maximum_iterations() const
    -> unsigned int
{
    return maximum_iterations;
}

template
<
typename RealType
>
void GradientDescent<RealType>::set_maximum_iterations(unsigned int n)
{
    maximum_iterations = n;
}

template
<
typename RealType
>
auto GradientDescent<RealType>::get_tolerance() const
    -> RealType
{
    return tolerance;
}

template
<
typename RealType
>
void GradientDescent<RealType>::set_tolerance(RealType tolerance_)
{
    tolerance = tolerance_;
}

template
<
typename RealType
>
auto GradientDescent<RealType>::get_initial_step() const
    -> RealType
{
    return step_length;
}

template
<
typename RealType
>
void GradientDescent<RealType>::set_initial_step(RealType initial_step_)
{
    step_length = initial_step_;
}

template
<
typename RealType
>
auto GradientDescent<RealType>::get_scale() const
    -> RealType
{
    return scale;
}

template
<
typename RealType
>
void GradientDescent<RealType>::set_scale(RealType scale_)
{
    scale = scale_;
}

// --------------------------- //
//  Perform Minimization Step  //
// --------------------------- //

template
<
typename RealType
>
template
<
typename VectorType,
typename MatrixType,
typename CostFunction
>
auto GradientDescent<RealType>::step(const VectorType &x,
                                     const VectorType &g,
                                     const MatrixType &B,
                                     CostFunction &J)
    -> VectorType
{
    if (step_count == 0)
        current_cost = J.cost_function(x);

    bool step_found = false;
    RealType gnorm = g.norm();
    RealType current_step_length = step_length;

    // Try current step length.
    VectorType dx = - step_length / gnorm * g;
    VectorType x_new  = x + dx;
    RealType new_cost;

    try
    {
        new_cost = J.cost_function(x_new);
    }
    catch(...)
    {
        new_cost = current_cost;
    }

    // Increase step length if first
    // step length is successful.
    if (new_cost < current_cost)
    {
        step_length *= scale;
        current_cost = new_cost;
        step_found   = true;
    }

    while (!step_found && step_length > minimum_step_length)
    {
        step_length /= scale;
        dx       = - step_length / gnorm * g;
        x_new    = x + dx;

        try
        {
            new_cost = J.cost_function(x_new);
        }
        catch (...)
        {
            new_cost = current_cost;
        }


        if (new_cost < current_cost)
        {
            current_cost = new_cost;
            step_found = true;
        }
    }

    step_count++;
    return dx;
}
