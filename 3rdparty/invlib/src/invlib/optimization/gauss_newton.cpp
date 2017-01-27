// ------------------------------- //
//  Constructors and Destructors   //
// ------------------------------- //

template
<
typename RealType,
typename Solver
>
GaussNewton<RealType, Solver>::GaussNewton(Solver solver_)
    : tolerance(1e-5), maximum_iterations(1000), solver(solver_)
{
    // Nothing to do here.
}

template
<
typename RealType,
typename Solver
>
GaussNewton<RealType, Solver>::GaussNewton(RealType tolerance_,
                                           unsigned int maximum_iterations_,
                                           Solver solver_)
    : tolerance(tolerance_), maximum_iterations(maximum_iterations_), solver(solver_)
{
    // Nothing to do here.
}

// -------------------------- //
//    Getters and Setters     //
// -------------------------- //

template
<
typename RealType,
typename Solver
>
auto GaussNewton<RealType, Solver>::get_maximum_iterations() const
    -> unsigned int
{
    return maximum_iterations;
}

template
<
typename RealType,
typename Solver
>
void GaussNewton<RealType, Solver>::set_maximum_iterations(unsigned int n)
{
    maximum_iterations = n;
}

template
<
typename RealType,
typename Solver
>
auto GaussNewton<RealType, Solver>::get_tolerance() const
    -> RealType
{
    return tolerance;
}

template
<
typename RealType,
typename Solver
>
void GaussNewton<RealType, Solver>::set_tolerance(RealType tolerance_)
{
    tolerance = tolerance_;
}

// --------------------------- //
//  Perform Minimization Step  //
// --------------------------- //

template
<
typename RealType,
typename Solver
>
template
<
typename VectorType,
typename MatrixType,
typename CostFunction
>
auto GaussNewton<RealType, Solver> ::step(const VectorType &,
                                          const VectorType &g,
                                          const MatrixType &B,
                                          const CostFunction &)
    -> VectorType
{
    try
    {
        VectorType dx = -1.0 * solver.solve(B, g);
        return dx;
    }
    catch (...)
    {
        std::runtime_error(
            "Linear System Solution Error in Gauss-Newton Method."
            );
        return VectorType();
    }
}
