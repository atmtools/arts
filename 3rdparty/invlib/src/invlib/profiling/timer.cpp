// -------------------------- //
//  Construction & Assignment //
// -------------------------- //

template <typename Base>
    template <typename T, typename>
Timer<Base>::Timer(T &&t) : Base(std::forward<T>(t))
{
    // Nothing to do here.
}

template <typename Base>
    template <typename T, typename>
auto Timer<Base>::operator=(T &&t)
    -> Timer &
{
    Base::operator=(std::forward<T>(t));
    return *this;
}

// ------------------------ //
//   Arithmetic Operations  //
// ------------------------ //

template <typename Base>
auto Timer<Base>::multiply(const VectorType &v) const
    -> VectorType
{
    auto t1 = steady_clock::now();
    auto w  = Base::multiply(v);
    auto t2 = steady_clock::now();
    multiply_mv_time += duration_cast<duration<double>>(t2 - t1);
    return w;
}

template <typename Base>
auto Timer<Base>::multiply(const MatrixType &B) const
    -> MatrixType
{
    auto t1 = steady_clock::now();
    auto C  = Base::multiply(B);
    auto t2 = steady_clock::now();
    multiply_mm_time += duration_cast<duration<double>>(t2 - t1);
    return C;
}

template <typename Base>
auto Timer<Base>::solve(const VectorType &v) const
    -> VectorType
{
    auto t1 = steady_clock::now();
    auto w  = Base::solve(v);
    auto t2 = steady_clock::now();
    solve_time += duration_cast<duration<double>>(t2 - t1);
    return w;
}

template <typename Base>
auto Timer<Base>::invert() const
    -> MatrixType
{
    auto t1 = steady_clock::now();
    auto B  = Base::invert();
    auto t2 = steady_clock::now();
    invert_time += duration_cast<duration<double>>(t2 - t1);
    return B;
}
