template <typename Matrix>
MatrixIdentity<Matrix>::MatrixIdentity() : c(1.0)
{
    // Nothing to do here.
}

template <typename Matrix>
MatrixIdentity<Matrix>::MatrixIdentity(RealType c_) : c(c_)
{
    // Nothing to do here.
}

template <typename Matrix>
auto MatrixIdentity<Matrix>::scale() const
    ->  RealType
{
    return c;
}

template <typename Matrix>
    template <typename T1>
auto MatrixIdentity<Matrix>::multiply(T1 &A) const
    -> T1
{
    A.scale(c);
    return A;
}

template <typename Matrix>
    template <typename T1>
auto MatrixIdentity<Matrix>::multiply(const T1 &A) const
    -> T1
{
    T1 B(A);
    B.scale(c);
    return B;
}

template <typename Matrix>
    template <typename T1>
auto MatrixIdentity<Matrix>::solve(const T1 &A) const
    -> T1
{
    T1 B(A);
    B.scale(1.0 / c);
    return B;
}

template <typename Matrix>
    template <typename T1>
auto MatrixIdentity<Matrix>::solve(T1 &A) const
    -> T1
{
    A.scale(1.0 / c);
    return A;
}

template <typename Matrix>
void MatrixIdentity<Matrix>::scale(RealType d)
{
    c *= d;
}

template <typename Matrix>
    template<typename T1>
auto MatrixIdentity<Matrix>::operator*(T1 &&A) const
    -> Product<T1>
{
    return Product<T1>(*this, A);
}
