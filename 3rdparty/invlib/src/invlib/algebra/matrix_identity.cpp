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
auto MatrixIdentity<Matrix>::diagonal() const
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

template
<
typename T1,
typename T2
>
template <typename T3>
auto MatrixProduct<MatrixIdentity<T1>, T2>::multiply(const T3 &t) const
    -> typename T3::ResultType
{
    using T3ResultType = typename T3::ResultType;

    T3ResultType u = remove_reference_wrapper(B).multiply(t);
    u.scale(A.scale());
    return u;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<MatrixIdentity<T1>, T2>::invert() const
    -> MatrixType
{
    MatrixType C = B.invert();
    C.scale(1.0 / A.scale());
    return C;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<MatrixIdentity<T1>, T2>::solve(const VectorType &u) const
    -> VectorType
{
    VectorType v = B.solve(u);
    v.scale(1.0 / A.scale());
    return v;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<MatrixIdentity<T1>, T2>::transpose() const
    -> MatrixType
{
    MatrixType C = B.transpose();
    C.scale(A.scale());
    return C;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<MatrixIdentity<T1>, T2>::diagonal() const
    -> VectorType
{
    VectorType v = B.diagonal();
    v.scale(A.scale());
    return v;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<MatrixIdentity<T1>, T2>::row(size_t i) const
    -> VectorType
{
    VectorType v = B.row(i);
    v.scale(A.scale());
    return v;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<MatrixIdentity<T1>, T2>::col(size_t i) const
    -> VectorType
{
    VectorType v = B.col(i);
    v.scale(A.scale());
    return v;
}

template
<
typename T1,
typename T2
>
MatrixProduct<MatrixIdentity<T1>, T2>::operator ResultType() const
{
    ResultType C = A.multiply(static_cast<const ResultType &>(B));
    return C;
}
