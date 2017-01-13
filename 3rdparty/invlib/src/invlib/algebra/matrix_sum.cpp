template<typename T1, typename T2>
MatrixSum<T1, T2>::MatrixSum(T1 Op1, T2 Op2)
        : A(Op1), B(Op2)
{
    // Nothing to do here.
}

template<typename T1, typename T2>
    template<typename T3>
auto MatrixSum<T1, T2>::multiply(const T3 &t) const
    -> typename T3::ResultType
{
    using T3ResultType = typename T3::ResultType;
    T3ResultType u = t;
    T3ResultType v = B.multiply(u);
    u.accumulate(A.multiply(u));
    return u;
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::multiply(const VectorType &v) const
    -> VectorType
{
    VectorType w = B.multiply(v);
    w.accumulate(A.multiply(v));
    return w;
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::multiply(const MatrixType &C) const
    -> MatrixType
{
    MatrixType D = A;
    D.accumulate(B);
    return D.multiply(C);
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::solve(const VectorType &v) const
    -> VectorType
{
    MatrixType C = A;
    C.accumulate(B);
    return C.solve(v);
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::invert() const
    -> MatrixType
{
    MatrixType C = A;
    C.accumulate(B);
    return C.invert();
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::diagonal() const
    -> VectorType
{
    VectorType diag = A.diagonal();
    diag.accumulate(B.diagonal());
    return diag;
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::row(size_t i) const
    -> VectorType
{
    VectorType r = A.row(i);
    r.accumulate(B.row(i));
    return r;
}

template<typename T1, typename T2>
auto MatrixSum<T1, T2>::col(size_t i) const
    -> VectorType
{
    VectorType c = A.col(i);
    c.accumulate(B.col(i));
    return c;
}

template<typename T1, typename T2>
    template<typename T3>
auto MatrixSum<T1, T2>::operator*(T3 &&C) const
    -> Product<T3>
{
    return Product<T3>(*this, C);
}

template<typename T1, typename T2>
    template <typename T3>
auto MatrixSum<T1, T2>::operator+(T3 &&C) const
    -> Sum<T3>
{
    return Sum<T3>(*this, C);
}

template<typename T1, typename T2>
    template <typename T3>
auto MatrixSum<T1, T2>::operator-(T3 &&C) const
    -> Difference<T3>
{
    return Difference<T3>(*this, C);
}

template<typename T1, typename T2>
MatrixSum<T1, T2>::operator ResultType() const
{
    ResultType C = A;
    C.accumulate(B);
    return C;
}
