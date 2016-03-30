template<typename T1, typename T2>
MatrixDifference<T1, T2>::MatrixDifference(T1 A_, T2 B_)
    : A(A_), B(B_)
{
    // Nothing to do here.
}

template<typename T1, typename T2>
auto MatrixDifference<T1, T2>::multiply(const VectorType &u) const
    -> VectorType
{
    VectorType v = A.multiply(u);
    VectorType w = B.multiply(u);
    v.subtract(w);
    return v;
}

template<typename T1, typename T2>
auto MatrixDifference<T1, T2>::multiply(const MatrixType &C) const
    -> MatrixType
{
    MatrixType D = A;
    D.subtract(B);
    return D.multiply(C);
}

template<typename T1, typename T2>
auto MatrixDifference<T1, T2>::solve(const VectorType &v) const
    -> VectorType
{
    MatrixType C = A;
    C.subtract(B);
    return C.solve(v);
}

template<typename T1, typename T2>
auto MatrixDifference<T1, T2>::invert() const
    -> MatrixType
{
    MatrixType C = A;
    C.subtract(B);
    return C.invert();
}

template<typename T1, typename T2>
    template<typename T3>
auto MatrixDifference<T1, T2>::operator*(T3 &&C) const
    -> Product<T3>
{
    return Product<T3>(*this, C);
}

template<typename T1, typename T2>
    template<typename T3>
auto MatrixDifference<T1, T2>::operator+(T3 &&C) const
    -> Sum<T3>
{
    return Sum<T3>(*this, C);
}

template<typename T1, typename T2>
    template<typename T3>
auto MatrixDifference<T1, T2>::operator-(const T3 &&C) const
    -> Difference<T3>
{
    return Difference<T3>(*this, C);
}

template<typename T1, typename T2>
MatrixDifference<T1, T2>::operator ResultType() const
{
    ResultType C = A;
    C.subtract(B);
    return C;
}
