#include <type_traits>

// ---------------------  //
//  Matrix Product Class  //
// ---------------------  //

template
<
typename T1,
typename T2
>
template <typename T3>
auto MatrixProduct<T1, T2>::multiply(const T3 &t) const
    -> typename T3::ResultType
{
    using T3ResultType = typename T3::ResultType;

    T3ResultType u = remove_reference_wrapper(B).multiply(t);
    T3ResultType v = remove_reference_wrapper(A).multiply(u);
    return v;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<T1, T2>::invert() const
    -> MatrixType
{
    MatrixType D = this->operator ResultType();
    MatrixType E = D.invert();
    return E;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<T1, T2>::solve(const VectorType &u) const
    -> VectorType
{
    VectorType v = B.solve(u);
    VectorType w = A.solve(v);
    return w;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<T1, T2>::transpose() const
    -> MatrixType
{
    MatrixType C = A.multiply((MatrixType) B);
    MatrixType D = C.transpose();
    return D;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<T1, T2>::diagonal() const
    -> VectorType
{
    size_t m = A.rows();
    VectorType diag; diag.resize(m);
    for (size_t i = 0; i < m; i++)
    {
        diag(i) = dot(A.row(i), B.col(i));
    }
    return diag;
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<T1, T2>::row(size_t i) const
    -> VectorType
{
    return B.transpose_multiply(A.row(i));
}

template
<
typename T1,
typename T2
>
auto MatrixProduct<T1, T2>::col(size_t i) const
    -> VectorType
{
    return A.multiply(remove_reference_wrapper(B).col(i));
}

template
<
typename T1,
typename T2
>
MatrixProduct<T1, T2>::operator ResultType() const
{
    ResultType C = A.multiply(static_cast<const ResultType &>(B));
    return C;
}
