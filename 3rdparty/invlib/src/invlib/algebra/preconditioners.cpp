template <typename VectorType>
template <typename T>
JacobianPreconditioner<VectorType>::JacobianPreconditioner(const T & M)
{
    diag = M.diagonal();
    diag.element_invert();
}

template <typename VectorType>
auto JacobianPreconditioner<VectorType>::operator()(const VectorType &v) const
    -> VectorType
{
    VectorType w; w.resize(v.rows());
    w = diag.element_multiply(v);
    return w;
}
