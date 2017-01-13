template <typename Real>
VectorArchetype<Real>::VectorArchetype(const VectorArchetype<Real> &v)
    : n(v.rows())
{
    data = std::unique_ptr<Real[]>(new Real[n]);
    std::copy(&v.data[0], &v.data[n], &data[0]);
}

template <typename Real>
auto VectorArchetype<Real>::operator=(const VectorArchetype<Real> &v)
    -> VectorArchetype &
{
    n = v.rows();
    data = std::unique_ptr<Real[]>(new Real[n]);
    std::copy(&v.data[0], &v.data[n], &data[0]);
    return *this;
}

template <typename Real>
auto VectorArchetype<Real>::get_block(unsigned int i,
                                      unsigned int di) const
    -> VectorArchetype
{
    assert((di > 0) &&(i + di <= n));

    VectorArchetype block; block.resize(di);
    int data_start = i;
    int data_end = data_start + di;

    std::copy(&data[data_start], &data[data_end], &block.data[0]);

    return block;
}

template <typename Real>
void VectorArchetype<Real>::resize(unsigned int i)
{
    n = i;
    data = std::unique_ptr<Real[]>(new Real[n]);
}

template <typename Real>
Real & VectorArchetype<Real>::operator()(unsigned int i)
{
    return data[i];
}

template <typename Real>
Real VectorArchetype<Real>::operator()(unsigned int i) const
{
    return data[i];
}


template <typename Real>
unsigned int VectorArchetype<Real>::rows() const
{
    return n;
}

template <typename Real>
Real* VectorArchetype<Real>::data_pointer(int i)
{
    return data.get() + i;
}

template <typename Real>
const Real* VectorArchetype<Real>::data_pointer(int i) const
{
    return data.get() + i;
}

template <typename Real>
void VectorArchetype<Real>::accumulate(RealType c)
{
    for (unsigned int i = 0; i < n; i++)
    {
	(*this)(i) += c;
    }
}

template <typename Real>
void VectorArchetype<Real>::accumulate(const VectorArchetype<Real> &v)
{
    assert(n == v.rows());

    for (unsigned int i = 0; i < n; i++)
    {
	(*this)(i) += v(i);
    }
}

template <typename Real>
void VectorArchetype<Real>::subtract(const VectorArchetype<Real>& v)
{
    assert(n == v.rows());

    for (unsigned int i = 0; i < n; i++)
    {
	(*this)(i) -= v(i);
    }
}

template <typename Real>
void VectorArchetype<Real>::scale(Real c)
{
    for (unsigned int i = 0; i < n; i++)
    {
	(*this)(i) *= c;
    }
}

template <typename Real>
Real VectorArchetype<Real>::norm() const
{
    return sqrt(dot(*this, *this));
}

template <typename Real>
auto VectorArchetype<Real>::element_multiply(const VectorArchetype & v) const
    -> VectorArchetype
{
    VectorArchetype w{}; w.resize(rows());
    for (size_t i = 0; i < rows(); i++)
    {
        w(i) = operator()(i) * v(i);
    }
    return w;
}

template <typename Real>
void VectorArchetype<Real>::element_invert()
{
    for (size_t i = 0; i < rows(); i++)
    {
        operator()(i) = 1.0 / operator()(i);
    }
}

template <typename Real>
Real dot(const VectorArchetype<Real> &v, const VectorArchetype<Real> &w)
{
    Real sum = 0.0;
    for (unsigned int i = 0; i < v.rows(); i++)
    {
	sum += v(i) * w(i);
    }
    return sum;
}

template <typename Real>
std::ostream & operator<<(std::ostream &out, const VectorArchetype<Real>& v)
{
    out << "[";
    for (unsigned int i = 0; i < v.rows()-1; i++)
    {
	out << v(i) << ", ";
    }
    out << v(v.rows() - 1) << "] " << std::endl;
}
