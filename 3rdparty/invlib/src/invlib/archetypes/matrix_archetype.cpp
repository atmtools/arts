#include <iostream>

template <typename Real>
MatrixArchetype<Real>::MatrixArchetype(const MatrixArchetype<Real> &A)
    : m(A.rows()), n(A.cols())
{
    data = std::unique_ptr<Real[]>(new Real[m * n]);
    std::copy(&A.data[0], &A.data[n*m], &data[0]);
}

template <typename Real>
MatrixArchetype<Real>& MatrixArchetype<Real>::operator=(const MatrixArchetype &A)
{
    m = A.rows();
    n = A.cols();

    data = std::unique_ptr<Real[]>(new Real[m * n]);
    std::copy(&A.data[0], &A.data[n*m], &data[0]);
}

template <typename Real>
void MatrixArchetype<Real>::resize(unsigned int i, unsigned int j)
{
    m = i;
    n = j;
    data = std::unique_ptr<Real[]>(new Real[m * n]);
}

template <typename Real>
Real & MatrixArchetype<Real>::operator()(unsigned int i, unsigned int j)
{
    assert((0 <= i) && (i < m));
    assert((0 <= j) && (j < n));

    return data[i * n + j];
}

template <typename Real>
Real MatrixArchetype<Real>::operator()(unsigned int i, unsigned int j) const
{
    assert((0 <= i) && (i < m));
    assert((0 <= j) && (j < n));

    return data[i * n + j];
}

template <typename Real>
unsigned int MatrixArchetype<Real>::cols() const
{
    return n;
}

template <typename Real>
unsigned int MatrixArchetype<Real>::rows() const
{
    return m;
}

template <typename Real>
void MatrixArchetype<Real>::accumulate(const MatrixArchetype &B)
{
    for (unsigned int i = 0; i < m; i++)
    {
	for (unsigned int j = 0; j < n; j++)
	{
	    data[n * i + j] += B(i,j);
	}
    }
}

template <typename Real>
void MatrixArchetype<Real>::subtract(const MatrixArchetype &B)
{
    for (unsigned int i = 0; i < m; i++)
    {
	for (unsigned int j = 0; j < n; j++)
	{
	    data[n * i + j] -= B(i,j);
	}
    }
}

template <typename Real>
auto MatrixArchetype<Real>::multiply(const MatrixArchetype<Real> &B) const
    -> MatrixArchetype
{
    assert(n == B.rows());

    MatrixArchetype<Real> C; C.resize(m, B.cols());

    for (unsigned int h = 0; h < m; h++)
    {
	for (unsigned int i = 0; i < B.cols(); i++)
	{
	    Real sum = 0.0;
	    for (unsigned int j = 0; j < n; j++)
	    {
		sum += (*this)(h, j) * B(j, i);
	    }
	    C(h, i) = sum;
	}
    }
    return C;
}

template <typename Real>
auto MatrixArchetype<Real>::multiply(const VectorType &v) const
    -> VectorType
{
    assert(n == v.rows());

    VectorType w;
    w.resize(m);

    for (unsigned int i = 0; i < m; i++)
    {
	Real sum = 0.0;
	for (unsigned int j = 0; j < n; j++)
	{
	    sum += (*this)(i, j) * v(j);
	}
	w(i) = sum;
    }
    return w;
}

template <typename Real>
void MatrixArchetype<Real>::scale(Real c)
{
    for (unsigned int i = 0; i < m; i++)
    {
	for (unsigned int j = 0; j < n; j++)
	{
	    (*this)(i, j) *= c;
	}
    }
}

template<typename Real>
auto MatrixArchetype<Real>::solve(const VectorType& v) const
    -> VectorType
{
    assert(n == m);

    MatrixType QR = this->QR();
    VectorType w; w.resize(n);
    w = QR.backsubstitution(v);
    return w;
}

template<typename Real>
auto MatrixArchetype<Real>::invert() const
    -> MatrixType
{
    assert(n == m);

    MatrixType QR = this->QR();
    VectorType v; v.resize(n);
    VectorType w; w.resize(n);

    for (unsigned int i = 0; i < n; i++)
    {
        v(i) = 0.0;
    }

    MatrixType B; B.resize(n, n);

    for (unsigned int i = 0; i < n; i++)
    {
        v(i) = 1.0;
        w = QR.backsubstitution(v);

        for (unsigned int j = 0; j < n; j++)
        {
            B(j, i) = w(j);
        }
        v(i) = 0.0;
    }
    return B;
}

template <typename Real>
auto MatrixArchetype<Real>::QR() const
    -> MatrixType
{
    assert(n == m);

    MatrixArchetype<Real> QR; QR.resize(n,n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            QR(i,j) = 0.0;

    for (int i = 0; i < n; i++)
    {
        // Q Matrix.
        for (int j = 0; j < i; j++)
        {
            Real q_sum = 0.0;
            for (int k = 0; k < j; k++)
            {
                q_sum += QR(i, k) * QR(k, j);
            }
            QR(i, j) = (*this)(i, j) - q_sum;
        }

        // Q Matrix.
        for (int j = 0; j < i; j++)
        {
            Real q_sum = 0.0;
            for (int k = 0; k < j; k++)
            {
                q_sum += QR(j, k) * QR(k, i);
            }
            if (QR(j,j) != 0.0)
                QR(j, i) = ((*this)(j, i) - q_sum) / QR(j,j);
        }

        Real diag_sum = 0.0;
        for (int k = 0; k < i; k++)
            diag_sum += QR(i, k) * QR(k, i);
        QR(i,i) = (*this)(i,i) - diag_sum;
    }
    return QR;
}

template <typename Real>
auto MatrixArchetype<Real>::backsubstitution(const VectorType &b) const
    -> VectorType
{
    assert(n == m);

    VectorType c; c.resize(n);
    VectorType d; d.resize(n);

    for (int i = 0; i < n; i++)
    {
        RealType sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += (*this)(i, j) * c(j);
        }
        c(i) = (b(i) - sum) / (*this)(i,i);
    }

    for (int i = n - 1; i >= 0; i--)
    {
        RealType sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += (*this)(i, j) * d(j);
        }
        d(i) = c(i) - sum;
    }
    return d;
}

template<typename Real>
auto MatrixArchetype<Real>::transpose()
    -> MatrixType
{
    MatrixType B{}; B.resize(n, m);

    for (unsigned int i = 0; i < m; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            B(j,i) = (*this)(i,j);
        }
    }
    return B;
}

template <typename Real>
auto MatrixArchetype<Real>::transpose_multiply(const MatrixArchetype<Real> &B) const
    -> MatrixArchetype
{
    assert(m == B.rows());

    MatrixArchetype<Real> C; C.resize(n, B.cols());

    for (unsigned int h = 0; h < n; h++)
    {
	for (unsigned int i = 0; i < B.cols(); i++)
	{
	    Real sum = 0.0;
	    for (unsigned int j = 0; j < m; j++)
	    {
		sum += (*this)(j, h) * B(j, i);
	    }
	    C(h, i) = sum;
	}
    }
    return C;
}

template <typename Real>
auto MatrixArchetype<Real>::transpose_multiply(const VectorType &v) const
    -> VectorType
{
    assert(m == v.rows());

    VectorType w;
    w.resize(n);

    for (unsigned int i = 0; i < n; i++)
    {
	Real sum = 0.0;
	for (unsigned int j = 0; j < m; j++)
	{
	    sum += (*this)(j, i) * v(j);
	}
	w(i) = sum;
    }
    return w;
}

template <typename Real>
std::ostream & operator<<(std::ostream & out, const MatrixArchetype<Real> &A)
{
    for (unsigned int i = 0; i < A.rows(); i++)
    {
        for (unsigned int j = 0; j < A.cols(); j++)
        {
            out << A(i, j) << " ";
        }
        out << std::endl;
    }
    return out;
}
