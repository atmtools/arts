// -------------- //
//  Matrix Class  //
// -------------- //

template<typename Base>
Matrix<Base>::Matrix(const Matrix &t)
    : Base(t)
{
    // Nothing to do here.
}

template<typename Base>
Matrix<Base>::Matrix(Matrix &&t)
    : Base(t)
{
    // Nothing to do here.
}

template<typename Base>
auto Matrix<Base>::operator=(const Matrix &t)
    -> Matrix &
{
    Base::operator=(t);
    return (*this);
}

template<typename Base>
auto Matrix<Base>::operator=(Matrix &&t)
    -> Matrix &
{
    Base::operator=(t);
    return (*this);
}

template<typename Base>
    template<typename T, typename>
Matrix<Base>::Matrix(T &&t)
    : Base(std::forward<T>(t))
{
    // Nothing to do here.
}

template<typename Base>
    template<typename T, typename>
auto Matrix<Base>::operator=(T &&t)
    -> Matrix &
{
    Base::operator=(std::forward<T>(t));
    return *this;
}

template<typename Base>
auto Matrix<Base>::begin()
    -> ElementIterator
{
    return ElementIterator(this);
}

template<typename Base>
auto Matrix<Base>::end()
    -> ElementIterator
{
    return ElementIterator(this, this->rows(), this->cols());
}

template <typename Base>
    template <typename T>
void Matrix<Base>::accumulate(const T &B)
{
    this->Base::accumulate(B);
}

template <typename Base>
void Matrix<Base>::accumulate(const MatrixIdentity<Matrix> &B)
{
    for (unsigned int i = 0; i < this->rows(); i++)
    {
	(*this)(i,i) += B.scale();
    }
}

template <typename Base>
void Matrix<Base>::accumulate(const MatrixZero &Z)
{
    // <Insert obscure method to add zeros to elements here>
}

template <typename Base>
    template <typename T1>
void Matrix<Base>::operator+=(T1 &&Z)
{
    this->accumulate(std::forward<T1>(Z));
}

template <typename Base>
    template<typename T>
auto Matrix<Base>::operator+(T &&B) const
    -> Sum<T>
{
    return Sum<T>(*this, B);
}

template <typename Base>
    template <typename T>
auto Matrix<Base>::operator-(T &&B) const -> Difference<T>
{
    return Difference<T>(*this, B);
}

template <typename Base>
    template<typename T>
auto Matrix<Base>::operator*(T &&B) const
    -> Product<T>
{
    return Product<T>(*this, B);
}

// ----------------------- //
//  ElementIterator Class  //
// ----------------------- //

template<typename Base>
Matrix<Base>::ElementIterator::ElementIterator(MatrixType *M_)
    : M(M_), i(0), j(0), k(0), m(M_->rows()), n(M_->cols())
{
    // Nothing to do here.
}

template<typename Base>
Matrix<Base>::ElementIterator::ElementIterator(MatrixType *M_,
                                               unsigned int i_,
                                               unsigned int j_)
    : M(M_), i(i_), j(j_), k(i_ * j_), m(M_->rows()), n(M_->cols())
{
    // Nothing to do here.
}

template<typename Base>
auto Matrix<Base>::ElementIterator::operator*()
    -> RealType&
{
    return M->operator()(i,j);
}

template<typename Base>
typename Matrix<Base>::ElementIterator& Matrix<Base>::ElementIterator::operator++()
{
    k++;
    i = k / n;
    j = k % n;
    return *this;
}

template<typename Base>
auto Matrix<Base>::ElementIterator::operator!=(ElementIterator it)
    -> bool
{
    return (k != it.k);
}
