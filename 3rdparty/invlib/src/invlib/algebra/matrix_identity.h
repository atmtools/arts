/** \file algebra/matrix_identity.h
 *
 * \brief Symbolic identity matrix.
 *
 * Contains the class invlib::Matrixidentity which implements a generic,
 * symbolic identity matrix which can be used to form algebraic expressions.
 *
 */

#ifndef ALGEBRA_MATRIX_IDENTITY_H
#define ALGEBRA_MATRIX_IDENTITY_H

#include <iostream>
#include "invlib/traits.h"
#include "invlib/algebra/matrix.h"

namespace invlib
{

// -------------------- //
// Forward Declarations //
// -------------------- //

template
<
typename T1,
typename T2
>
class MatrixProduct;

template
<
typename T1,
typename T2
>
class MatrixSum;

// --------------------- //
// Class MatrixIdentity  //
// --------------------- //

/** \brief Symbolic matrix identity.
 *
 * A class template representing a (scaled) matrix identity matrix.
 *
 * \tparam RealType The floating point type used for scalars.
 * \tparam MatrixType The underlying MatrixType type that is used.
 *
 */
template
<
typename Matrix
>
class MatrixIdentity
{

public:

    /*! The basic scalar type. */
    using RealType   = typename Matrix::RealType;
    /*! The basic vector type  */
    using VectorType = typename Matrix::VectorType;
    /*! The basic matrix type. */
    using MatrixType = Matrix;
    /*! The basic matrix type. */
    using ResultType = Matrix;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Create true matrix identity (no scaling) */
    MatrixIdentity();

    /*! Create scaled matrix identity */
    MatrixIdentity(RealType c_);

    MatrixIdentity(const MatrixIdentity &) = default;
    MatrixIdentity(MatrixIdentity &&)      = default;

    MatrixIdentity & operator=(const MatrixIdentity &) = default;
    MatrixIdentity & operator=(MatrixIdentity &&)      = default;

    /*! Return scaling factor of matrix */
    RealType scale() const;

    /*! Return diagonal element of the identity matrix. */
    RealType diagonal() const;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    template <typename T1>
    T1 multiply(T1 &A) const;

    template <typename T1>
    T1 multiply(const T1 &A) const;

    template <typename T1>
    T1 solve(T1 &A) const;

    template <typename T1>
    T1 solve(const T1 &A) const;

    void scale(RealType c);

    // --------------------- //
    // Arrithmetic Operators //
    // --------------------- //

    template <typename T1>
    using Product = MatrixProduct<MatrixIdentity, T1>;

    template<typename T1>
    Product<T1> operator*(T1 &&A) const;

private:

    RealType c;

};

/** \brief Multiplication by a scalar.
 *
 * Overload of the * operator for multiplication of an algebraic
 * expression by a scalar.
 *
 * \tparam T1 The type of the algebraic expression.
 *
 * \param c The scaling factor.
 * \param B The algebraic expression to be scaled.
 * \return A matrix product proxy object with a scaled identity matrix and
 * the given algebraic expression as operands.
 */
template
<
typename T1,
typename RealType = typename decay<T1>::RealType,
typename MatrixType = typename decay<T1>::MatrixType,
typename = disable_if< is_same<decay<T1>, MatrixIdentity<MatrixType> > >
>
auto operator*(double c, T1&& B)
    -> MatrixProduct<MatrixIdentity<MatrixType>, T1>

{
    using I = MatrixIdentity<MatrixType>;
    using P = typename I::template Product<T1>;
    return P(I(c), B);
}

/** \brief Scale identity matrix.
 *
 * Overload of the * operator for multiplication of an identity
 * matrix by a scalar.
 *
 * \tparam T1 The underlying matrix type of the identity matrix object.
 *
 * \param c The scaling factor.
 * \param B The identity matrix object to be scaled.
 * \return A scaled identity matrix object.
 *
 */
template
<
typename T1
>
auto operator*(double c,
               const MatrixIdentity<T1>& A)
    -> MatrixIdentity<T1>
{
    MatrixIdentity<T1> B(A);
    B.scale(c);
    return B;
}

/** \brief Identity matrix inverse.
 *
 * Compute the inverse of a (scaled) identity matrix.
 *
 * \tparam T1 The underlying matrix type of the identity matrix object.
 * \param B The identity matrix to be inverted.
 * \return The inverted identity matrix object.
 *
 */
template
<
typename T1
>
auto inv(const MatrixIdentity<T1> &A)
    -> MatrixIdentity<T1>
{
    return MatrixIdentity<T1>(1.0 / A.scale());
}

// ---------------------------------------------- //
// Partial Specialization of Matrix Product Class //
// ---------------------------------------------- //

template
<
typename T1,
typename T2
>
class MatrixProduct<MatrixIdentity<T1>, T2>
{

public:

    /*! The basic scalar type. */
    using RealType   = typename decay<T2>::RealType;
    /*! The basic vector type  */
    using VectorType = typename decay<T2>::VectorType;
    /*! The basic matrix type. */
    using MatrixType = typename decay<T2>::MatrixType;
    /*! The type of the result of the expression */
    using ResultType = typename decay<T2>::ResultType;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    MatrixProduct(const MatrixIdentity<T1> & A_, T2 B_ )
        : A(A_), B(B_) {}
    MatrixProduct(const MatrixProduct &) = default;
    MatrixProduct(MatrixProduct &&)      = default;
    MatrixProduct & operator=(const MatrixProduct &) = default;
    MatrixProduct & operator=(MatrixProduct &&)      = delete;
    ~MatrixProduct() = default;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    template <typename T3>
    auto multiply(const T3 &u) const -> typename T3::ResultType;

    MatrixType invert() const;
    VectorType solve(const VectorType &v) const;
    MatrixType transpose()   const;
    VectorType diagonal()    const;
    VectorType row(size_t i) const;
    VectorType col(size_t i) const;

    // --------------------- //
    //  Algebraic Operators  //
    // --------------------- //

    template <typename T3>
    using NestedProduct = typename decay<T2>::template Product<T3>;
    template <typename T3>
    using Product = MatrixProduct<MatrixIdentity<T1>, NestedProduct<T3>>;

    template <typename T3>
    Product<T3> operator*(T3 &&C) const
    {
	return Product<T3>(A, B * std::forward<T3>(C));
    }

    template <typename T3>
    using Sum = MatrixSum<MatrixProduct , T3>;
    template <typename T3>
    Sum<T3> operator+(T3 &&C) const
    {
	return Sum<T3>(*this, C);
    }

    template <typename T3>
    using Difference = MatrixDifference<MatrixProduct , T3>;
    template <typename T3>
    Difference<T3> operator-(T3 &&C)
    {
	return Difference<T3>(*this, C);
    }

    operator ResultType() const;

private:

    MatrixIdentity<T1> A;
    T2 B;
};

#include "matrix_identity.cpp"

}      // namespace invlib

#endif // ALGEBRA_MATRIX_IDENTITY_H
