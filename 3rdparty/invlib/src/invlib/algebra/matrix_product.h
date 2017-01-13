/** \file algebra/matrix_product.h
 *
 * \brief Proxy class for computing the product of two algebraic expressions.
 *
 */

#ifndef ALGEBRA_MATRIX_PRODUCT_H
#define ALGEBRA_MATRIX_PRODUCT_H

#include <iostream>
#include <utility>
#include <invlib/traits.h>

namespace invlib
{

// -------------------- //
// Forward Declarations //
// -------------------- //

template
<
typename Base
>
class Vector;

template
<
typename T1,
typename T2
>
class MatrixSum;

template
<
typename T1,
typename T2
>
class MatrixDifference;

// ---------------------  //
//  Matrix Product Class  //
// ---------------------  //

/** \brief Proxy class for computing the product of two algebraic expressions.
 *
 * The MatrixProduct class template provides a proxy class template for the
 * product of two algebraic expressions. The class is basically a placeholder
 * that holds references to the matrices involved in the product. It is also
 * possible to build nested algebraic expression in which the operands may be
 * other proxy classes. In this case, copies of the operands are held by
 * the MatrixProduct object.
 *
 * The result of the product must be either of matrix type or vector type,
 * depending of the result type of the right hand operator. The computation
 * of the product is invoked, when the product is converted to its result type.
 *
 *
 * The class expects the right hand
 * operator to be convertible to either a vector or a matrix and the left
 * hand operator to provide the two member functions
 *
 * - VectorType multiply(const VectorType &) const
 * - MatrixType multiply(const MatrixType &) const
 *
 * \tparam T1 The type of the left hand side operator
 * \tparam T2 the type of the right hand side operator
 * \tparam Matrix The underlying matrix type used.
 *
 */
template
<
typename T1,
typename T2
>
class MatrixProduct
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

    /*! \name Constructors and Destructors
     */

    //@{

    /*!
     * Create a proxy object for the matrix product \f$C = A  B\f$  where
     * \f$A\f$ is the left hand operator and  \f$B\f$ is the right hand
     * operator.
     */
    MatrixProduct(T1 A_, T2 B_ )
        : A(A_), B(B_) {}

    /*! Default copy constructor.
     *
     * Proxy objects are lightweight and can be efficiently copied. However,
     * note that there is no caching in place to reuse the results from identical
     * arithmetic operations.
     *
     * \todo Enable caching.
     *
     */
    MatrixProduct(const MatrixProduct &) = default;

    /*! Default move constuctor. */
    MatrixProduct(MatrixProduct &&)      = default;

    /*! Default assignment operator. */
    MatrixProduct & operator=(const MatrixProduct &) = default;
    MatrixProduct & operator=(MatrixProduct &&)      = delete;

    /*! Default destructor. */
    ~MatrixProduct() = default;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    ///@}
    /*! \name Nested Evaluation
     *
     * These functions are used for the evaluation of the product when
     * it is nested in another algebraic expression.
     *
     */
    ///@{

    /*! Multiply product from the right.
     *
     * If a MatrixProduct \f$C = A B\f$ object is multiplied by a matrix or
     * vector expression, the matrix or vector expression is evaluated first
     * to its result type and then this intermediate result is successively
     * multiplied with B and A.
     *
     * The result type of the the product of a matrix product and a matrix or
     * vector expression is the result type of the matrix or vector expression.
     *
     * Each of the matrices A and B must thus provide multiply methods for
     * multiplying them from the right with the result type.
     *
     * \tparam T3 The type of the matrix or vector expression being multiplied
     * from the right.
     * \return The product of the matrix product \f$A B\f$ and the matrix or
     * vector expression \f$t\f$
     */
    template <typename T3>
    auto multiply(const T3 &u) const -> typename T3::ResultType;

    /*! Inverse of a matrix product.
     *
     * The matrix inverse of a matrix product \f$C = AB\f$ is computed by first
     * evaluating the product and then inverting it. The inversion requires the
     * full evaluate of the left hand operand. The product is then computed by
     * calling the multiply(const MatrixType &) member function of the left hand
     * operand. The inverse is computed by calling the invert() member function
     * of the evaluated product of the two operands.
     *
     *  \return The inverse of the matrix product.
     */
    MatrixType invert() const;

    /*! Solve linear system.
     *
     * Solve the linear system \f$AB w = v\f$ corresponding to the MatrixProduct
     * object. Instead of assembling the complete system, the system is solved
     * by first solving the system corresponding to \f$A\f$ and then the one
     * corresponding to \f$B\f$. This is computationally advantageous to
     * first assembling the system \f$C = AB\f$ and the solving the full system.
     *
     * The solution of the system is delegated to the operand types by calling
     * the solve(const VectorType &) member function.
     *
     * \return The solution \f$ w = (AB)^{-1}v \f$ of the linear system
     * given by \f$C = AB\f$
     */
    VectorType solve(const VectorType &v) const;

    /*! Transpose of the matrix product.
     *
     * Computes the transpose of the matrix product. To this end first
     * the product is evaluated to a MatrixType object and then the
     * member function transpose() of this result is called to compute
     * the transpose.
     *
     * Except for the result matrix, computation of the intermediate matrix
     * requires the allocation of an additional matrix.
     *
     * \return The transpose \f$ (AB)^T = B^T A^T \f$ of the matrix product
     * \f$ AB \f$.
     */
    MatrixType transpose() const;

    /*! Return the vector representing the diagonal of the matrix.
     *
     * Creates the diagonal vector computing the dot products of the rows in the
     * left operand with the corresponding columns of the right operand.
     *
     * \return The vector representing the diagonal of the product.
     */
    VectorType diagonal() const;

    /*! Return row \p i of the matrix product.
     *
     * Computes the ith row of the matrix product by multiplying the transpose of the
     * right-hand operand with the ith row of the left-hand operand.
     *
     * \arg i  The index of the row to compute.
     * \return The ith row of the matrix product.
     */
    VectorType row(size_t i) const;

    /*! Return column \p i of the matrix product.
     *
     * Computes the ith column of the matrix product by the
     * left-hand operand with the ith column of the right-hand operand.
     *
     * \arg i  The index of the column to compute.
     * \return The ith row of the matrix product.
     */
    VectorType col(size_t i) const;

    ///@}
    /*! \name Algebraic Operators
     *
     * The algebraic operator defined on the MatrixProduct class allow to
     * construct nested expressiong from MatrixProducts.
     */
    ///@{

    /*!
     * Type of the product of the right hand operand with a given object
     *  of type T3.
     */
    template <typename T3>
    using NestedProduct = typename decay<T2>::template Product<T3>;

    /*!
     * Type of the product of the the matrix product with another given object
     *  of type T3.
     */
    template <typename T3>
    using Product = MatrixProduct<T1, NestedProduct<T3>>;

    /*! Create a nested matrix product.
     *
     * Creates a proxy object for the product of the given matrix product
     * and another algebraic object C of type T3.
     *
     * \return The proxy object representing the multiplication.
     */
    template <typename T3>
    Product<T3> operator*(T3 &&C) const
    {
	return Product<T3>(A, B * std::forward<T3>(C));
    }

    /*!
     * Type of the sum of a MatrixProduct object and another given object
     * of type T3.
     */
    template <typename T3>
    using Sum = MatrixSum<MatrixProduct , T3>;

    /*! Create a nested matrix sum.
     *
     * Creates a proxy object for the sum of the given matrix product
     * and another algebraic object C of type T3.
     *
     * \return The proxy object representing the summation.
     */
    template <typename T3>
    Sum<T3> operator+(T3 &&C) const
    {
	return Sum<T3>(*this, C);
    }

    /*!
     * Type of the difference of a MatrixProduct object and another given object
     * of type T3.
     */
    template <typename T3>
    using Difference = MatrixDifference<MatrixProduct , T3>;

    /*! Create a nested matrix difference.
     *
     * Creates a proxy object for the difference of the given matrix product
     * and another algebraic object C of type T3.
     *
     * \return The proxy object representing the difference.
     */
    template <typename T3>
    Difference<T3> operator-(T3 &&C)
    {
	return Difference<T3>(*this, C);
    }

    /*! Evaluate the product.
     *
     * Evaluates the given product. To this end the right hand operator
     * is converted to ResultType and the result is then computed by
     * calling the multiply(const ResultType &) member function of the left
     * hand argument.
     *
     * \return The result of the product.
     *
     */
    operator ResultType() const;

    //@}

private:

    T1 A;
    T2 B;

};

#include "matrix_product.cpp"

}      // namespace invlib

#endif // ALGEBRA_MATRIX_PRODUCT_H
