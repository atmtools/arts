/** \file algebra/matrix.h
 *
 * \brief Contains Matrix class template for symbolic computations on generic
 * matrix types.
 *
 */

#ifndef ALGEBRA_MATRIX_H
#define ALGEBRA_MATRIX_H

#include <utility>

#include "invlib/algebra/matrix_product.h"
#include "invlib/algebra/matrix_sum.h"
#include "invlib/algebra/matrix_difference.h"
#include "invlib/algebra/matrix_identity.h"
#include "invlib/algebra/matrix_zero.h"
#include "invlib/traits.h"

namespace invlib
{

template <typename Base>
class Vector;


// -------------- //
//  Matrix Class  //
// -------------- //

/**
 * \brief Wrapper class for symbolic matrix computations.
 *
 * The Matrix class template provides an abstraction layer for delayin matrix
 * computations on a given base matrix type Base and a corresponding vector type
 * VectorType. The following delayed operations are provided:
 *
 * - Addition: operator+()
 * - Subtraction: operator-()
 * - Multiplication: operator*()
 * - Computing the inverse matrix: inv()
 * - Solving the system for a given vector: solve()
 * - Transpose: transp(const Matrix&)
 *
 * All delayed operations return proxy classes representing the computation,
 * which contain only reference to the matrices and vectors holding the actual
 * data. The execution of the computation is triggered when the proxy class
 * is converted to the given base matrix type or vector. During the computation
 * the operations are forwarded to the base matrix type.
 *
 * \tparam Base The base matrix type.
 * \tparam VectorType The corresponding vector type.
 */
template
<
typename Base
>
class Matrix : public Base
{

public:

    // ------------------- //
    //  Element Iterator   //
    // ------------------- //

    class ElementIterator;

    /*!
     *\return An element iterator object pointing to the first element
     *in the matrix.
     */
    ElementIterator begin();

    /*!
     *\return An element iterator pointing to the end of the matrix.
     */
    ElementIterator end();

    // -------------- //
    //  Type Aliases  //
    // -------------- //

    /*! The basic scalar type. */
    using RealType   = typename Base::RealType;
    /*! The basic vector type  */
    using VectorType = Vector<typename Base::VectorType>;
    /*! The basic matrix type. */
    using MatrixType = Matrix;
    /*!
     * Result type of an algebraic expression with Matrix as right hand
     * operator
     */
    using ResultType = Matrix<typename Base::ResultType>;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    Matrix() = default;

    Matrix(const Matrix &A);
    Matrix(Matrix &&A);

    Matrix & operator=(const Matrix &A);
    Matrix & operator=(Matrix &&A);

    /*! Perfect forwarding copy constructor. If the provided argument can
     * is a base class of Matrix or a Base object can be constructed from it,
     * the call is forwarded to the copy constructor of Base using perfect
     * forwarding.
     *
     * \param t The object to be copied from.
     */
    template
    <
    typename T,
    typename = enable_if_either<is_base<decay<T>, Matrix>,
                                is_constructible<Base, T>>
    >
    Matrix(T &&t);

    /*! Perfect forwarding assignment operator. If the provided argument can
     * is a base class of Matrix or a Base object can be constructed from it,
     * the call is forwarded to the copy constructor of Base using perfect
     * forwarding.
     *
     * \param t The object to be copied from.
     */
    template
    <
    typename T,
    typename = enable_if_either<is_base<decay<T>, Matrix>,
                                is_assignable<Base, T>>
    >
    Matrix & operator=(T &&);

    template <typename T>
    void accumulate(const T& B);

    /*! Accumulate identity matrix into this matrix.
     *
     * Overload for the multiply(...) method provided by the base matrix
     * type to handle symbolic identity matrices.
     *
     * \tparam The scalar type of the given IdenityMatrix object.
     * \param  The identity matrix to accumulate into this matrix.
     */
    void accumulate(const MatrixIdentity<Matrix> &B);

    /*! Accumulate zero matrix into this matrix.
     *
     * \param The MatrixZero object to accumulate in to this matrix.
     */
    void accumulate(const MatrixZero &Z);

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    /*! Just a convenience wrapper for the accumulate member function
     * of the base type.
     */
    template <typename T1>
    void operator+=(T1 &&Z);

    /*!
     * Proxy type template for the sum of a matrix and another algebraic expression
     * of given type.
     */
    template <typename T>
        using Sum = MatrixSum<const Matrix&, T>;

    /*! Create sum arithmetic expression.
     *
     * \tparam T1 The type of the object to add to this matrix.
     * \return An algebraic expression object representing the sum of
     * this matrix and the provided argument.
     */
    template<typename T>
	Sum<T> operator+(T &&B) const;

    /*
     * Proxy template type for the difference of a matrix and another 
     * algebraic expression of given type.
     */
    template <typename T>
    using Difference = MatrixDifference<const Matrix&, T>;

    /*! Create difference arithmetic expression.
     *
     * \tparam T1 The type of the object to subtract from this matrix.
     * \return An algebraic expression object representing the difference of
     * this matrix and the provided argument.
     */
    template <typename T>
    Difference<T> operator-(T &&C) const;

    /*
     * Proxy template type for the product of a matrix and another 
     * algebraic expression of given type.
     */
    template <typename T>
    using Product = MatrixProduct<const Matrix&, T>;

    /*! Create product arithmetic expression.
     *
     * \tparam T1 The type of the object to multiply this matrix with.
     * \return An algebraic expression object representing the product of
     * this matrix and the provided argument.
     */
    template<typename T>
    Product<T> operator*(T &&B) const;

};

// ---------------------- //
//  ElementIterator Class //
// ---------------------- //

/**
 * \brief Iterator for matrix element access.
 *
 * Provides an interface for element-wise access of matrix elements. Accesses
 * elements using <tt>Base::operator()(unsigned int, unsigned int)</tt>. Assumes
 * that the number of rows and columns in the matrix can be obtained from the
 * member functions <tt>rows()</tt> and <tt>cols()</tt>, respectively. Assumes
 * indexing to starting at 0.
 *
 * \tparam Base The base matrix type.
 * \tparam VectorType The corresponding vector type.
 */
template
<
typename Base
>
class Matrix<Base>::ElementIterator
{
public:

    ElementIterator() = default;
    ElementIterator(MatrixType* M_);
    ElementIterator(MatrixType* M_, unsigned int i, unsigned int j);

    RealType& operator*();
    RealType& operator++();

    bool operator!=(ElementIterator it);

private:

    MatrixType *M;
    unsigned int i, j, k, n, m;

};

#include "matrix.cpp"

}      // namespace invlib

#endif // ALGEBRA_MATRIX_H

