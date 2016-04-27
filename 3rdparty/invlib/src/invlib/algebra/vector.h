/** \file algebra/vector.h
 *
 * \brief Contains Vector class template for symbolic computations on generic
 * vector types.
 *
 */

#ifndef ALGEBRA_VECTOR_H
#define ALGEBRA_VECTOR_H

#include "invlib/algebra/matrix.h"
#include "invlib/algebra/matrix_sum.h"
#include "invlib/traits.h"
#include <utility>
#include <iostream>

namespace invlib
{

/**
 * \brief Wrapper class for symbolic computations involving vectors.
 *
 * The Vector class provides an abstraction layer to delay computations
 * involving vectors. The following delayed operations are provided:
 *
 * - Addition: operator+()
 * - Subtraction: operator-()
 *
 * Those operations return a proxy type representing the computation, that
 * can be combined with the other matrix operations. The computation is delayed
 * until the resulting proxy object is converted to a matrix or a vector.
 *
 * \tparam Base The undelying vector type to be used.
 *
 */
template
<
typename Base
>
class Vector : public Base
{
public:

    // ------------------- //
    //  Element Iterator   //
    // ------------------- //

    class ElementIterator;

    /*!
     *\return An element iterator object pointing to the first element
     *in the vector.
     */
    ElementIterator begin();

    /*!
     *\return An element iterator pointing to the end of the vector.
     */
    ElementIterator end();

    struct INVALID_TYPE
    {
        using VectorType = Vector;
    };

    // -------------- //
    //  Type Aliases  //
    // -------------- //

    /* template <typename T1> */
    /* using Product = LEFT1_VECTOR_MULTIPLY_NOT1_SUPPORTED; */

    using BaseType   = Base;
    /*! The basic scalar type. */
    using RealType   = typename Base::RealType;
    /*! The basic vector type  */
    using VectorType = Vector;
    /*! The basic matrix type. */
    using MatrixType = Matrix<typename Base::MatrixType>;
    /*! The type of the result of the expression */
    using ResultType = Vector;
    /*! The type of the result of the expression */
    template<typename T1>
    using Product = INVALID_TYPE;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Default constructor. */
    Vector() = default;

    /*! Default copy constructor.
     *
     * Call the Base copy constructor, which should perform a deep copy of
     * the provided vector v.
     *
     * /param v The vector v to be copied from.
     */
    Vector(const Vector& v) = default;

    /*! Default move constructor.
     *
     * Call the Base move constructor, which should move the data references
     * from the vector v into this Vector object.
     *
     * /param v The vector v to be moved from.
     */
    Vector(Vector&& v)      = default;

    /*! Default assignment operator.
     *
     * Call the Base assignment operator, which should copy the values of
     * the provided vector v into this object.
     *
     * \param v The vector v to be assigned from.
     */
    Vector& operator=(const Vector& v) = default;

    /*! Default move assignment operator.
     *
     * Call the Base move assignment operator, which should copy the values of
     * the provided vector v into this object.
     *
     * \param v The vector v to be moved from.
     */
    Vector& operator=(Vector&& v) = default;

    /*! Copy base object.
     *
     * Copies the provided base object into the vector.
     *
     * \param v The vector to be copied.
     */
    //Vector(const Base& v);

    /*! Move base object into vector.
     *
     * Moves the provided base object into the vector. In this way results
     * from arithmetic operations on the base type can be moved directly
     * into a Vector object in order to avoid expensive copying.
     *
     * \param v The temporary vector of base type to be moved from.
     */
    //Vector(Base&& v);
    template <typename T,
    typename = disable_if<is_same<decay<T>, Vector>>>
    Vector(T &&t) : Base(std::forward<T>(t)) {}

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    /*! Just a convenience wrapper for the accumulate member function
     * of the base type.
     */
    template <typename T1>
    void operator+=(T1 &&);

    /*! Proxy type for the sum of two vectors. */
    template <typename T1>
        using Sum = MatrixSum<const Vector &, T1>;

    /*! Create sum arithmetic expression.
     *
     * \tparam T1 The type of the object to add to this vector.
     * \return An algebraic expression object representing the sum of
     * this vector and the provided argument.
     */
    template<typename T1>
    Sum<T1> operator+(T1 &&v) const;

    /*! Proxy type for the difference of two vectors. */
    template <typename T1>
    using Difference = MatrixDifference<const Vector &, T1>;

    /*! Create difference algebraic expression.
     *
     * \tparam The type of object to subtract from the vector.
     * \return An algebraic expression object representing the difference
     * of this vector and the given argument.
     */
    template <typename T1>
    Difference<T1> operator-(T1 &&v) const;

};

template <typename T>
class Id
{
public:
    using type = T;
};

/*! Dot product of two vectors.
 *
 * The dot product is computed by callling the dot(T1 t) member function
 * of the first argument vector with the second argument vector are argument.
 *
 * \tparam Base The fundamental vector type.
 * \return The dot product of the two vectors.
 */
template
<
typename T1,
typename T2,
typename VectorType = typename T1::VectorType
>
auto dot(const T1 &v,const T2 &w)
    -> typename VectorType::RealType;

/**
 * \brief Iterator for element-wise acces.
 *
 * Iterates over the elements in the vector using
 * <tt>operator()(unsigned int)</tt> for element access. Assumes
 * the length of the vector can be obtained usin <tt>cols()</tt> and
 * that indexing starts at 0.
 *
 * \tparam Base The undelying vector type to be used.
 *
 */
template
<
typename Base
>
class Vector<Base>::ElementIterator
{
public:

    ElementIterator() = default;
    ElementIterator(VectorType* v_);
    ElementIterator(VectorType* v_, unsigned int k);

    RealType& operator*();
    RealType& operator++();
    bool operator!=(ElementIterator it);

private:

    VectorType *v;
    unsigned int k, n;
};

#include "vector.cpp"

}      // namespace invlib

#endif // ALGEBRA_VECTOR_H
