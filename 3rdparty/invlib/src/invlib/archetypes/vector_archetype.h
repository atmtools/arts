/**
 * \file archtypes/vector_archetype.h
 *
 * \brief Contains VectorArchetype class, which is an archetype for the basic
 * vector type used.
 *
 */
#ifndef ARCHETYPE_VECTOR_ARCHETYPE
#define ARCHETYPE_VECTOR_ARCHETYPE

#include <cassert>
#include <cmath>
#include <iterator>
#include <memory>

namespace invlib
{

// -------------------- //
// Forward Declarations //
// -------------------- //

template
<
typename RealType
>
class MatrixArchetype;

template
<
typename RealType
>
class VectorArchetype;

template
<
typename RealType
>
RealType dot(const VectorArchetype<RealType>&, const VectorArchetype<RealType>&);

// ------------------------  //
//   Vector Archetype Class  //
// ------------------------  //

/*! Vector archetype.
 *
 * Implements a straight forward dense vector class to verify the
 * generic matrix algebra and illustrate the required interface
 * of the fundamental vector type.
 *
 * \tparam The floating point type used to represent scalars.
 */
template
<
typename Real
>
class VectorArchetype
{
public:

    /*! The floating point type used to represent scalars. */
    using RealType   = Real;
    /*! The fundamental vector type used for the matrix algebra.*/
    using VectorType = VectorArchetype<RealType>;
    /*! The fundamental matrix type used for the matrix algebra.*/
    using MatrixType = MatrixArchetype<RealType>;
    /*! The result type of multiplying an algebraic expression with this
     * matrix from the right.
     */
    using ResultType = VectorArchetype<RealType>;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Default constructor. */
    VectorArchetype() = default;

    /*! Copy constructor
     *
     * Performs a deep copy of the vector elements so that the constructed
     * vector is identical but independent of the vector given to the
     * constructor
     *
     */
    VectorArchetype(const VectorArchetype &);

    /*! Move constructor
     *
     * Should move the data references to the data holding the vector
     * elements into this object.
     */
    VectorArchetype(VectorArchetype &&) = default;

    /*! Assignment operator.
     *
     * Performs a deep copy of the vector elements so that the constructed
     * vector is identical but independent of the vector given to the
     * constructor
     *
     */
    VectorArchetype& operator=(const VectorArchetype &);

    /*! Move assignment operator.
     *
     * Move a given vector into this vector by deep copy of its element.
     *
     */
    VectorArchetype& operator=(VectorArchetype &&) = default;
    ~VectorArchetype() = default;

    VectorArchetype get_block(unsigned int i,
                              unsigned int di) const;

    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    /*! Resize vector.
     *
     * Resize the vector to an \f$i\f$ dimensional vector.
     *
     * \param i Number of rows of the resized matrix.
     */
    void resize(unsigned int i);

    /*! Element access.
     *
     * \param i Index of the element to access.
     */
    RealType & operator()(unsigned int i);

    /*! Read-only element access.
     *
     * \param i Index of the element to access.
     */
    RealType operator()(unsigned int i) const;

    /*! Number of rows of the vector
     *
     * \return The number of rows (dimension) of the vector.
     */
    unsigned int rows() const;

    RealType * data_pointer(int i = 0);
    const RealType * data_pointer(int i = 0) const;

    // ------------ //
    //  Arithmetic  //
    // ------------ //

    /*! Accumulate into vector.
     *
     * Element-wise add another vector to this vector.
     *
     * \param v The vector to accumate into this one.
     */
    void accumulate(const VectorArchetype &v);

    /*! Accumulate into vector.
     *
     * Add scalar to all elements in vector. This function is required
     * if the diagonal of a sum involving an identity matrix is to be
     * computed, which is the case if a Jacobian preconditioner is used
     * for the Levenberg-Marquardt method with an identity damping matrix.
     *
     * \param v The scalar to add to each element.
     */
    void accumulate(RealType c);

    /*! Subtract from vector.
     *
     * Element-wise subtract another vector from this one.
     *
     * \param v The vector to subtract from this one.
     */
    void subtract(const VectorArchetype &v);

    /*! Scale vector.
     *
     * Multiply each element by a scalar factor.
     *
     * \param c The factor c to scale the vector with.
     */
    void scale(RealType c);

    /*! Elementwise product of this vector and another vector.
     *
     * \param v The vector to elementwise multiply this vector with.
     */
    VectorArchetype element_multiply(const VectorArchetype & v) const;

    /*! Elementwise inverse of the vector.
     *
     * Sets all elements in the vector to their reciprocal.
     */
    void element_invert();

    /*! Dot product of two vectors
     *
     * \return The dot product \f$ \sum_{i=1}^n v_i w_i \f$ of the
     * two given vectors.
     */
    friend RealType dot<RealType>(const VectorArchetype&, const VectorArchetype&);

    /*! Euclidean norm of a vector.
    *
    * \return The Euclidean norm of this vector.
    */
    RealType norm() const;

private:

    unsigned int n = 0;
    std::unique_ptr<RealType[]> data;

};

/*! Stream vector to string */
template <typename RealType>
std::ostream & operator<<(std::ostream &, const VectorArchetype<RealType>&);

#include "vector_archetype.cpp"

}      // namespace invlib
#endif // ARCHETYPES_VECTOR_ARCHETYPE_H
