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

// -------------------- //
// Forward Declarations //
// -------------------- //

template
<
typename Real
>
class MatrixArchetype;

template
<
typename Real
>
class VectorArchetype;

template
<
typename Real
>
Real dot(const VectorArchetype<Real>&, const VectorArchetype<Real>&);

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
    using VectorType = VectorArchetype<Real>;
    /*! The fundamental matrix type used for the matrix algebra.*/
    using MatrixType = MatrixArchetype<Real>;
    /*! The result type of multiplying an algebraic expression with this
     * matrix from the right.
     */
    using ResultType = VectorArchetype<Real>;

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
    VectorArchetype(VectorArchetype &&)            = default;

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
    Real & operator()(unsigned int i);

    /*! Read-only element access.
     *
     * \param i Index of the element to access.
     */
    Real   operator()(unsigned int i) const;

    /*! Number of rows of the vector
     *
     * \return The number of rows (dimension) of the vector.
     */
    unsigned int rows() const;

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
    void scale(Real c);

    /*! Dot product of two vectors
     *
     * \return The dot product \f$ \sum_{i=1}^n v_i w_i \f$ of the
     * two given vectors.
     */
    friend Real dot<>(const VectorArchetype&, const VectorArchetype&);

    /*! Euclidean norm of a vector.
    *
    * \return The Euclidean norm of this vector.
    */
    Real norm();

private:

    unsigned int n;
    std::unique_ptr<Real[]> data;

};

/*! Stream vector to string */
template <typename Real>
std::ostream & operator<<(std::ostream &, const VectorArchetype<Real>&);

#include "vector_archetype.cpp"

#endif // ARCHETYPES_VECTOR_ARCHETYPE_H
