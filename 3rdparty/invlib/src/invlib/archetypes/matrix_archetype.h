/**
 * \file archtypes/dense_matrix.h
 *
 * \brief Contains the MatrixArchetype class, which is an archetype for the basic
 * matrix type used.
 *
 */
#ifndef ARCHETYPES_MATRIX_ARCHETYPE
#define ARCHETYPES_MATRIX_ARCHETYPE

#include "invlib/archetypes/vector_archetype.h"
#include <memory>
#include <iterator>

// ------------------------  //
//   Matrix Archetype Class  //
// ------------------------  //

/*! Matrix archtype.
 *
 * Implements a straight forward dense matrix class to verify the
 * generic matrix algebra and illustrate the interface to the
 * fundamental matrix type.
 *
 * \tparam The floating point type used to represent scalars.
 */
template
<
typename Real
>
class MatrixArchetype
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
    using ResultType = MatrixArchetype<Real>;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Default constructor. */
    MatrixArchetype() = default;

    /*! Copy constructor.
     *
     * The copy constructor should implement a deep copy of the argument, i.e.
     * after the call the constructed object hold a matrix identical to the one
     * passed as argument, but that is completely independent of the argument
     */
    MatrixArchetype(const MatrixArchetype &);
    MatrixArchetype(MatrixArchetype &&);

    /*! Assignment operator.
     *
     * The assignment operator should also perform a deep copy of the argument
     * matrix, such that the assigned to object after the call is identical to
     * the provided argument but independent.
     */
    MatrixArchetype& operator=(const MatrixArchetype &);
    MatrixArchetype& operator=(MatrixArchetype &&);

    /*!
     * Frees all resources occupied by the matrix object.
     */
    ~MatrixArchetype() = default;

    MatrixArchetype get_block(unsigned int i,
                              unsigned int j,
                              unsigned int di,
                              unsigned int dj) const;
    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    /*! Resize matrix.
     *
     * Resize the matrix to a \f$i \times j \f$ matrix.
     *
     * \param i Number of rows of the resized matrix.
     * \param j Number of columns of the resized matrix.
     */
    void resize(unsigned int i, unsigned int j);

    /*! Element access.
     * \param i The row of the element to access.
     * \param j The column of the element to access.
     * \return Lvalue reference to the matrix element in row i
     * and column j.
     */
    inline Real & operator()(unsigned int, unsigned int);

    /*! Constant element access.
     * \param i The row of the element to access.
     * \param j The column of the element to access.
     * \return The value of the matrix element in row i and column j.
     */
    inline Real operator()(unsigned int, unsigned int) const;

    /*! Number of columns of the matrix.
     * \return The number of columns of the matrix.
     */
    unsigned int cols() const;

    /*! Number of rows of the matrix.
     * \return The number of rows of the matrix.
     */
    unsigned int rows() const;

    /*! Raw pointer to the matrix data. Needed for MPI testing. */
    RealType * data_pointer();
    const RealType * data_pointer() const;

    // ------------ //
    //  Arithmetic  //
    // ------------ //

    /*! Accumulate into matrix.
     *
     * Add the elements of the given matrix to this matrix object.
     *
     * \param The matrix B to acummulate into this matrix.
     */
    void accumulate(const MatrixArchetype &B);

    /*! Subtract from matrix.
     *
     * Subtract the elements of the given matrix from this matrix object.
     *
     * \param The matrix B to subtract from this matrix.
     */
    void subtract(const MatrixArchetype &);

    /*! Matrix-matrix product.
     *
     * Compute the matrix product \f$C = A B\f$ of this matrix \f$A\f$
     * and the given matrix \f$B\f$.
     *
     * \param B The right hand operator for the multiplication
     * \return The matrix containing the result \f$C\f$ of the matrix
     * matrix product.
     */
    MatrixArchetype multiply(const MatrixArchetype &) const;

    /*! Matrix-vector product.
     *
     * Compute the matrix product \f$w = A v\f$ of this matrix \f$A\f$
     * and the given vector \f$v\f$.
     *
     * \param v The vector to be multiplied from the right.
     * \return The matrix containing the result \f$w\f$ of the matrix-vector
     *  product.
     */
    VectorType multiply(const VectorType &v) const;

    /*! Solve linear system.
     *
     * \return The solution \f$w\f$ of the linear system \f$Aw = v\f$
     * corresponding to this matrix \f$A\f$
     */
    VectorType solve(const VectorType &) const;

    /*! Invert matrix.
     *
     * \return The inverse \f$A^{-1}\f$ of this matrix \f$A\f$.
     */
    MatrixType invert() const;

    /*! Scale matrix.
     *
     * \param c The scalar to scale the matrix with.
     */
    void scale(Real c);

    /*! Transpose matrix.
     *
     * \return A Matrix \f$C = A^T\f$ containing the transpose of the matrix
     * \f$A\f$.
     */
    MatrixType transpose();

    /*! Combined transpose and multiply by matrix.
     *
     * Provided for efficiency reasons, since the transpose operation and the
     * product can be efficiently combined.
     *
     * \return The product \f$C = A^T B\f$ of this matrix \f$A\f$ and the provided
     * argument \f$B\f$.
     */
    MatrixType transpose_multiply(const MatrixType&) const;

    /*! Combined transpose and multiply by vector.
     *
     * Provided for efficiency reasons, since the transpose operation and the
     * product can be efficiently combined.
     *
     * \return The product \f$w = A^T v\f$ of this matrix \f$A\f$ and the provided
     * vector \f$v\f$.
     */
    VectorType transpose_multiply(const VectorType&) const;

    /*! Combined transpose and multiply by part of vector.
     *
     * Needed only when matrices are to be distributed over MPI processes. Computed
     * the contribution of the local matrix to the total product of the distributed
     * matrix.
     *
     * \param block_start The start of the block in the vector v corresponding
     * to the local matrix.
     * \param block_length The size of the block. Must be m.
     */
    VectorType transpose_multiply_block(const VectorType&,
                                        int block_start,
                                        int block_length) const;

    /*! Helper function to compute the QR decomposition of a matrix. Not
     * required by the matrix algebra interface.
     */
    MatrixType QR() const;

    /*! Helper function to backsubstitute a given QR decomposition. Not
     * required by the matrix algebra interface.
     */
    VectorType backsubstitution(const VectorType &) const;

    // --------------------------------- //
    //  Diagonal, Row adn Column Access  //
    // --------------------------------- //

    /*! Return the diagonal of the matrix as a vector.
     *
     *  Required for using the CG solver with Jacobian preconditioner.
     */
    VectorType diagonal() const;

    /*! Return the ith row of the matrix as a vector.
     *
     *  Required for using the CG solver with Jacobian preconditioner.
     */
    VectorType row(size_t i) const;

    /*! Return the ith column of the matrix as a vector.
     *
     *  Required for using the CG solver with Jacobian preconditioner.
     */
    VectorType col(size_t i) const;

private:

    unsigned int m = 0;
    unsigned int n = 0;
    std::unique_ptr<Real[]> data;

};

/*! Stream vector to string */
template <typename Real>
std::ostream & operator<<(std::ostream &, const MatrixArchetype<Real>&);

#include "matrix_archetype.cpp"

#endif
