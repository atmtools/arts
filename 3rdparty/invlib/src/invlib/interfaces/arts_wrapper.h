#ifndef INTERFACES_ARTS_WRAPPER_H
#define INTERFACES_ARTS_WRAPPER_H

#ifndef HAVE_CONFIG_H
    #define HAVE_CONFIG_H (1)
#endif

#include <iostream>

#include "matpack.h"
#include "matpackI.h"
#include "matpackII.h"
#include "lin_alg.h"

#include "invlib/traits.h"

using invlib::disable_if;
using invlib::is_same;
using invlib::decay;

class ArtsMatrix;
class ArtsVector;
template <typename MatrixType> class ArtsMatrixReference;

// --------------//
//  Arts Vector  //
// ------------- //

class ArtsVector : public Vector
{
public:

    using RealType = Numeric;
    using VectorType = ArtsVector;
    using MatrixType = ArtsMatrix;
    using ResultType = ArtsVector;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    ArtsVector() : Vector() {}
    ArtsVector(const ArtsVector &A) = default;
    ArtsVector & operator=(const ArtsVector &A) = default;
    ArtsVector(const Vector &v) : Vector(v) {}

    ArtsVector(ArtsVector &&A);
    ArtsVector & operator=(ArtsVector &&A);

    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    Index rows() const;

    Numeric operator()(Index i) const;
    Numeric& operator()(Index i);

    Numeric * data_pointer();
    const Numeric * data_pointer() const;

    void accumulate(const ArtsVector& w);
    void subtract(const ArtsVector& w);

    void scale(Numeric c);

    Numeric norm() const;

};

Numeric dot(const ArtsVector& v, const ArtsVector& w);

// --------------//
//  Arts Matrix  //
// ------------- //

/** \brief Arts dense matrix interace wrapper.
 *
 * Simple wrapper class providing an interface to the ARTS matrix class.
 *
 */
class ArtsMatrix : public Matrix
{
public:

    using RealType = Numeric;
    using VectorType = ArtsVector;
    using MatrixType = ArtsMatrix;
    using ResultType = ArtsMatrix;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    ArtsMatrix() : Matrix() {}
    ArtsMatrix (const ArtsMatrix &A)            = default;
    ArtsMatrix & operator=(const ArtsMatrix &A) = default;

    ArtsMatrix (const Matrix &A);

    ArtsMatrix(ArtsMatrix &&A);
    ArtsMatrix & operator=(ArtsMatrix &&A);

    template <typename ArtsType>
    ArtsMatrix(const ArtsMatrixReference<ArtsType> & A);

    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    Index rows() const;
    Index cols() const;

    RealType & operator()(Index i, Index j);
    RealType operator()(Index i, Index j) const;

    RealType * data_pointer();

    // ------------ //
    //  Arithmetic  //
    // ------------ //

    void accumulate(const ArtsMatrix& B);
    void subtract(const ArtsMatrix& B);

    ArtsMatrix multiply(const ArtsMatrix &B) const;
    ArtsVector multiply(const ArtsVector &v) const;

    ArtsMatrix transpose_multiply(const ArtsMatrix &B) const;
    ArtsVector transpose_multiply(const ArtsVector &v) const;
    ArtsVector transpose_multiply_block(const ArtsVector &v,
                                        unsigned int start,
                                        unsigned int extent) const;

    VectorType solve(const VectorType& v) const;
    ArtsMatrix invert() const;

    void scale(Numeric c);

    ArtsMatrix transpose() const;

};

// ------------------------//
//  Arts Matrix Reference  //
// ----------------------- //

/** \brief Wrapper for reference to Arts matrices.
 *
 * Lightweight wrapper object that wraps around an
 * existing (!) Arts matrix. Avoids copying of input
 * and output matrices.
 *
 * \tparam The Arts matrix type to wrap around, i.e.
 * Matrix or Sparse.
 */
template <typename ArtsType>
class ArtsMatrixReference
{
public:

    using RealType = Numeric;
    using VectorType = ArtsVector;
    using MatrixType = ArtsMatrix;
    using ResultType = ArtsMatrix;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    ArtsMatrixReference() = delete;

    template <typename T>
    ArtsMatrixReference(T & A_) : A(A_) {}

    ArtsMatrixReference(const ArtsMatrixReference &) = default;
    ArtsMatrixReference(ArtsMatrixReference&&) = delete;

    ArtsMatrixReference & operator=(ArtsMatrixReference&)
        = default;
    ArtsMatrixReference & operator=(ArtsMatrixReference &&)
        = delete;

    operator ArtsType   & () const {return A.get();}

    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    Index rows() const;
    Index cols() const;

    RealType operator()(unsigned int i, unsigned int j) const;

    // ------------ //
    //  Arithmetic  //
    // ------------ //

    ArtsVector multiply(const ArtsVector &v) const;
    ArtsVector transpose_multiply(const ArtsVector &v) const;
    ArtsMatrix multiply(const ArtsMatrix &B) const;
    ArtsMatrix transpose_multiply(const ArtsMatrix &v) const;

private:

    std::reference_wrapper<ArtsType> A;

};

#include "arts_wrapper.cpp"

#endif // INTERFACES_ARTS_WRAPPER_H
