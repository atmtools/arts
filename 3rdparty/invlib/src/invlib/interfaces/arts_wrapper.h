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
class ArtsSparse;

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

    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    Index rows() const;
    Index cols() const;

    RealType & operator()(Index i, Index j);
    RealType operator()(Index i, Index j) const;

    // ------------ //
    //  Arithmetic  //
    // ------------ //

    void accumulate(const ArtsMatrix& B);
    void accumulate(const ArtsSparse& B);

    void subtract(const ArtsMatrix& B);

    ArtsMatrix multiply(const ArtsMatrix &B) const;
    ArtsVector multiply(const ArtsVector &v) const;

    ArtsMatrix transpose_multiply(const ArtsMatrix &B) const;
    ArtsVector transpose_multiply(const ArtsVector &v) const;

    VectorType solve(const VectorType& v) const;
    ArtsMatrix invert() const;

    void scale(Numeric c);

    ArtsMatrix transpose() const;

};

/** \brief Arts dense matrix interace wrapper.
 *
 * Simple wrapper class providing an interface to the ARTS matrix class.
 *
 */
class ArtsSparse
{
public:

    using RealType = Numeric;
    using VectorType = ArtsVector;
    using MatrixType = ArtsMatrix;
    using ResultType = ArtsMatrix;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    ArtsSparse() = delete;
    ArtsSparse(const Sparse& A_) : A(A_) {}

    ArtsSparse(const ArtsSparse&) = default;
    ArtsSparse(ArtsSparse&&)      = delete;

    ArtsSparse & operator=(const ArtsSparse&) = delete;
    ArtsSparse & operator=(ArtsSparse &&)     = delete;

    // ----------------- //
    //   Manipulations   //
    // ----------------- //

    Index rows() const;
    Index cols() const;

    // ------------ //
    //  Arithmetic  //
    // ------------ //

    ArtsVector multiply(const ArtsVector &v) const;
    ArtsVector transpose_multiply(const ArtsVector &v) const;
    ArtsMatrix multiply(const ArtsMatrix &B) const;

    operator ArtsMatrix() const;

private:

    const Sparse& A;

};

#include "arts_wrapper.cpp"

#endif // INTERFACES_ARTS_WRAPPER_H
