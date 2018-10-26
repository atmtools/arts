/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   complex.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-12-19
  
  \brief  A class implementing complex numbers for ARTS.
*/

#ifndef complex_h
#define complex_h

#include <complex>
#include "matpackI.h"

typedef std::complex<Numeric> Complex;


std::complex<float> operator+ (const double& d, const std::complex<float>& c);
std::complex<float> operator* (const double& d, const std::complex<float>& c);

std::complex<float> operator+ (const std::complex<float>& c, const double& d);
std::complex<float> operator* (const std::complex<float>& c, const double& d);

std::complex<double> operator+ (const float& f, const std::complex<double>& c);
std::complex<double> operator* (const float& f, const std::complex<double>& c);

std::complex<double> operator+ (const std::complex<double>& c, const float& d);
std::complex<double> operator* (const std::complex<double>& c, const float& d);

// Declare existence of the global joker object:
extern const Joker joker;

// Declare the existence of class ConstComplexMatrixView:
class ConstComplexIterator1D;

// Declare the existence of class ComplexVectorView:
class ComplexVectorView;

// Declare the existence of class ConstComplexVectorView:
class ConstComplexVectorView;

// Declare the existence of class ConstMatrixView:
class ConstComplexMatrixView;

// Eigen library interactions:
typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ComplexMatrixType;
typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> ComplexStrideType;
typedef Eigen::Map<ComplexMatrixType, 0, ComplexStrideType> ComplexMatrixViewMap;
typedef Eigen::Map<const ComplexMatrixType, 0, ComplexStrideType> ComplexConstMatrixViewMap;

/** The iterator class for sub vectors. This takes into account the
 *  defined stride. */
class ComplexIterator1D {
public:
    /** Default constructor. */
    ComplexIterator1D() : mx(NULL), mstride(0) { /* Nothing to do here. */ }
    
    /** Explicit constructor. */
    ComplexIterator1D(Complex *x, Index stride) : mx(x), mstride(stride)
    { /* Nothing to do here. */ }
    
    // Operators:
    
    /** Prefix increment operator. */
    ComplexIterator1D& operator++()
    { mx += mstride; return *this; }
    
    /** Dereferencing. */
    Complex& operator*() const { return *mx; }
    
    /** Not equal operator, needed for algorithms like copy. */
    bool operator!=(const ComplexIterator1D& other) const
    { if (mx != other.mx) return true; else return false; }
    
    friend void copy(ConstComplexIterator1D origin,
                     const ConstComplexIterator1D& end,
                     ComplexIterator1D target);
    
private:
    /** Current position. */
    Complex *mx;
    /** Stride. */
    Index mstride;
};

/** The constant iterator class for sub vectors. This takes into
 *  account the defined stride. */
class ConstComplexIterator1D {
public:
    /** Default constructor. */
    ConstComplexIterator1D() : mx(NULL), mstride(0)
    { /* Nothing to do here. */ }
    
    /** Explicit constructor. */
    ConstComplexIterator1D(Complex *x, Index stride) : mx(x), mstride(stride)
    { /* Nothing to do here. */ }
    
    // Operators:
    /** Prefix increment operator. */
    ConstComplexIterator1D& operator++()
    { mx += mstride; return *this; }
    
    /** Dereferencing. */
    const Complex& operator*() const { return *mx; }
    
    /** Not equal operator, needed for algorithms like copy. */
    bool operator!=(const ConstComplexIterator1D& other) const
    { if (mx != other.mx) return true; else return false; }
    
    friend void copy(ConstComplexIterator1D origin,
                     const ConstComplexIterator1D& end,
                     ComplexIterator1D target);
private:
    /** Current position. */
    const Complex *mx;
    /** Stride. */
    Index mstride;
};

// Declare the complex vector class:
class ComplexVector;

// Declare the ComplexMatrixView class
class ComplexMatrixView;

/** A constant view of a ComplexVector.
 * 
 Together with the derived class ComplexVectorView this contains the main
 implementation of a ComplexVector. The class ComplexVector is just a special
 case of a ComplexVectorView which also allocates storage. */
class ConstComplexVectorView {
public:
    // Typedef for compatibility with STL
    typedef ConstComplexIterator1D const_iterator;
    
    // Member functions:
    bool empty() const;
    Index nelem() const;
    Complex sum() const;
    
    // Const index operators:
    /** Plain const index operator. */
    Complex operator[](Index n) const
    { // Check if index is valid:
        assert( 0<=n );
        assert( n<mrange.mextent );
        return get(n);
    }
    
//     /** Get element implementation without assertions. */
//     Complex get(Index n) const
//     {
//         return *( mdata +
//         mrange.mstart +
//         n*mrange.mstride );
//     }
//     
//     /** Get element implementation without assertions. */
//     Numeric get_real(Index n) const { return get(n).real(); }
//     
//     /** Get element implementation without assertions. */
//     Numeric get_imag(Index n) const { return get(n).imag(); }
    
    /** Get element implementation without assertions. */
    const Complex& get(Index n) const
    {
      return *( mdata +
      mrange.mstart +
      n*mrange.mstride );
    }
    
    /** Get element implementation without assertions. */
    const Numeric& get_real(Index n) const
    { return reinterpret_cast<const Numeric (&)[2]>(get(n))[0]; }
    
    /** Get element implementation without assertions. */
    const Numeric& get_imag(Index n) const
    { return reinterpret_cast<const Numeric (&)[2]>(get(n))[1]; }
    
    ConstComplexVectorView operator[](const Range& r) const;
    friend Complex operator*(const ConstComplexVectorView& a, const ConstComplexVectorView& b);
    
    // Functions returning iterators:
    ConstComplexIterator1D begin() const;
    ConstComplexIterator1D end() const;
    
    // Conversion to 1 column matrix:
    operator ConstComplexMatrixView() const;
    
    //! Destructor
    virtual ~ConstComplexVectorView() {}
    
    // Friends:
    friend class ComplexVectorView;
    friend class ConstComplexIterator2D;
    friend class ConstComplexMatrixView;
    friend void mult (ComplexVectorView,
                      const ConstComplexMatrixView &,
                      const ConstComplexVectorView &);
    
    friend void diagonalize( ComplexMatrixView,
                             ComplexVectorView,
                             const ConstComplexMatrixView&);
    
    friend ComplexConstMatrixViewMap MapToEigen(const ConstComplexVectorView&);
    friend ComplexConstMatrixViewMap MapToEigenCol(const ConstComplexVectorView&);
    friend ComplexMatrixViewMap MapToEigen(ComplexVectorView&);
    friend ComplexMatrixViewMap MapToEigenCol(ComplexVectorView&);
    
    // A special constructor, that allows to make a ConstVectorView of a scalar.
    ConstComplexVectorView(const Complex& a);
    
protected:
    // Constructors:
    ConstComplexVectorView();
    ConstComplexVectorView(Complex *data, const Range& range);
    ConstComplexVectorView(Complex *data, const Range& p, const Range& n);
    
    // Data members:
    // -------------
    /** The range of mdata that is actually used. */
    Range mrange;
    /** Pointer to the plain C array that holds the data */
    Complex *mdata;
};

/** The ComplexVectorView class.
 * 
 This contains the main implementation of a complex vector. The class     
 ComplexVector is just a special case of subvector which also allocates
 storage.
 
 Unfortunately, names of element functions of derived classes hide
 the names of the original class, even if the arguments are
 different. This means that we have to redefine those element
 functions that can have different arguments, for example the
 constant index operators and iterators. */
class ComplexVectorView : public ConstComplexVectorView {
public:
    ComplexVectorView (const ComplexVector&);
    ComplexVectorView (ComplexVector& v);
    
    // Typedef for compatibility with STL
    typedef ComplexIterator1D iterator;
    
    // Const index operators:
    /** Plain const index operator. Has to be redifined here, because the
     one from ConstVectorView is hidden. */                           
     Complex operator[](Index n) const
     { return ConstComplexVectorView::operator[](n); }
     
     /** Get element implementation without assertions. */
     const Complex& get(Index n) const
     { return ConstComplexVectorView::get(n); }
     const Numeric& get_real(Index n) const
     { return ConstComplexVectorView::get_real(n); }
     const Numeric& get_imag(Index n) const
     { return ConstComplexVectorView::get_imag(n); }
     
     ConstComplexVectorView operator[](const Range& r) const;
     
     /** Plain Index operator. */
     Complex& operator[](Index n)
     { // Check if index is valid:
         assert( 0<=n );
         assert( n<mrange.mextent );
         return get(n);
     }
     
     /** Get element implementation without assertions. */
     Complex& get(Index n)
     {
         return *( mdata +
         mrange.mstart +
         n*mrange.mstride );
     }
     
     /** Get element implementation without assertions. */
     Numeric& get_real(Index n) 
     { return reinterpret_cast<Numeric (&)[2]>(get(n))[0]; }
     
     /** Get element implementation without assertions. */
     Numeric& get_imag(Index n)
     { return reinterpret_cast<Numeric (&)[2]>(get(n))[1]; }
     
     ComplexVectorView operator[](const Range& r);
     
     // Constant iterators:
     ConstComplexIterator1D begin() const;
     ConstComplexIterator1D end() const;
     // ComplexIterators:
     ComplexIterator1D begin();
     ComplexIterator1D end();
     
     // Assignment operators:
     ComplexVectorView& operator=(const ConstComplexVectorView& v);
     ComplexVectorView& operator=(const ComplexVectorView& v);
     ComplexVectorView& operator=(const ComplexVector& v);
     ComplexVectorView& operator=(const Array<Complex>& v); 
     ComplexVectorView& operator=(Complex x);
     ComplexVectorView& operator=(const ConstVectorView& v);
     ComplexVectorView& operator=(const VectorView& v);
     ComplexVectorView& operator=(const Vector& v);
     ComplexVectorView& operator=(const Array<Numeric>& v); 
     ComplexVectorView& operator=(Numeric x);
     
     // Other operators:
     ComplexVectorView operator*=(Complex x);
     ComplexVectorView operator/=(Complex x);
     ComplexVectorView operator+=(Complex x);
     ComplexVectorView operator-=(Complex x);
     ComplexVectorView operator*=(Numeric x);
     ComplexVectorView operator/=(Numeric x);
     ComplexVectorView operator+=(Numeric x);
     ComplexVectorView operator-=(Numeric x);
     
     ComplexVectorView operator*=(const ConstComplexVectorView& x);
     ComplexVectorView operator/=(const ConstComplexVectorView& x);
     ComplexVectorView operator+=(const ConstComplexVectorView& x);
     ComplexVectorView operator-=(const ConstComplexVectorView& x);
     ComplexVectorView operator*=(const ConstVectorView& x);
     ComplexVectorView operator/=(const ConstVectorView& x);
     ComplexVectorView operator+=(const ConstVectorView& x);
     ComplexVectorView operator-=(const ConstVectorView& x);
     
     // Conversion to 1 column matrix:
     operator ComplexMatrixView();
     // Conversion to a plain C-array
     const Complex *get_c_array() const;
     Complex *get_c_array();
     
     // Vector
     Vector real()
     {
         Vector A(nelem());
         for(Index i=0;i<nelem();i++)
           A[i]=get_real(i);
         return A;
     };
     Vector imag()
     {
         Vector A(nelem());
         for(Index i=0;i<nelem();i++)
           A[i]=get_imag(i);
         return A;
     };
     
     //! Destructor
     virtual ~ComplexVectorView() {}
     
     // Friends:
     friend class ConstComplexIterator2D;
     friend class ComplexIterator2D;
     friend class ComplexMatrixView;
     
     // A special constructor, that allows to make a ComplexVectorView of a scalar.
     ComplexVectorView(Complex& a);
     
     
protected:
    // Constructors:
    ComplexVectorView();
    ComplexVectorView(Complex *data, const Range& range);
    ComplexVectorView(Complex *data, const Range& p, const Range& n);
};

/** The row iterator class for sub matrices. This takes into account the
 defined row stride. The iterator points to a row of the matrix,  *
 which acts just like a ComplexVectorView. */
class ComplexIterator2D {
public:
    // Constructors:
    /** Default constructor. */
    ComplexIterator2D() : msv(), mstride(0)  { /* Nothing to do here. */ }
    
    /** Explicit constructor. */
    ComplexIterator2D(const ComplexVectorView& x, Index stride) : msv(x), mstride(stride)
    { /* Nothing to do here. */ }
    
    // Operators:
    /** Prefix increment operator. */
    ComplexIterator2D& operator++() { msv.mdata += mstride; return *this; }
    
    /** Not equal operator, needed for algorithms like copy. */
    bool operator!=(const ComplexIterator2D& other) const
    { if ( msv.mdata + msv.mrange.mstart !=
        other.msv.mdata + other.msv.mrange.mstart )
    return true;
    else
        return false;
    }
    
    /** The -> operator is needed, so that we can write i->begin() to get
     the 1D iterators. */                                             
     ComplexVectorView * operator->() { return &msv; }
     
     /** Dereferencing. */
     ComplexVectorView& operator*() { return msv; }
     
private:
    /** Current position. */
    ComplexVectorView msv;
    /** Row stride. */
    Index mstride;
};

/** The const row iterator class for sub matrices. This takes into account the
 defined row stride. The iterator points to a row of the matrix,  *
 which acts just like a ComplexVectorView. */
class ConstComplexIterator2D {
public:
    // Constructors:
    /** Default constructor. */
    ConstComplexIterator2D() : msv(), mstride(0) { /* Nothing to do here. */ }
    
    /** Explicit constructor. */
    ConstComplexIterator2D(const ConstComplexVectorView& x, Index stride)
    : msv(x), mstride(stride)
    { /* Nothing to do here. */ }
    
    // Operators:
    /** Prefix increment operator. */
    ConstComplexIterator2D& operator++() { msv.mdata += mstride; return *this; }
    
    /** Not equal operator, needed for algorithms like copy. */
    bool operator!=(const ConstComplexIterator2D& other) const
    { if ( msv.mdata + msv.mrange.mstart !=
        other.msv.mdata + other.msv.mrange.mstart )
    return true;
    else
        return false;
    }
    
    /** The -> operator is needed, so that we can write i->begin() to get
     t he 1D iterators. */                                             
     const ConstComplexVectorView* operator->() const { return &msv; }
     
     /** Dereferencing. */
     const ConstComplexVectorView& operator*() const { return msv; }
     
private:
    /** Current position. */
    ConstComplexVectorView msv;
    /** Row stride. */
    Index mstride;
};

/** The ComplexVector class. This is a subvector that also allocates storage
 automatically, and deallocates it when it is destroyed. We take  *
 all the functionality from ComplexVectorView. Additionally defined in
 this class are:
 
 1. Constructors and destructors (allocating memory).
 2. Assignment operator
 3. Assignment operator from scalar.
 4. Resize function.
 */
class ComplexVector : public ComplexVectorView {
public:
    // Constructors:
    ComplexVector();
    explicit ComplexVector(Index n);
    ComplexVector(Index n, Complex fill);
    ComplexVector(Index n, Numeric fill);
    ComplexVector(Complex start, Index extent, Complex stride);
    ComplexVector(Complex start, Index extent, Numeric stride);
    ComplexVector(Numeric start, Index extent, Complex stride);
    ComplexVector(Numeric start, Index extent, Numeric stride);
    ComplexVector(const ConstComplexVectorView& v);
    ComplexVector(const ComplexVector& v);
    ComplexVector(const Vector& v);
    ComplexVector(const std::vector<Complex>&);
    ComplexVector(const std::vector<Numeric>&);
    
    // Assignment operators:
    ComplexVector& operator=(ComplexVector v);
    ComplexVector& operator=(const Array<Complex>& v);
    ComplexVector& operator=(Complex x);
    
    // Resize function:
    void resize(Index n);
    
    // Swap function:
    friend void swap(ComplexVector& v1, ComplexVector& v2);
    
    // Destructor:
    virtual ~ComplexVector();
};

// Declare class ComplexMatrix:
class ComplexMatrix;


/** A constant view of a ComplexMatrix.
 * 
 This, together with the derived class Complex*MatrixView, contains the   *
 main implementation of a ComplexMatrix. It defines the concepts of
 ComplexMatrixView. Plus additionally the recursive subrange operator,
 which makes it possible to create a ComplexMatrixView from a subrange of
 a ComplexMatrixView.
 
 The class ComplexMatrix is just a special case of a ComplexMatrixView
 which also allocates storage. */
class ConstComplexMatrixView {
public:
    // Typedef for compatibility with STL
    typedef ConstComplexIterator2D const_iterator;
    
    // Member functions:
    bool empty() const;
    Index nrows() const;
    Index ncols() const;
    
    // Const index operators:
    /** Plain const index operator. */
    Complex operator()(Index r, Index c) const
    { // Check if indices are valid:
        assert( 0<=r );
        assert( 0<=c );
        assert( r<mrr.mextent );
        assert( c<mcr.mextent );
        
        return get(r, c);
    }
    
    /** Get element implementation without assertions. */
    Complex get(Index r, Index c) const
    {
        return *( mdata +
        mrr.mstart +
        r*mrr.mstride +
        mcr.mstart +
        c*mcr.mstride );
    }
    
    /** Get element implementation without assertions. */
    Numeric get_real(Index r, Index c) const { return get(r, c).real(); }
    
    /** Get element implementation without assertions. */
    Numeric get_imag(Index r, Index c) const { return get(r, c).imag(); }
    
    ConstComplexMatrixView operator()(const Range& r, const Range& c) const;
    ConstComplexVectorView operator()(const Range& r, Index c) const;
    ConstComplexVectorView operator()(Index r, const Range& c) const;
    
    // Functions returning iterators:
    ConstComplexIterator2D begin() const;
    ConstComplexIterator2D end() const;
    
    // View on diagonal complex vector
    ConstComplexVectorView diagonal() const;
    
    //! Destructor
    virtual ~ConstComplexMatrixView() {}
    
    // Friends:
    friend class ComplexMatrixView;
    friend class ConstComplexVectorView;
    friend ConstComplexMatrixView transpose(ConstComplexMatrixView m);
    friend void mult (ComplexVectorView,
                      const ConstComplexMatrixView &,
                      const ConstComplexVectorView &);
    friend void mult (ComplexMatrixView,
                      const ConstComplexMatrixView &,
                      const ConstComplexMatrixView &);
    friend void mult (ComplexMatrixView,
                      const ConstMatrixView &,
                      const ConstComplexMatrixView &);
    friend void mult (ComplexMatrixView,
                      const ConstComplexMatrixView &,
                      const ConstMatrixView &);
    
    friend void inv( ComplexMatrixView,
                     const ConstComplexMatrixView&);
    friend void diagonalize( ComplexMatrixView,
                             ComplexVectorView,
                             const ConstComplexMatrixView&);
    
    friend ComplexConstMatrixViewMap MapToEigen(const ConstComplexMatrixView&);
    friend ComplexMatrixViewMap MapToEigen(ComplexMatrixView&);
    
protected:
    // Constructors:
    ConstComplexMatrixView();
    ConstComplexMatrixView(Complex *data, const Range& r, const Range& c);
    ConstComplexMatrixView(Complex *data,
                    const Range& pr, const Range& pc,
                    const Range& nr, const Range& nc);
    
    // Data members:
    // -------------
    /** The row range of mdata that is actually used. */
    Range mrr;
    /** The column range of mdata that is actually used. */
    Range mcr;
    /** Pointer to the plain C array that holds the data */
    Complex *mdata;
};

/** The ComplexMatrixView class
 * 
 This contains the main implementation of a ComplexMatrix. It defines    
 the concepts of ComplexMatrixView. Plus additionally the recursive
 subrange operator, which makes it possible to create a ComplexMatrixView
 from a subrange of a ComplexMatrixView. 
 
 The class ComplexMatrix is just a special case of a ComplexMatrixView
 which also allocates storage. */
class ComplexMatrixView : public ConstComplexMatrixView {
public:
    // Typedef for compatibility with STL
    typedef ComplexIterator2D iterator;
    
    // Const index operators:
    /** Plain const index operator. Has to be redefined here, since it is
     h iden by the non-const operator of the derived class. */         
     Complex operator()(Index r, Index c) const
     { return ConstComplexMatrixView::operator()(r,c); }
     
     /** Get element implementation without assertions. */
     Complex get(Index r, Index c) const
     { return ConstComplexMatrixView::get(r,c); }
     Numeric get_real(Index r, Index c) const
     { return ConstComplexMatrixView::get_real(r,c); }
     Numeric get_imag(Index r, Index c) const
     { return ConstComplexMatrixView::get_imag(r,c); }
     
     ConstComplexMatrixView operator()(const Range& r, const Range& c) const;
     ConstComplexVectorView operator()(const Range& r, Index c) const;
     ConstComplexVectorView operator()(Index r, const Range& c) const;
     // Index Operators:
     /** Plain index operator. */
     Complex& operator()(Index r, Index c)
     { // Check if indices are valid:
         assert( 0<=r );
         assert( 0<=c );
         assert( r<mrr.mextent );
         assert( c<mcr.mextent );
         
         return get(r, c);
     }
     
     /** Get element implementation without assertions. */
     Complex& get(Index r, Index c)
     {
         return *( mdata +
         mrr.mstart +
         r*mrr.mstride +
         mcr.mstart +
         c*mcr.mstride );
     }
     
     /** Get element implementation without assertions. */
     Numeric& get_real(Index r, Index c)
     { return reinterpret_cast<Numeric (&)[2]>(get(r, c))[0]; }
     
     /** Get element implementation without assertions. */
     Numeric& get_imag(Index r, Index c)
     { return reinterpret_cast<Numeric (&)[2]>(get(r, c))[1]; }
     
     ComplexMatrixView operator()(const Range& r, const Range& c);
     ComplexVectorView operator()(const Range& r, Index c);
     ComplexVectorView operator()(Index r, const Range& c);
     
     // Functions returning const iterators:
     ConstComplexIterator2D begin() const;
     ConstComplexIterator2D end() const;
     // Functions returning iterators:
     ComplexIterator2D begin();
     ComplexIterator2D end();
     
     // Assignment operators:
     ComplexMatrixView& operator=(const ConstComplexMatrixView& v);
     ComplexMatrixView& operator=(const ComplexMatrixView& v);
     ComplexMatrixView& operator=(const ComplexMatrix& v);
     ComplexMatrixView& operator=(const ConstComplexVectorView& v);
     ComplexMatrixView& operator=(Complex x);
     
     // Other operators:
     ComplexMatrixView& operator*=(Complex x);
     ComplexMatrixView& operator/=(Complex x);
     ComplexMatrixView& operator+=(Complex x);
     ComplexMatrixView& operator-=(Complex x);
     ComplexMatrixView& operator*=(Numeric x);
     ComplexMatrixView& operator/=(Numeric x);
     ComplexMatrixView& operator+=(Numeric x);
     ComplexMatrixView& operator-=(Numeric x);
     
     ComplexMatrixView& operator*=(const ConstComplexMatrixView& x);
     ComplexMatrixView& operator/=(const ConstComplexMatrixView& x);
     ComplexMatrixView& operator+=(const ConstComplexMatrixView& x);
     ComplexMatrixView& operator-=(const ConstComplexMatrixView& x);
     
     ComplexMatrixView& operator*=(const ConstMatrixView& x);
     ComplexMatrixView& operator/=(const ConstMatrixView& x);
     ComplexMatrixView& operator+=(const ConstMatrixView& x);
     ComplexMatrixView& operator-=(const ConstMatrixView& x);
     
     ComplexMatrixView& operator*=(const ConstComplexVectorView& x);
     ComplexMatrixView& operator/=(const ConstComplexVectorView& x);
     ComplexMatrixView& operator+=(const ConstComplexVectorView& x);
     ComplexMatrixView& operator-=(const ConstComplexVectorView& x);
     
     // Conversion to a plain C-array
     const Complex *get_c_array() const;
     Complex *get_c_array();
     
     // Matrix
     Matrix real()
     {
         Matrix A(nrows(),ncols());
         for(Index i=0;i<nrows();i++)
             for(Index j=0;j<ncols();j++)
                 A(i,j)=get_real(i,j);
         return A;
     };
     Matrix imag()
     {
         Matrix A(nrows(),ncols());
         for(Index i=0;i<nrows();i++)
             for(Index j=0;j<ncols();j++)
                 A(i,j)=get_imag(i,j);
         return A;
     };
     
     //! Destructor
     virtual ~ComplexMatrixView() {}
     
     // Friends:
     friend class ComplexVectorView;
     friend ConstComplexMatrixView transpose(ConstComplexMatrixView m);
     friend ComplexMatrixView transpose(ComplexMatrixView m);
     
protected:
    // Constructors:
    ComplexMatrixView();
    ComplexMatrixView(Complex *data, const Range& r, const Range& c);
    ComplexMatrixView(Complex *data,
               const Range& pr, const Range& pc,
               const Range& nr, const Range& nc);
};

/** The ComplexMatrix class. This is a ComplexMatrixView that also allocates storage
 automatically, and deallocates it when it is destroyed. We take  *
 all the functionality from ComplexMatrixView. Additionally defined here
 are: 
 
 1. Constructors and destructor.
 2. Assignment operator from scalar.
 3. Resize function. */
class ComplexMatrix : public ComplexMatrixView {
public:
    // Constructors:
    ComplexMatrix();
    ComplexMatrix(Index r, Index c);
    ComplexMatrix(Index r, Index c, Complex fill);
    ComplexMatrix(Index r, Index c, Numeric fill);
    ComplexMatrix(const ConstComplexMatrixView& v);
    ComplexMatrix(const ComplexMatrix& v);
    
    // Assignment operators:
    ComplexMatrix& operator=(ComplexMatrix x);
    ComplexMatrix& operator=(Complex x);
    ComplexMatrix& operator=(const ConstComplexVectorView& v);
    
    // Resize function:
    void resize(Index r, Index c);
    
    // Swap function:
    friend void swap(ComplexMatrix& m1, ComplexMatrix& m2);
    
    // Destructor:
    virtual ~ComplexMatrix();
    
    Complex *get_raw_data() { return mdata; }
};

// Function declarations:
// ----------------------

ConstComplexMatrixView transpose(ConstComplexMatrixView m);

ComplexMatrixView transpose(ComplexMatrixView m);

void copy(ConstComplexIterator1D origin,
          const ConstComplexIterator1D& end,
          ComplexIterator1D target);

void copy(Complex x,
          ComplexIterator1D target,
          const ComplexIterator1D& end);

void copy(ConstComplexIterator2D origin,
          const ConstComplexIterator2D& end,
          ComplexIterator2D target);

void copy(Complex x,
          ComplexIterator2D target,
          const ComplexIterator2D& end);

void mult( ComplexVectorView y,
           const ConstComplexMatrixView& M,
           const ConstComplexVectorView& x );


void mult( ComplexMatrixView A,
           const ConstComplexMatrixView& B,
           const ConstComplexMatrixView& C );
void mult( ComplexMatrixView A,
           const ConstMatrixView& B,
           const ConstComplexMatrixView& C );
void mult( ComplexMatrixView A,
           const ConstComplexMatrixView& B,
           const ConstMatrixView& C );
void mult( ComplexMatrixView A,
           const ConstMatrixView& B,
           const ConstMatrixView& C );

Complex operator*(const ConstComplexVectorView& a, const ConstComplexVectorView& b);

std::ostream& operator<<(std::ostream& os, const ConstComplexVectorView& v);

std::ostream& operator<<(std::ostream& os, const ConstComplexMatrixView& v);

// Converts constant matrix to constant eigen map
ComplexConstMatrixViewMap MapToEigen(const ConstComplexMatrixView& A);
// Converts constant vector to constant eigen row-view
ComplexConstMatrixViewMap MapToEigen(const ConstComplexVectorView& A);
// Converts constant vector to constant eigen row-view
ComplexConstMatrixViewMap MapToEigenRow(const ConstComplexVectorView& A);
// Converts constant vector to constant eigen column-view
ComplexConstMatrixViewMap MapToEigenCol(const ConstComplexVectorView& A);
// Converts matrix to eigen map
ComplexMatrixViewMap MapToEigen(ComplexMatrixView& A);
// Converts vector to eigen map row-view
ComplexMatrixViewMap MapToEigen(ComplexVectorView& A);
// Converts vector to eigen map row-view
ComplexMatrixViewMap MapToEigenRow(ComplexVectorView& A);
// Converts vector to eigen map column-view
ComplexMatrixViewMap MapToEigenCol(ComplexVectorView& A);

typedef Array<ComplexVector> ArrayOfComplexVector;
typedef Array<ComplexMatrix> ArrayOfComplexMatrix;

#endif

