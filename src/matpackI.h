/* Copyright (C) 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

/**
  Implementation of Matrix, Vector, and such stuff.

  A VectorView consists of the data, which is stored in a continuous piece
  of memory, and a selection, specified by start, extend, and
  stride. A Vector is a VectorView which also allocates its memory
  automatically. 

  VectorViews can not be generated directly, they only result from
  operations on Vectors, such as using the index operator with a Range
  object. However, you can store them, like:

  VectorView A = B[Range(0,3)]

  A VectorView acts like a reference to the selected region in the
  parent matrix. Functions that operate on an existing matrix (i.e.,
  they do not use resize) should take VectorView x as arguement, rather
  than Vector& x. That has the advantage that they can be called
  either with a VectorView or Vector. E.g., if you have:

  void fill_with_junk(VectorView x);
  Vector v;

  then you can call the function in these two ways:

  fill_with_junk(v);
  fill_with_junk(v[Range(1,3)]) 

  Assignment (=) copies the data from one Vector or VectorView to
  another one. Dimensions must agree. Only exceptions are the copy
  constructors which automatically set the dimensions to
  match.

  Things work in the same way for the type Matrix. 

  There exist operators *=, /=, +=, and -= to multiply  (divide,...)
  by a scalar. Plain operators *,... do not exist, because they would
  result in the creation of temporaries and therefor be inefficient.  

  However, you can use * to compute the scalar product. This is
  efficient, since the return value is just a scalar. 

  There is a constructor for vector filling it with a sequence of
  values. 

  Matrices:

  You can extract sub matrices (MatrixView) using Range objects. You
  can also extract rows and columns this way.

  transpose(A) Returns a special MatrixView that is the transpose of the
  original. The original is not changed by this!

  mult(A,B,C) computes A = B*C
  Note that the order is different from MTL (output first)!

  A VectorView or Vector can be taken in the place of a nx1
  matrix. That means, Vectors are interpreted as column
  vectors. Hence, you can compute:

  Vector a(10),b(20);
  Matrix M(10,20);

  mult(a,M,b);    // a = M*b

  but also, by using transpose:

  mult(transpose(b),transpose(a),M);    // b^t = a^t * M

  See the section about Matrices and Vectors in the ARTS user guide
  for more details.
  
  \author Stefan Buehler
  \date   2001-06-12
 */

#ifndef matpackI_h
#define matpackI_h

#include <iomanip>
#include "arts.h"

/** The Joker class.

    This class is used by Vector and Matrix in connection with Range
    to implement Matlab-like subranges of vectors and matrices.

    This class has no members. We just need a special type to indicate
    the joker. There is a global joker object defined somewhere:
    
    Joker joker;
*/
class Joker {
  // Nothing here.
};

// Declare existence of the global joker object:
extern Joker joker;

// Declare the existence of class VectorView:
class VectorView;

/** The range class. 

    This is used to specifiy a range of a vector. In general, a range is
    given by a start index, an extent, and a stride. The entire vector
    would be:
    start = 0, range = # elements, stride = 1

    Stride specifies the stepsize of the vector. A stride of 2 means
    only every second element. This is particularly important in
    connection with matrices.

    There are a number of special constructors for this class, of
    particular interest should be those using jokers, which provide a
    Matlab-like functionality.
*/
class Range{
public:
  // Constructors:
  Range(Index start, Index extent, Index stride=1);
  Range(Index start, Joker joker, Index stride=1);
  Range(Joker joker, Index stride=1);
  Range(Index max_size, const Range& r);
  Range(const Range& p, const Range& n);

  // Friends:
  friend class ConstVectorView;
  friend class VectorView;
  friend class Vector;
  friend class ConstMatrixView;
  friend class MatrixView;
  friend class Matrix;
  friend class Iterator2D;
  friend class Iterator3D;
  friend class ConstIterator2D;
  friend class ConstIterator3D;
  friend class SparseMatrixView;
  friend class ConstTensor3View;
  friend class Tensor3View;
  friend class Tensor3;
  friend std::ostream& operator<<(std::ostream& os, const SparseMatrixView& v);


private:
  /** The start index. */
  Index mstart;
  /** The number of elements. -1 means extent to the end of the
      vector. */
  Index mextent;
  /** The stride. Can be positive or negative. */
  Index mstride;
};

/** The iterator class for sub vectors. This takes into account the
    defined stride. */
class Iterator1D {
public:
  // Constructors:
  Iterator1D();
  Iterator1D(const Iterator1D& o);
  Iterator1D(Numeric *x, Index stride);

  // Operators:
  Iterator1D& operator++();
  Numeric& operator*() const;
  bool operator!=(const Iterator1D& other) const;

private:
  /** Current position. */
  Numeric *mx;
  /** Stride. */
  Index mstride;
};

/** The constant iterator class for sub vectors. This takes into
    account the defined stride. */
class ConstIterator1D {
public:
  // Constructors:
  ConstIterator1D();
  ConstIterator1D(const ConstIterator1D& o);
  ConstIterator1D(Numeric *x, Index stride);

  // Operators:
  ConstIterator1D& operator++();
  const Numeric& operator*() const;
  bool operator!=(const ConstIterator1D& other) const;

private:
  /** Current position. */
  const Numeric *mx;
  /** Stride. */
  Index mstride;
};

// Declare the vector class:
class Vector;

// Declare the MatrixView class
class MatrixView;

/** A constant view of a Vector.

    Together with the derived class VectorView this contains the main
    implementation of a Vector. The class Vector is just a special
    case of a VectorView which also allocates storage. */
class ConstVectorView {
public:

  // Member functions:
  Index nelem() const;
  Numeric sum() const;

  // Const index operators:
  Numeric  operator[](Index n) const;
  ConstVectorView operator[](const Range& r) const;

  // Functions returning iterators:
  ConstIterator1D begin() const;
  ConstIterator1D end() const;
  
  // Conversion to 1 column matrix:
  operator ConstMatrixView() const;

  // Friends:
  friend class VectorView;
  friend class ConstIterator2D;
  friend class ConstMatrixView;
  friend class ConstTensor3View;

protected:
  // Constructors:
  ConstVectorView();
  ConstVectorView(Numeric *data, const Range& range);
  ConstVectorView(Numeric *data, const Range& p, const Range& n);

  // Data members:
  // -------------
  /** The range of mdata that is actually used. */
  Range mrange;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The VectorView class.

    This contains the main implementation of a vector. The class
    Vector is just a special case of subvector which also allocates
    storage. 

    Unfortunately, names of element functions of derived classes hide
    the names of the original class, even if the arguments are
    different. This means that we have to redefine those element
    functions that can have different arguments, for example the
    constant index operators and iterators. */
class VectorView : public ConstVectorView {
public:

  // Const index operators:
  Numeric  operator[](Index n) const;
  ConstVectorView operator[](const Range& r) const;
  // Index Operators:
  Numeric& operator[](Index n);
  VectorView operator[](const Range& r);

  // Constant iterators:
  ConstIterator1D begin() const;
  ConstIterator1D end() const;
  // Iterators:
  Iterator1D begin();
  Iterator1D end();
  
  // Assignment operators:
  VectorView operator=(const ConstVectorView& v);
  VectorView operator=(const VectorView& v);
  VectorView operator=(const Vector& v);
  VectorView operator=(const Array<Numeric>& v);
  VectorView operator=(Numeric x);

  // Other operators:
  VectorView operator*=(Numeric x);
  VectorView operator/=(Numeric x);
  VectorView operator+=(Numeric x);
  VectorView operator-=(Numeric x);

  VectorView operator*=(const ConstVectorView& x);
  VectorView operator/=(const ConstVectorView& x);
  VectorView operator+=(const ConstVectorView& x);
  VectorView operator-=(const ConstVectorView& x);

  // Conversion to 1 column matrix:
  operator MatrixView();

  // Friends:
  friend class ConstIterator2D;
  friend class Iterator2D;
  friend class MatrixView;
  friend class Tensor3View;

protected:
  // Constructors:
  VectorView();
  VectorView(Numeric *data, const Range& range);
  VectorView(Numeric *data, const Range& p, const Range& n);
};

/** The row iterator class for sub matrices. This takes into account the
    defined row stride. The iterator points to a row of the matrix,
    which acts just like a VectorView. */
class Iterator2D {
public:
  // Constructors:
  Iterator2D();
  Iterator2D(const Iterator2D& o);
  Iterator2D(const VectorView& x, Index stride);

  // Operators:
  Iterator2D& operator++();
  bool operator!=(const Iterator2D& other) const;
  VectorView* const operator->();
  VectorView& operator*();
  
private:
  /** Current position. */
  VectorView msv;
  /** Row stride. */
  Index mstride;
};

/** The const row iterator class for sub matrices. This takes into account the
    defined row stride. The iterator points to a row of the matrix,
    which acts just like a VectorView. */
class ConstIterator2D {
public:
  // Constructors:
  ConstIterator2D();
  ConstIterator2D(const ConstIterator2D& o);
  ConstIterator2D(const ConstVectorView& x, Index stride);

  // Operators:
  ConstIterator2D& operator++();
  bool operator!=(const ConstIterator2D& other) const;
  const ConstVectorView* operator->() const;
  const ConstVectorView& operator*() const;

private:
  /** Current position. */
  ConstVectorView msv;
  /** Row stride. */
  Index mstride;
};

/** The Vector class. This is a subvector that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from VectorView. Additionally defined in
    this class are:

    1. Constructors and destructors (allocating memory).
    2. Assignment operator from scalar.
    3. Resize function.

    It seems that we need no assignment operator from VectorView or
    Vector, because the one of VectorView is used. */
class Vector : public VectorView {
public:
  // Constructors:
  Vector();
  explicit Vector(Index n);
  Vector(Index n, Numeric fill);
  Vector(Numeric start, Index extent, Numeric stride);
  Vector(const ConstVectorView& v);
  Vector(const Vector& v);

  // Assignment operators:
  //  Vector& operator=(VectorView x);
  Vector& operator=(const Vector& v);
  Vector& operator=(const Array<Numeric>& v);
  Vector& operator=(Numeric x);

  // Resize function:
  void resize(Index n);

  // Destructor:
  ~Vector();
};

// Declare class Matrix:
class Matrix;


/** A constant view of a Matrix.

    This, together with the derived class MatrixView, contains the
    main implementation of a Matrix. It defines the concepts of
    MatrixView. Plus additionally the recursive subrange operator,
    which makes it possible to create a MatrixView from a subrange of
    a MatrixView.

    The class Matrix is just a special case of a MatrixView
    which also allocates storage. */
class ConstMatrixView {
public:
  // Member functions:
  Index nrows() const;
  Index ncols() const;

  // Const index operators:
  Numeric  operator()(Index r, Index c) const;
  ConstMatrixView operator()(const Range& r, const Range& c) const;
  ConstVectorView operator()(const Range& r, Index c) const;
  ConstVectorView operator()(Index r, const Range& c) const;

  // Functions returning iterators:
  ConstIterator2D begin() const;
  ConstIterator2D end() const;
  
  // Friends:
  friend class ConstVectorView;
  friend class MatrixView;
  friend class ConstIterator3D;
  friend class ConstTensor3View;
  friend ConstMatrixView transpose(ConstMatrixView m);

protected:
  // Constructors:
  ConstMatrixView();
  ConstMatrixView(Numeric *data, const Range& r, const Range& c);
  ConstMatrixView(Numeric *data,
                  const Range& pr, const Range& pc,
                  const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The row range of mdata that is actually used. */
  Range mrr;
  /** The column range of mdata that is actually used. */
  Range mcr;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The MatrixView class

    This contains the main implementation of a Matrix. It defines
    the concepts of MatrixView. Plus additionally the recursive
    subrange operator, which makes it possible to create a MatrixView
    from a subrange of a MatrixView. 

    The class Matrix is just a special case of a MatrixView
    which also allocates storage. */
class MatrixView : public ConstMatrixView {
public:

  // Const index operators:
  Numeric  operator()(Index r, Index c) const;
  ConstMatrixView operator()(const Range& r, const Range& c) const;
  ConstVectorView operator()(const Range& r, Index c) const;
  ConstVectorView operator()(Index r, const Range& c) const;
  // Index Operators:
  Numeric& operator()(Index r, Index c);
  MatrixView operator()(const Range& r, const Range& c);
  VectorView operator()(const Range& r, Index c);
  VectorView operator()(Index r, const Range& c);

  // Functions returning const iterators:
  ConstIterator2D begin() const;
  ConstIterator2D end() const;
  // Functions returning iterators:
  Iterator2D begin();
  Iterator2D end();
  
  // Assignment operators:
  MatrixView& operator=(const ConstMatrixView& v);
  MatrixView& operator=(const MatrixView& v);
  MatrixView& operator=(const Matrix& v);
  MatrixView& operator=(const ConstVectorView& v);
  MatrixView& operator=(Numeric x);

  // Other operators:
  MatrixView& operator*=(Numeric x);
  MatrixView& operator/=(Numeric x);
  MatrixView& operator+=(Numeric x);
  MatrixView& operator-=(Numeric x);

  MatrixView& operator*=(const ConstMatrixView& x);
  MatrixView& operator/=(const ConstMatrixView& x);
  MatrixView& operator+=(const ConstMatrixView& x);
  MatrixView& operator-=(const ConstMatrixView& x);

  MatrixView& operator*=(const ConstVectorView& x);
  MatrixView& operator/=(const ConstVectorView& x);
  MatrixView& operator+=(const ConstVectorView& x);
  MatrixView& operator-=(const ConstVectorView& x);

  // Friends:
  friend class VectorView;
  friend class Iterator3D;
  friend class Tensor3View;
  friend ConstMatrixView transpose(ConstMatrixView m);
  friend MatrixView transpose(MatrixView m);

protected:
  // Constructors:
  MatrixView();
  MatrixView(Numeric *data, const Range& r, const Range& c);
  MatrixView(Numeric *data,
             const Range& pr, const Range& pc,
             const Range& nr, const Range& nc);
};

/** The Matrix class. This is a MatrixView that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from MatrixView. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operator from scalar.
    3. Resize function. */
class Matrix : public MatrixView {
public:
  // Constructors:
  Matrix();
  Matrix(Index r, Index c);
  Matrix(Index r, Index c, Numeric fill);
  Matrix(const ConstMatrixView& v);
  Matrix(const Matrix& v);

  // Assignment operators:
  Matrix& operator=(const Matrix& x);
  Matrix& operator=(Numeric x);
  Matrix& operator=(const ConstVectorView& v);

  // Resize function:
  void resize(Index r, Index c);

  // Destructor:
  ~Matrix();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator1D origin,
                 const ConstIterator1D& end,
                 Iterator1D target);

inline void copy(Numeric x,
                 Iterator1D target,
                 const Iterator1D& end);

inline void copy(ConstIterator2D origin,
                 const ConstIterator2D& end,
                 Iterator2D target);

inline void copy(Numeric x,
                 Iterator2D target,
                 const Iterator2D& end);



// Declare the existance of class Array:
template<class base>
class Array;

/** An array of vectors. */
typedef Array<Vector> ArrayOfVector;

/** An array of matrices. */
typedef Array<Matrix> ArrayOfMatrix;


// Functions for Range:
// --------------------

/* Explicit constructor. 

  \param Start must be >= 0.

  \param Extent also. Although internally negative extent means "to the end",
  this can not be created this way, only with the joker. Zero
  extent is allowed, though, which corresponds to an empty range.

  \param Stride can be anything. It can be omitted, in which case the
  default value is 1. */
inline Range::Range(Index start, Index extent, Index stride) :
  mstart(start), mextent(extent), mstride(stride)
{
  // Start must be >= 0:
  assert( 0<=mstart );
  // Extent also. Although internally negative extent means "to the end",
  // this can not be created this way, only with the joker. Zero
  // extent is allowed, though, which corresponds to an empty range.
  assert( 0<=mextent );
  // Stride can be anything except 0.
  // SAB 2001-09-21: Allow 0 stride.
  //  assert( 0!=mstride);
}

/** Constructor with joker extent. Depending on the sign of stride,
    this means "to the end", or "to the beginning". */
inline Range::Range(Index start, Joker joker, Index stride) :
  mstart(start), mextent(-1), mstride(stride)
{
  // Start must be >= 0:
  assert( 0<=mstart );
}

/** Constructor with just a joker. This means, take everything. You
    can still optionally give a stride, though. This constructor is
    just shorter notation for Range(0,joker) */
inline Range::Range(Joker joker, Index stride) :
  mstart(0), mextent(-1), mstride(stride)
{
  // Nothing to do here.
}

/** Constructor which converts a range with joker to an explicit
    range.

    \param max_size The maximum allowed size of the vector. 
    \param r The new range, with joker. */
inline Range::Range(Index max_size, const Range& r) :
  mstart(r.mstart),
  mextent(r.mextent),
  mstride(r.mstride)
{
  // Start must be >= 0:
  assert( 0<=mstart );
  // ... and < max_size:
  assert( mstart<max_size );

  // Stride must be != 0:
  assert( 0!=mstride);
  
  // Convert negative extent (joker) to explicit extent
  if ( mextent<0 )
    {
      if ( 0<mstride )
        mextent = 1 + (max_size-1-mstart)/mstride;
      else
        mextent = 1 + (0-mstart)/mstride;
    }
  else
    {
#ifndef NDEBUG
      // Check that extent is ok:
      Index fin = mstart+(mextent-1)*mstride;
      assert( 0   <= fin );
      assert( fin <  max_size );
#endif
    }
}

/** Constructor of a new range relative to an old range. The new range
    may contain -1 for the stride, which acts as a joker.

    \param p Previous range.
    \param n New range. */
inline Range::Range(const Range& p, const Range& n) :
  mstart(p.mstart + n.mstart*p.mstride),
  mextent(n.mextent),
  mstride(p.mstride*n.mstride)
{
  // We have to juggle here a bit with previous, new, and resulting
  // quantities. I.e.;
  // p.mstride: Previous stride
  // n.mstride: New stride (as specified)
  // mstride:   Resulting stride (old*new)

  // Get the previous final element:
  Index prev_fin = p.mstart+(p.mextent-1)*p.mstride;

  // Resulting start must be >= previous start:
  assert( p.mstart<=mstart );
  // and <= prev_fin:
  assert( mstart<=prev_fin );

  // Resulting stride must be != 0:
  assert( 0!=mstride);

  // Convert negative extent (joker) to explicit extent
  if ( mextent<0 )
    {
      if ( 0<mstride )
        mextent = 1 + (prev_fin-mstart)/mstride;
      else
        mextent = 1 + (p.mstart-mstart)/mstride;
    }
  else
    {
#ifndef NDEBUG
      // Check that extent is ok:
      Index fin = mstart+(mextent-1)*mstride;
      assert( p.mstart <= fin      );
      assert( fin      <= prev_fin );
#endif
    }

}


// Functions for Iterator1D
// ------------------------

/** Default constructor. */
inline Iterator1D::Iterator1D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator1D::Iterator1D(const Iterator1D& o) :
  mx(o.mx), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator1D::Iterator1D(Numeric *x, Index stride) :
  mx(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline Iterator1D& Iterator1D::operator++()
{
  mx += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool Iterator1D::operator!=(const Iterator1D& other) const
{
  if (mx!=other.mx)
    return true;
  else
    return false;
}

/** Dereferencing. */
inline Numeric& Iterator1D::operator*() const
{
  return *mx;
}

// Functions for ConstIterator1D
// -----------------------------

/** Default constructor. */
inline ConstIterator1D::ConstIterator1D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator1D::ConstIterator1D(const ConstIterator1D& o) :
  mx(o.mx), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator1D::ConstIterator1D(Numeric *x, Index stride) :
  mx(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline ConstIterator1D& ConstIterator1D::operator++()
{
  mx += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool ConstIterator1D::operator!=(const ConstIterator1D& other) const
{
  if (mx!=other.mx)
    return true;
  else
    return false;
}

/** Dereferencing. */
inline const Numeric& ConstIterator1D::operator*() const
{
  return *mx;
}

// Functions for Iterator2D
// ------------------------

/** Default constructor. */
inline Iterator2D::Iterator2D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator2D::Iterator2D(const Iterator2D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator2D::Iterator2D(const VectorView& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline Iterator2D& Iterator2D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool Iterator2D::operator!=(const Iterator2D& other) const
{
  if ( msv.mdata + msv.mrange.mstart !=
       other.msv.mdata + other.msv.mrange.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline VectorView* const Iterator2D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline VectorView& Iterator2D::operator*()
{
  return msv;
}

// Functions for ConstIterator2D
// -----------------------------

/** Default constructor. */
inline ConstIterator2D::ConstIterator2D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator2D::ConstIterator2D(const ConstIterator2D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator2D::ConstIterator2D(const ConstVectorView& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline ConstIterator2D& ConstIterator2D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool ConstIterator2D::operator!=(const ConstIterator2D& other) const
{
  if ( msv.mdata + msv.mrange.mstart !=
       other.msv.mdata + other.msv.mrange.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline const ConstVectorView* ConstIterator2D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstVectorView& ConstIterator2D::operator*() const
{
  return msv;
}


// Functions for ConstVectorView:
// ------------------------------

/** Returns the number of elements.  The names `size' and `length'
    are already used by STL functions returning size_t. To avoid
    confusion we choose the name `nelem'. This is also more
    consistent with `nrow' and `ncol' for matrices.
    
    The value range of long, which is used to store the index is on a
    PC from -2147483648 to 2147483647. This means that a 15GB large
    array of float can be addressed with this index. So the extra bit
    that size_t has compared to long is not needed. */
inline Index ConstVectorView::nelem() const
{
  return mrange.mextent;
}

/** The sum of all elements of a Vector. */
inline Numeric ConstVectorView::sum() const
{
  Numeric s=0;
  ConstIterator1D i = begin();
  const ConstIterator1D e = end();

  for ( ; i!=e; ++i )
    s += *i;

  return s;
}

/** Plain const index operator. */
inline Numeric ConstVectorView::operator[](Index n) const
{
  // Check if index is valid:
  assert( 0<=n );
  assert( n<mrange.mextent );
  return *( mdata +
            mrange.mstart +
            n*mrange.mstride );
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Vector. This allows
    correct recursive behavior.  */
inline ConstVectorView ConstVectorView::operator[](const Range& r) const
{
  return ConstVectorView(mdata, mrange, r);
}

/** Return const iterator to first element. */
inline ConstIterator1D ConstVectorView::begin() const
{
  return ConstIterator1D(mdata+mrange.mstart, mrange.mstride);
}

/** Return const iterator behind last element. */
inline ConstIterator1D ConstVectorView::end() const
{
  return ConstIterator1D( mdata +
                          mrange.mstart +
                          (mrange.mextent)*mrange.mstride,
                          mrange.mstride );
}

/** Conversion to const 1 column matrix. */
inline ConstVectorView::operator ConstMatrixView() const
{
  return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Vector. */
inline ConstVectorView::ConstVectorView() :
  mrange(0,0), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
inline ConstVectorView::ConstVectorView(Numeric *data,
                                        const Range& range) :
  mrange(range),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub ranges from
    sub ranges. That means that the new range has to be interpreted
    relative to the original range. The new range may contain -1 for
    the extent which acts as a joker. However, the used Range
    constructor converts this to an explicit range, consistent with
    the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
inline ConstVectorView::ConstVectorView(Numeric *data,
                                        const Range& p,
                                        const Range& n) :
  mrange(p,n),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the vector. The iterators know which part of the vector
    is `active', and also the stride. */
inline std::ostream& operator<<(std::ostream& os, const ConstVectorView& v)
{
  ConstIterator1D i=v.begin();
  const ConstIterator1D end=v.end();

  if ( i!=end )
    {
      os << setw(3) << *i;
      ++i;
    }
  for ( ; i!=end ; ++i )
    {
      os << "\n" << setw(3) << *i;
    }

  return os;
}


// Functions for VectorView:
// ------------------------

/** Plain const index operator. Has to be redifined here, because the
    one from ConstVectorView is hidden. */
inline Numeric VectorView::operator[](Index n) const
{
  return ConstVectorView::operator[](n);
}

/** Const index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  Has to be redifined here, because the
    one from ConstVectorView is hidden. */
inline ConstVectorView VectorView::operator[](const Range& r) const
{
  return ConstVectorView::operator[](r);
}

/** Plain Index operator. */
inline Numeric& VectorView::operator[](Index n)
{
  // Check if index is valid:
  assert( 0<=n );
  assert( n<mrange.mextent );
  return *( mdata +
            mrange.mstart +
            n*mrange.mstride );
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  */
inline VectorView VectorView::operator[](const Range& r)
{
  return VectorView(mdata, mrange, r);
}

/** Return const iterator to first element. Has to be redefined here,
    since it is hiden by the non-const operator of the derived
    class.*/
inline ConstIterator1D VectorView::begin() const
{
  return ConstVectorView::begin();
}

/** Return const iterator behind last element. Has to be redefined
    here, since it is hiden by the non-const operator of the derived
    class.*/
inline ConstIterator1D VectorView::end() const
{
  return ConstVectorView::end();
}

/** Return iterator to first element. */
inline Iterator1D VectorView::begin()
{
  return Iterator1D(mdata+mrange.mstart, mrange.mstride);
}

/** Return iterator behind last element. */
inline Iterator1D VectorView::end()
{
  return Iterator1D( mdata +
                     mrange.mstart +
                     (mrange.mextent)*mrange.mstride,
                     mrange.mstride );
}

/** Assignment operator. This copies the data from another VectorView
    to this VectorView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this VectorView by
    setting its range. */
inline VectorView VectorView::operator=(const ConstVectorView& v)
{
  //  cout << "Assigning VectorView from ConstVectorView.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);
  
  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assignment from VectorView to VectorView. This is a tricky
    one. The problem is that since VectorView is derived from
    ConstVectorView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline VectorView VectorView::operator=(const VectorView& v)
{
  //  cout << "Assigning VectorView from VectorView.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);
  
  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assignment from Vector. This is important to avoid a bug when
    assigning a Vector to a VectorView. */
inline VectorView VectorView::operator=(const Vector& v)
{
  //  cout << "Assigning VectorView from Vector.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);
  
  copy( v.begin(), v.end(), begin() );

  return *this;
}

/** Assigning a scalar to a VectorView will set all elements to this
    value. */
inline VectorView VectorView::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
inline VectorView VectorView::operator*=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i *= x;
  return *this;
}

/** Division by scalar. */
inline VectorView VectorView::operator/=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i /= x;
  return *this;
}

/** Addition of scalar. */
inline VectorView VectorView::operator+=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i += x;
  return *this;
}

/** Subtraction of scalar. */
inline VectorView VectorView::operator-=(Numeric x)
{
  const Iterator1D e=end();
  for ( Iterator1D i=begin(); i!=e ; ++i )
    *i -= x;
  return *this;
}

/** Element-vise multiplication by another vector. */
inline VectorView VectorView::operator*=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i *= *s;
  return *this;
}

/** Element-vise division by another vector. */
inline VectorView VectorView::operator/=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i /= *s;
  return *this;
}

/** Element-vise addition of another vector. */
inline VectorView VectorView::operator+=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i += *s;
  return *this;
}

/** Element-vise subtraction of another vector. */
inline VectorView VectorView::operator-=(const ConstVectorView& x)
{
  assert( nelem()==x.nelem() );

  ConstIterator1D s=x.begin();

  Iterator1D i=begin();
  const Iterator1D e=end();

  for ( ; i!=e ; ++i,++s )
    *i -= *s;
  return *this;
}

/** Conversion to 1 column matrix. */
inline VectorView::operator MatrixView()
{
  return MatrixView(mdata,mrange,Range(mrange.mstart,1));
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Vector. */
inline VectorView::VectorView() :
  ConstVectorView()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
inline VectorView::VectorView(Numeric *data,
                            const Range& range) :
  ConstVectorView(data,range)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub ranges from
    sub ranges. That means that the new range has to be interpreted
    relative to the original range. The new range may contain -1 for
    the extent which acts as a joker. However, the used Range
    constructor converts this to an explicit range, consistent with
    the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
inline VectorView::VectorView(Numeric *data,
                              const Range& p,
                              const Range& n) :
  ConstVectorView(data,p,n)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors. */
inline void copy(ConstIterator1D origin,
                 const ConstIterator1D& end,
                 Iterator1D target)
{
  for ( ; origin!=end ; ++origin,++target )
    *target = *origin;
}

/** Copy a scalar to all elements. */
inline void copy(Numeric x,
                 Iterator1D target,
                 const Iterator1D& end)
{
  for ( ; target!=end ; ++target )
    *target = x;
}


// Functions for Vector:
// ---------------------

/** Default constructor. */
inline Vector::Vector() 
{
  // Nothing to do here
}

/** Constructor setting size. */
inline Vector::Vector(Index n) :
  VectorView( new Numeric[n],
             Range(0,n))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Vector::Vector(Index n, Numeric fill) :
  VectorView( new Numeric[n],
             Range(0,n))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x=mdata; x<mdata+n; ++x )
    *x = fill;
}

/** Constructor filling with values. 

    Examples:

    Vector v(1,5,1);  // 1, 2, 3, 4, 5
    Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
    Vector v(5,5,-1); // 5, 4, 3, 2, 1
*/
inline Vector::Vector(Numeric start, Index extent, Numeric stride) :
  VectorView( new Numeric[extent],
             Range(0,extent))
{
  // Fill with values:
  Numeric x = start;
  Iterator1D       i=begin();
  const Iterator1D e=end();
  for ( ; i!=e; ++i )
    {
      *i = x;
      x += stride;
    }
}

/** Copy constructor from VectorView. This automatically sets the size
    and copies the data. The vector created will have start zero and
    stride 1, independent on how these parameters are set for the
    original. So, what is copied is the data, not the shape
    of the selection. */
inline Vector::Vector(const ConstVectorView& v) :
  VectorView( new Numeric[v.nelem()],
              Range(0,v.nelem()))
{
  copy(v.begin(),v.end(),begin());
}

/** Copy constructor from Vector. This is important to override the
    automatically generated shallow constructor. We want deep copies!  */
inline Vector::Vector(const Vector& v) :
  VectorView( new Numeric[v.nelem()],
              Range(0,v.nelem()))
{
  copy(v.begin(),v.end(),begin());
}

/** Assignment from another Vector. Important to avoid segmentation
    fault for 
    x = Vector(n);
 */
inline Vector& Vector::operator=(const Vector& v)
{
  //  cout << "Assigning VectorView from Vector View.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.mrange.mextent);
  copy( v.begin(), v.end(), begin() );
  return *this;
}

/** Assignment operator from Array<Numeric>. This copies the data from
    a Array<Numeric> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Numeric> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy.

    Assignment operators are not inherited, so we have to make an
    explicit call to:

    VectorView VectorView::operator=(const Array<Numeric>& v)

    here. */  
inline Vector& Vector::operator=(const Array<Numeric>& x)
{
  VectorView::operator=(x);
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */  
inline Vector& Vector::operator=(Numeric x)
{
  VectorView::operator=(x);
  return *this;
}

// /** Assignment operator from VectorView. Assignment operators are not
//     inherited. */  
// inline Vector& Vector::operator=(const VectorView v)
// {
//   cout << "Assigning Vector from Vector View.\n";
//  // Check that sizes are compatible:
//   assert(mrange.mextent==v.mrange.mextent);
//   VectorView::operator=(v);
//   return *this;
// }

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Vector is not
    initialized, so it will contain random values.  */
inline void Vector::resize(Index n)
{
  if ( mrange.mextent != n )
    {
      delete mdata;
      mdata = new Numeric[n];
      mrange.mstart = 0;
      mrange.mextent = n;
      mrange.mstride = 1;
    }
}

/** Destructor for Vector. This is important, since Vector uses new to
    allocate storage. */
inline Vector::~Vector()
{
  delete [] mdata;
}


// Functions for ConstMatrixView:
// ------------------------------

/** Returns the number of rows. */
inline Index ConstMatrixView::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstMatrixView::ncols() const
{
  return mcr.mextent;
}

/** Plain const index operator. */
inline Numeric ConstMatrixView::operator()(Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  return *( mdata +
            mrr.mstart +
            r*mrr.mstride +
            mcr.mstart +
            c*mcr.mstride );
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior.  */
inline ConstMatrixView ConstMatrixView::operator()(const Range& r,
                                                   const Range& c) const
{
  return ConstMatrixView(mdata, mrr, mcr, r, c);
}

/** Const index operator returning a column as an object of type
    ConstVectorView.

    \param r A range of rows.
    \param c Index of selected column */
inline ConstVectorView ConstMatrixView::operator()(const Range& r, Index c) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata + mcr.mstart + c*mcr.mstride,
                         mrr, r);
}

/** Const index operator returning a row as an object of type
    ConstVectorView.

    \param r Index of selected row.
    \param c Range of columns */
inline ConstVectorView ConstMatrixView::operator()(Index r, const Range& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return ConstVectorView(mdata + mrr.mstart + r*mrr.mstride,
                         mcr, c);
}

/** Return const iterator to first row. */
inline ConstIterator2D ConstMatrixView::begin() const
{
  return ConstIterator2D(ConstVectorView(mdata+mrr.mstart,
                                         mcr),
                         mrr.mstride);
}

/** Return const iterator behind last row. */
inline ConstIterator2D ConstMatrixView::end() const
{
  return ConstIterator2D( ConstVectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
                                          mcr),
                          mrr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstMatrixView::ConstMatrixView() :
  mrr(0,0,1), mcr(0,0,1), mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
inline ConstMatrixView::ConstMatrixView(Numeric *data,
                                        const Range& rr,
                                        const Range& cr) :
  mrr(rr),
  mcr(cr),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
inline ConstMatrixView::ConstMatrixView(Numeric *data,
                                        const Range& pr, const Range& pc,
                                        const Range& nr, const Range& nc) :
  mrr(pr,nr),
  mcr(pc,nc),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the matrix. The iterators know which part of the matrix
    is `active', and also the strides in both directions. This
    function is a bit more complicated than necessary to illustrate
    the concept, because the formating should look nice. This means
    that the first row, and the first element in each row, have to be
    treated individually. */
inline std::ostream& operator<<(std::ostream& os, const ConstMatrixView& v)
{
  // Row iterators:
  ConstIterator2D ir=v.begin();
  const ConstIterator2D end_row=v.end();

  if ( ir!=end_row )
    {
      ConstIterator1D ic =  ir->begin();
      const ConstIterator1D end_col = ir->end();

      if ( ic!=end_col )
        {
          os << setw(3) << *ic;
          ++ic;
        }
      for ( ; ic!=end_col ; ++ic )
        {
          os << " " << setw(3) << *ic;
        }
      ++ir;
    }
  for ( ; ir!=end_row ; ++ir )
    {
      ConstIterator1D ic =  ir->begin();
      const ConstIterator1D end_col = ir->end();

      os << "\n";
      if ( ic!=end_col )
        {
          os << setw(3) << *ic;
          ++ic;
        }
      for ( ; ic!=end_col ; ++ic )
        {
          os << " " << setw(3) << *ic;
        }
    }

  return os;
}


// Functions for MatrixView:
// -------------------------

/** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline Numeric MatrixView::operator()(Index r, Index c) const
{
  return ConstMatrixView::operator()(r,c);
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline ConstMatrixView MatrixView::operator()(const Range& r, const Range& c) const
{
  return ConstMatrixView::operator()(r,c);  
}

/** Const index operator returning a column as an object of type
    ConstVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r A range of rows.
    \param c Index of selected column */
inline ConstVectorView MatrixView::operator()(const Range& r, Index c) const
{
  return ConstMatrixView::operator()(r,c);
}

/** Const index operator returning a row as an object of type
    ConstVectorView. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.

    \param r Index of selected row.
    \param c Range of columns */
inline ConstVectorView MatrixView::operator()(Index r, const Range& c) const
{
  return ConstMatrixView::operator()(r,c);
}

/** Plain index operator. */
inline Numeric& MatrixView::operator()(Index r, Index c)
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<mrr.mextent );
  assert( c<mcr.mextent );

  return *( mdata +
            mrr.mstart +
            r*mrr.mstride +
            mcr.mstart +
            c*mcr.mstride );
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Matrix. This allows correct
    recursive behavior.  */
inline MatrixView MatrixView::operator()(const Range& r, const Range& c)
{
  return MatrixView(mdata, mrr, mcr, r, c);
}

/** Index operator returning a column as an object of type VectorView.

    \param r A range of rows.
    \param c Index of selected column */
inline VectorView MatrixView::operator()(const Range& r, Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return VectorView(mdata + mcr.mstart + c*mcr.mstride,
                    mrr, r);
}

/** Index operator returning a row as an object of type VectorView.

    \param r Index of selected row.
    \param c Range of columns */
inline VectorView MatrixView::operator()(Index r, const Range& c)
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r <  mrr.mextent );

  return VectorView(mdata + mrr.mstart + r*mrr.mstride,
                    mcr, c);
}

/** Return const iterator to first row. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator2D MatrixView::begin() const
{
  return ConstMatrixView::begin();
}

/** Return const iterator behind last row. */
inline ConstIterator2D MatrixView::end() const
{
  return ConstMatrixView::end();
}

/** Return iterator to first row. */
inline Iterator2D MatrixView::begin()
{
  return Iterator2D(VectorView(mdata+mrr.mstart, mcr),
                    mrr.mstride);
}

/** Return iterator behind last row. */
inline Iterator2D MatrixView::end()
{
  return Iterator2D( VectorView(mdata + mrr.mstart + (mrr.mextent)*mrr.mstride,
                                mcr),
                     mrr.mstride );
}

/** Assignment operator. This copies the data from another MatrixView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
inline MatrixView& MatrixView::operator=(const ConstMatrixView& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from MatrixView to MatrixView. This is a tricky
    one. The problem is that since MatrixView is derived from
    ConstMatrixView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline MatrixView& MatrixView::operator=(const MatrixView& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Matrix. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline MatrixView& MatrixView::operator=(const Matrix& m)
{
  // Check that sizes are compatible:
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a vector. This copies the data from a VectorView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
inline MatrixView& MatrixView::operator=(const ConstVectorView& v)
{
  // Check that sizes are compatible:
  assert( mrr.mextent==v.nelem() );
  assert( mcr.mextent==1         );
  //  dummy = ConstMatrixView(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstMatrixView dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Assigning a scalar to a MatrixView will set all elements to this
    value. */
inline MatrixView& MatrixView::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Multiplication by scalar. */
inline MatrixView& MatrixView::operator*=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c *= x;
    }
  return *this;
}

/** Division by scalar. */
inline MatrixView& MatrixView::operator/=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c /= x;
    }
  return *this;
}

/** Addition of scalar. */
inline MatrixView& MatrixView::operator+=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c += x;
    }
  return *this;
}

/** Subtraction of scalar. */
inline MatrixView& MatrixView::operator-=(Numeric x)
{
  const Iterator2D er=end();
  for ( Iterator2D r=begin(); r!=er ; ++r )
    {
      const Iterator1D ec = r->end();
      for ( Iterator1D c = r->begin(); c!=ec ; ++c )
        *c -= x;
    }
  return *this;
}

/** Element-vise multiplication by another Matrix. */
inline MatrixView& MatrixView::operator*=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c *= *sc;
    }
  return *this;
}

/** Element-vise division by another Matrix. */
inline MatrixView& MatrixView::operator/=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c /= *sc;
    }
  return *this;
}

/** Element-vise addition of another Matrix. */
inline MatrixView& MatrixView::operator+=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c += *sc;
    }
  return *this;
}

/** Element-vise subtraction of another Matrix. */
inline MatrixView& MatrixView::operator-=(const ConstMatrixView& x)
{
  assert(nrows()==x.nrows());
  assert(ncols()==x.ncols());
  ConstIterator2D  sr = x.begin();
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sr )
    {
      ConstIterator1D  sc = sr->begin(); 
      Iterator1D        c = r->begin();
      const Iterator1D ec = r->end();
      for ( ; c!=ec ; ++c,++sc )
        *c -= *sc;
    }
  return *this;
}

/** Element-vise multiplication by a Vector (acting like a 1-column Matrix). */
inline MatrixView& MatrixView::operator*=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin(); 
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c *= *sc;
    }
  return *this;
}

/** Element-vise division by a Vector (acting like a 1-column Matrix). */
inline MatrixView& MatrixView::operator/=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin(); 
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c /= *sc;
    }
  return *this;
}

/** Element-vise addition of a Vector (acting like a 1-column Matrix). */
inline MatrixView& MatrixView::operator+=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin(); 
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c += *sc;
    }
  return *this;
}

/** Element-vise subtraction of a Vector (acting like a 1-column Matrix). */
inline MatrixView& MatrixView::operator-=(const ConstVectorView& x)
{
  assert(nrows()==x.nelem());
  assert(ncols()==1);
  ConstIterator1D  sc = x.begin(); 
  Iterator2D        r = begin();
  const Iterator2D er = end();
  for ( ; r!=er ; ++r,++sc )
    {
      Iterator1D        c = r->begin();
      *c -= *sc;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Matrix. */
inline MatrixView::MatrixView() :
  ConstMatrixView()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
inline MatrixView::MatrixView(Numeric *data,
                            const Range& rr,
                            const Range& cr) :
  ConstMatrixView(data, rr, cr)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
inline MatrixView::MatrixView(Numeric *data,
                            const Range& pr, const Range& pc,
                            const Range& nr, const Range& nc) :
  ConstMatrixView(data,pr,pc,nr,nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors.

    Origin, end, and target are 2D iterators, marking rows in a
    matrix. For each row the 1D iterator is obtained and used to copy
    the elements. 
*/
inline void copy(ConstIterator2D origin,
                 const ConstIterator2D& end,
                 Iterator2D target)
{
  for ( ; origin!=end ; ++origin,++target )
    {
      ConstIterator1D       o = origin->begin();
      const ConstIterator1D e = origin->end();
      Iterator1D            t = target->begin();
      for ( ; o!=e ; ++o,++t )
        *t = *o;
    }
}

/** Copy a scalar to all elements. */
inline void copy(Numeric x,
                 Iterator2D target,
                 const Iterator2D& end)
{
  for ( ; target!=end ; ++target )
    {
      Iterator1D       t = target->begin();
      const Iterator1D e = target->end();
      for ( ; t!=e ; ++t )
        *t = x;
    }
}


// Functions for Matrix:
// ---------------------

/** Default constructor. */
inline Matrix::Matrix() :
  MatrixView::MatrixView()
{
  // Nothing to do here. However, note that the default constructor
  // for MatrixView has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized. 
}

/** Constructor setting size. This constructor has to set the stride
    in the row range correctly! */
inline Matrix::Matrix(Index r, Index c) :
  MatrixView( new Numeric[r*c],
             Range(0,r,c),
             Range(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Matrix::Matrix(Index r, Index c, Numeric fill) :
  MatrixView( new Numeric[r*c],
              Range(0,r,c),
              Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x=mdata; x<mdata+r*c; ++x )
    *x = fill;
}

/** Copy constructor from MatrixView. This automatically sets the size
    and copies the data. */
inline Matrix::Matrix(const ConstMatrixView& m) :
  MatrixView( new Numeric[m.nrows()*m.ncols()],
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Matrix. This automatically sets the size
    and copies the data. */
inline Matrix::Matrix(const Matrix& m) :
  MatrixView( new Numeric[m.nrows()*m.ncols()],
             Range( 0, m.nrows(), m.ncols() ),
             Range( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty matrix, then it will have
  // 0 colunns. But this is used to initialize the stride of the row
  // Range! Thus, this method has to be consistent with the behaviour
  // of Range::Range. For now, Range::Range allows also stride 0.
  copy(m.begin(),m.end(),begin());
}

/** Assignment operator from another matrix. It is important that this
    operator exists. Otherwise the = operator seems to copy references
    instead of content in some cases. 

    The Behavior of this one is a bit special: If the size of the
    target Matrix is 0 then it will be automatically resized to match
    (this is needed to have the correct initialization for constructed
    classes that use the assignment operator to initialize their
    data). 
*/
inline Matrix& Matrix::operator=(const Matrix& m)
{
  //  cout << "Matrix copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "             n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid matrix, so we just
  // have to check one.
  if ( 0 == mrr.mextent )
    {
      // Adjust if previously empty.
      resize( m.mrr.mextent, m.mcr.mextent ); 
    }
  else
    {
      // Check that sizes are compatible:
      assert( mrr.mextent==m.mrr.mextent );
      assert( mcr.mextent==m.mcr.mextent );
    }

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment operator from scalar. Assignment operators also seem to
    be not inherited. */
inline Matrix& Matrix::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Assignment from a vector. This copies the data from a VectorView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
inline Matrix& Matrix::operator=(const ConstVectorView& v)
{
  // Check that sizes are compatible:
  assert( mrr.mextent==v.nelem() );
  assert( mcr.mextent==1         );
  //  dummy = ConstMatrixView(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstMatrixView dummy(v);
  copy( dummy.begin(), dummy.end(), begin() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Matrix is not
    initialized, so it will contain random values.*/
inline void Matrix::resize(Index r, Index c)
{
  assert( 0<=r );
  assert( 0<=c );

  if ( mrr.mextent!=r || mcr.mextent!=c )
    {
      delete mdata;
      mdata = new Numeric[r*c];

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
}

/** Destructor for Matrix. This is important, since Matrix uses new to
    allocate storage. */
inline Matrix::~Matrix()
{
//   cout << "Destroying a Matrix:\n"
//        << *this << "\n........................................\n";
  delete [] mdata;
}


// Some general Matrix Vector functions:

/** Scalar product. The two vectors may be identical. */
inline Numeric operator*(const ConstVectorView& a, const ConstVectorView& b)
{
  // Check dimensions:
  assert( a.nelem() == b.nelem() );

  const ConstIterator1D ae = a.end();
  ConstIterator1D       ai = a.begin();
  ConstIterator1D       bi = b.begin();

  Numeric res = 0;
  for ( ; ai!=ae ; ++ai, ++bi )
    res += (*ai) * (*bi);

  return res;
}

/** Matrix Vector multiplication. y = M*x. Note that the order is different
    from MTL, output comes first! Dimensions of y, M, and x must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping MatrixViews belonging
    to the same Matrix will lead to unpredictable results. */
inline void mult( VectorView y,
                  const ConstMatrixView& M,
                  const ConstVectorView& x )
{
  // Check dimensions:
  assert( y.nelem() == M.nrows() );
  assert( M.ncols() == x.nelem() );

  // Let's get a 1D iterator for y:
  const Iterator1D ye = y.end();
  Iterator1D       yi = y.begin();

  // ... a 2D iterator for M:
  ConstIterator2D  mi = M.begin();

  // ... and 1D iterators pointing to the start and end of x:
  const ConstIterator1D xs = x.begin();
  const ConstIterator1D xe = x.end();

  // This walks through the rows of y and M:
  for ( ; yi!=ye ; ++yi, ++mi )
    {
      // Compute the scalar product between this row of M and x: 
      Numeric dummy = 0;
      
      // Iterators for row of M:
      ConstIterator1D       ri = mi->begin();

      // ... and for x (always initialized to the start of x):
      ConstIterator1D       xi = xs;

      for ( ; xi!=xe ; ++xi, ++ri )
        {
          dummy += (*ri) * (*xi);
        }

      *yi = dummy;
    }
}

/** Matrix multiplication. A = B*C. Note that the order is different
    from MTL, output comes first! Dimensions of A, B, and C must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping MatrixViews belonging
    to the same Matrix will lead to unpredictable results. */
inline void mult( MatrixView A,
                  const ConstMatrixView& B,
                  const ConstMatrixView& C )
{
  // Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Let's get the transpose of C, so that we can use 2D iterators to
  // access the columns (= rows of the transpose).
  ConstMatrixView CT = transpose(C);

  const Iterator2D ae = A.end();
  Iterator2D       ai = A.begin();
  ConstIterator2D  bi = B.begin();

  // This walks through the rows of A and B:
  for ( ; ai!=ae ; ++ai, ++bi )
    {
      const Iterator1D ace = ai->end();
      Iterator1D       aci = ai->begin();
      ConstIterator2D  cti = CT.begin();

      // This walks through the columns of A with a 1D iterator, and
      // at the same time through the rows of CT, which are the columns of
      // C, with a 2D iterator:
      for ( ; aci!=ace ; ++aci, ++cti )
        {
          // The operator * is used to compute the scalar product
          // between rows of B and rows of C.transpose().
          *aci = (*bi) * (*cti);
        }
    }
}

/** Const version of transpose. */
inline ConstMatrixView transpose(ConstMatrixView m)
{
  return ConstMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
    transpose. The original is not changed! */
inline MatrixView transpose(MatrixView m)
{
  return MatrixView(m.mdata, m.mcr, m.mrr);
}

/** A generic transform function for vectors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    Although the matrix version of this can also be used for vectors,
    thanks to the automatic interpretation of a vector as a one column
    matrix, this one is slightly more efficient. However, the
    difference is very small (only a few percent). 

    The two views may be the same one, in which case the
    conversion happens in place. 

    \retval   y   the results of the function acting on each element of x
    \param    my_func a function (e.g., sqrt)
    \param    x   a vector */
inline void transform( VectorView y,
                       double (&my_func)(double),
                       ConstVectorView x )
{
  // Check dimensions:
  assert( y.nelem()==x.nelem() );

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();
  Iterator1D            yi = y.begin();
  for ( ; xi!=xe; ++xi, ++yi )
    *yi = my_func(*xi);
}

/** A generic transform function for matrices, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    This function can also be used for Vectors, because there is a
    conversion to MatrixView.

    The two Matrix views may be the same one, in which case the
    conversion happens in place. 

   \retval   y   the results of the function acting on each element of x
   \param    my_func a function (e.g., sqrt)
   \param    x   a matrix */
inline void transform( MatrixView y,
                       double (&my_func)(double),
                       ConstMatrixView x )
{
  // Check dimensions:
  assert( y.nrows()==x.nrows() );
  assert( y.ncols()==x.ncols() );

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();
  Iterator2D             ry = y.begin();
  for ( ; rx!=rxe; ++rx, ++ry )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();
      Iterator1D             cy = ry->begin();
      for ( ; cx!=cxe; ++cx, ++cy )
        *cy = my_func(*cx);
    }
}

/** Max function, vector version. */
inline Numeric max(const ConstVectorView& x)
{
  // Initial value for max:
  Numeric max = x[0];

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();

  for ( ; xi!=xe ; ++xi )
    {
      if ( *xi > max )
        max = *xi;
    }

  return max;
}

/** Max function, matrix version. */
inline Numeric max(const ConstMatrixView& x)
{
  // Initial value for max:
  Numeric max = x(0,0);

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();

  for ( ; rx!=rxe ; ++rx )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();

      for ( ; cx!=cxe ; ++cx )
        if ( *cx > max )
          max = *cx;
    }
  
  return max;
}

/** Min function, vector version. */
inline Numeric min(const ConstVectorView& x)
{
  // Initial value for min:
  Numeric min = x[0];

  const ConstIterator1D xe = x.end();
  ConstIterator1D       xi = x.begin();

  for ( ; xi!=xe ; ++xi )
    {
      if ( *xi < min )
        min = *xi;
    }

  return min;
}

/** Min function, matrix version. */
inline Numeric min(const ConstMatrixView& x)
{
  // Initial value for min:
  Numeric min = x(0,0);

  const ConstIterator2D rxe = x.end();
  ConstIterator2D        rx = x.begin();

  for ( ; rx!=rxe ; ++rx )
    {
      const ConstIterator1D cxe = rx->end();
      ConstIterator1D        cx = rx->begin();

      for ( ; cx!=cxe ; ++cx )
        if ( *cx < min )
          min = *cx;
    }
  
  return min;
}


#endif    // matpackI_h
