/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>
                      Wolfram-Andre Haas <wolhaas@hermes.fho-emden.de>

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
  Implementation of Tensors of Rank 5.

  Based on Tensor3 by Stefan Buehler.

  The five dimensions are called: shelf, book, page, row, column.

  \author Wolfram-Andre Haas
  \date   2002-03-01
 */

#ifndef matpackV_h
#define matpackV_h

#include <iomanip>
#include "matpackI.h"
#include "matpackIII.h"
#include "matpackIV.h"

/** The outermost iterator class for rank 5 tensors. This takes into
    account the defined strided. */
class Iterator5D {
public:
  // Constructors:
  Iterator5D();
  Iterator5D(const Iterator5D& o);
  Iterator5D(const Tensor4View& x, Index stride);

  // Operators:
  Iterator5D& operator++();
  bool operator!=(const Iterator5D& other) const;
  Tensor4View* const operator->();
  Tensor4View& operator*();

private:
  /** Current position. */
  Tensor4View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator5D. */
class ConstIterator5D {
public:
  // Constructors:
  ConstIterator5D();
  ConstIterator5D(const ConstIterator5D& o);
  ConstIterator5D(const ConstTensor4View& x, Index stride);

  // Operators:
  ConstIterator5D& operator++();
  bool operator!=(const ConstIterator5D& other) const;
  const ConstTensor4View* operator->() const;
  const ConstTensor4View& operator*()  const;

private:
  /** Current position. */
  ConstTensor4View msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor5:
class Tensor5;


/** A constant view of a Tensor5.

    This, together with the derived class Tensor5View, contains the
    main implementation of a Tensor5. It defines the concepts of
    Tensor5View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor5View from a subrange of
    a Tensor5View.

    The five dimensions of the tensor are called: shelf, book, page, row, column.

    The class Tensor5 is just a special case of a Tensor5View
    which also allocates storage. */
class ConstTensor5View {
public:
  // Member functions:
  Index nshelves() const;
  Index nbooks()   const;
  Index npages()   const;
  Index nrows()    const;
  Index ncols()    const;

  // Const index operators:
  ConstTensor5View operator()( const Range& s, const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor4View operator()( const Range& s, const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor4View operator()( const Range& s, const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor4View operator()( const Range& s, const Range& b, Index p,        const Range& r, const Range& c ) const;
  ConstTensor4View operator()( const Range& s, Index b,        const Range& p, const Range& r, const Range& c ) const;
  ConstTensor4View operator()( Index s,        const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& s, const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstTensor3View operator()( const Range& s, const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& s, const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& s, Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& s, Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& s, Index b,        Index p,        const Range& r, const Range& c ) const;
  ConstTensor3View operator()( Index s,        const Range& b, Index p,        const Range& r, const Range& c ) const;
  ConstTensor3View operator()( Index s,        const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( Index s,        const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( Index s,        Index b,        const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& s, const Range& b, Index p,        Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& s, Index b,        const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& s, Index b,        Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& s, Index b,        Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index s,        const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index s,        const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index s,        const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( Index s,        Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index s,        Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index s,        Index b,        Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( const Range& s, Index b,        Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index s,        const Range& b, Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index s,        Index b,        const Range& p, Index r,        Index c        ) const;
  ConstVectorView  operator()( Index s,        Index b,        Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( Index s,        Index b,        Index p,        Index r,        const Range& c ) const;

  Numeric          operator()( Index s,        Index b,        Index p,        Index r,        Index c        ) const;

  // Functions returning iterators:
  ConstIterator5D begin() const;
  ConstIterator5D end()   const;

  // Friends:
  friend class Tensor5View;

protected:
  // Constructors:
  ConstTensor5View();
  ConstTensor5View(Numeric *data,
                   const Range& s, const Range& b, const Range& p, const Range& r, const Range& c);
  ConstTensor5View(Numeric *data,
                   const Range& ps, const Range& pb, const Range& pp, const Range& pr, const Range& pc,
                   const Range& ns, const Range& nb, const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The shelf range of mdata that is actually used. */
  Range msr;
  /** The book range of mdata that is actually used. */
  Range mbr;
  /** The page range of mdata that is actually used. */
  Range mpr;
  /** The row range of mdata that is actually used. */
  Range mrr;
  /** The column range of mdata that is actually used. */
  Range mcr;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The Tensor5View class

    This contains the main implementation of a Tensor5. It defines
    the concepts of Tensor5View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor5View
    from a subrange of a Tensor5View.

    The class Tensor5 is just a special case of a Tensor5View
    which also allocates storage. */
class Tensor5View : public ConstTensor5View {
public:

  // Const index operators:
  ConstTensor5View operator()( const Range& s, const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor4View operator()( const Range& s, const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor4View operator()( const Range& s, const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor4View operator()( const Range& s, const Range& b, Index p,        const Range& r, const Range& c ) const;
  ConstTensor4View operator()( const Range& s, Index b,        const Range& p, const Range& r, const Range& c ) const;
  ConstTensor4View operator()( Index s,        const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& s, const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstTensor3View operator()( const Range& s, const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& s, const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& s, Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& s, Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& s, Index b,        Index p,        const Range& r, const Range& c ) const;
  ConstTensor3View operator()( Index s,        const Range& b, Index p,        const Range& r, const Range& c ) const;
  ConstTensor3View operator()( Index s,        const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( Index s,        const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( Index s,        Index b,        const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& s, const Range& b, Index p,        Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& s, Index b,        const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& s, Index b,        Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& s, Index b,        Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index s,        const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index s,        const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index s,        const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( Index s,        Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index s,        Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index s,        Index b,        Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( const Range& s, Index b,        Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index s,        const Range& b, Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index s,        Index b,        const Range& p, Index r,        Index c        ) const;
  ConstVectorView  operator()( Index s,        Index b,        Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( Index s,        Index b,        Index p,        Index r,        const Range& c ) const;

  Numeric          operator()( Index s,        Index b,        Index p,        Index r,        Index c        ) const;

  // Non-const index operators:

  Tensor5View operator()( const Range& s, const Range& b, const Range& p, const Range& r, const Range& c );

  Tensor4View operator()( const Range& s, const Range& b, const Range& p, const Range& r, Index c        );
  Tensor4View operator()( const Range& s, const Range& b, const Range& p, Index r,        const Range& c );
  Tensor4View operator()( const Range& s, const Range& b, Index p,        const Range& r, const Range& c );
  Tensor4View operator()( const Range& s, Index b,        const Range& p, const Range& r, const Range& c );
  Tensor4View operator()( Index s,        const Range& b, const Range& p, const Range& r, const Range& c );

  Tensor3View operator()( const Range& s, const Range& b, const Range& p, Index r,        Index c        );
  Tensor3View operator()( const Range& s, const Range& b, Index p,        const Range& r, Index c        );
  Tensor3View operator()( const Range& s, const Range& b, Index p,        Index r,        const Range& c );
  Tensor3View operator()( const Range& s, Index b,        const Range& p, Index r,        const Range& c );
  Tensor3View operator()( const Range& s, Index b,        const Range& p, const Range& r, Index c        );
  Tensor3View operator()( const Range& s, Index b,        Index p,        const Range& r, const Range& c );
  Tensor3View operator()( Index s,        const Range& b, Index p,        const Range& r, const Range& c );
  Tensor3View operator()( Index s,        const Range& b, const Range& p, Index r,        const Range& c );
  Tensor3View operator()( Index s,        const Range& b, const Range& p, const Range& r, Index c        );
  Tensor3View operator()( Index s,        Index b,        const Range& p, const Range& r, const Range& c );

  MatrixView  operator()( const Range& s, const Range& b, Index p,        Index r,        Index c        );
  MatrixView  operator()( const Range& s, Index b,        const Range& p, Index r,        Index c        );
  MatrixView  operator()( const Range& s, Index b,        Index p,        const Range& r, Index c        );
  MatrixView  operator()( const Range& s, Index b,        Index p,        Index r,        const Range& c );
  MatrixView  operator()( Index s,        const Range& b, Index p,        Index r,        const Range& c );
  MatrixView  operator()( Index s,        const Range& b, Index p,        const Range& r, Index c        );
  MatrixView  operator()( Index s,        const Range& b, const Range& p, Index r,        Index c        );
  MatrixView  operator()( Index s,        Index b,        const Range& p, const Range& r, Index c        );
  MatrixView  operator()( Index s,        Index b,        const Range& p, Index r,        const Range& c );
  MatrixView  operator()( Index s,        Index b,        Index p,        const Range& r, const Range& c );

  VectorView  operator()( const Range& s, Index b,        Index p,        Index r,        Index c        );
  VectorView  operator()( Index s,        const Range& b, Index p,        Index r,        Index c        );
  VectorView  operator()( Index s,        Index b,        const Range& p, Index r,        Index c        );
  VectorView  operator()( Index s,        Index b,        Index p,        const Range& r, Index c        );
  VectorView  operator()( Index s,        Index b,        Index p,        Index r,        const Range& c );

  Numeric&    operator()( Index s,        Index b,        Index p,        Index r,        Index c        );

  // Functions returning const iterators:
  ConstIterator5D begin() const;
  ConstIterator5D end()   const;
  // Functions returning iterators:
  Iterator5D begin();
  Iterator5D end();

  // Assignment operators:
  Tensor5View& operator=(const ConstTensor5View& v);
  Tensor5View& operator=(const Tensor5View& v);
  Tensor5View& operator=(const Tensor5& v);
  Tensor5View& operator=(Numeric x);

  // Other operators:
  Tensor5View& operator*=(Numeric x);
  Tensor5View& operator/=(Numeric x);
  Tensor5View& operator+=(Numeric x);
  Tensor5View& operator-=(Numeric x);

  Tensor5View& operator*=(const ConstTensor5View& x);
  Tensor5View& operator/=(const ConstTensor5View& x);
  Tensor5View& operator+=(const ConstTensor5View& x);
  Tensor5View& operator-=(const ConstTensor5View& x);

  // Friends:
  // friend class VectorView;
  // friend ConstTensor5View transpose(ConstTensor5View m);
  // friend Tensor5View transpose(Tensor5View m);

protected:
  // Constructors:
  Tensor5View();
  Tensor5View(Numeric *data,
              const Range& s, const Range& b, const Range& p, const Range& r, const Range& c);
  Tensor5View(Numeric *data,
              const Range& ps, const Range& pb, const Range& pp, const Range& pr, const Range& pc,
              const Range& ns, const Range& nb, const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor5 class. This is a Tensor5View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor5View. Additionally defined here
    are:

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor5 : public Tensor5View {
public:
  // Constructors:
  Tensor5();
  Tensor5(Index s, Index b, Index p, Index r, Index c);
  Tensor5(Index s, Index b, Index p, Index r, Index c, Numeric fill);
  Tensor5(const ConstTensor5View& v);
  Tensor5(const Tensor5& v);

  // Assignment operators:
  Tensor5& operator=(const Tensor5& x);
  Tensor5& operator=(Numeric x);

  // Resize function:
  void resize(Index s, Index b, Index p, Index r, Index c);

  // Destructor:
  ~Tensor5();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator5D origin,
                 const ConstIterator5D& end,
                 Iterator5D target);

inline void copy(Numeric x,
                 Iterator5D target,
                 const Iterator5D& end);


// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor5. */
typedef Array<Tensor5> ArrayOfTensor5;


// Functions for Iterator5D
// ------------------------

/** Default constructor. */
inline Iterator5D::Iterator5D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator5D::Iterator5D(const Iterator5D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator5D::Iterator5D(const Tensor4View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here.
}

/** Prefix increment operator. */
inline Iterator5D& Iterator5D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool Iterator5D::operator!=(const Iterator5D& other) const
{
  if ( msv.mdata +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mbr.mstart +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 4D iterators. */
inline Tensor4View* const Iterator5D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline Tensor4View& Iterator5D::operator*()
{
  return msv;
}

// Functions for ConstIterator5D
// -----------------------------

/** Default constructor. */
inline ConstIterator5D::ConstIterator5D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator5D::ConstIterator5D(const ConstIterator5D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator5D::ConstIterator5D(const ConstTensor4View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here.
}

/** Prefix increment operator. */
inline ConstIterator5D& ConstIterator5D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool ConstIterator5D::operator!=(const ConstIterator5D& other) const
{
  if ( msv.mdata +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mbr.mstart +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 4D iterators. */
inline const ConstTensor4View* ConstIterator5D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstTensor4View& ConstIterator5D::operator*() const
{
  return msv;
}

// Functions for ConstTensor5View:
// ------------------------------

/** Returns the number of shelfs. */
inline Index ConstTensor5View::nshelves() const
{
  return msr.mextent;
}

/** Returns the number of books. */
inline Index ConstTensor5View::nbooks() const
{
  return mbr.mextent;
}

/** Returns the number of pages. */
inline Index ConstTensor5View::npages() const
{
  return mpr.mextent;
}

/** Returns the number of rows. */
inline Index ConstTensor5View::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor5View::ncols() const
{
  return mcr.mextent;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor5. This allows
    correct recursive behavior.  */
inline ConstTensor5View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  return ConstTensor5View(mdata,
                          msr, mbr, mpr, mrr, mcr,
                          s,   b,   p,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
inline ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     const Range& p,
                                                     const Range& r,
                                                     Index c       ) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ConstTensor4View(mdata +
                          mcr.mstart + c * mcr.mstride,
                          msr, mbr, mpr, mrr,
                          s,   b,   p,   r);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
inline ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     const Range& p,
                                                     Index r,
                                                     const Range& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return ConstTensor4View(mdata +
                          mrr.mstart + r * mrr.mstride,
                          msr, mbr, mpr, mcr,
                          s,   b,   p,   c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
inline ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     Index p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that p is valid:
  assert( 0 <= p );
  assert( p < mpr.mextent );

  return ConstTensor4View(mdata +
                          mpr.mstart + p * mpr.mstride,
                          msr, mbr, mrr, mcr,
                          s,   b,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
inline ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                                     Index b,
                                                     const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that b is valid:
  assert( 0 <= b );
  assert( b < mbr.mextent );

  return ConstTensor4View(mdata +
                          mbr.mstart + b * mbr.mstride,
                          msr, mpr, mrr, mcr,
                          s,   p,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
inline ConstTensor4View ConstTensor5View::operator()(Index s,
                                                     const Range& b,
                                                     const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that s is valid:
  assert( 0 <= s );
  assert( s < msr.mextent );

  return ConstTensor4View(mdata +
                          msr.mstart + s * msr.mstride,
                          mbr, mpr, mrr, mcr,
                          b,   p,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     const Range& p,
                                                     Index r,
                                                     Index c       ) const
{
  // Check that r and c is valid:
  assert( 0 <= r );
  assert( 0 <= c );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstTensor3View(mdata +
                          mrr.mstart + r * mrr.mstride +
                          mcr.mstart + c * mcr.mstride,
                          msr, mbr, mpr,
                          s,   b,   p);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     Index p,
                                                     const Range& r,
                                                     Index c       ) const
{
  // Check that p and c are valid:
  assert( 0 <= p );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstTensor3View(mdata +
                          mpr.mstart + p * mpr.mstride +
                          mcr.mstart + c * mcr.mstride,
                          msr, mbr, mrr,
                          s,   b,   r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(const Range& s,
                                                     const Range& b,
                                                     Index p,
                                                     Index r,
                                                     const Range& c) const
{
  // Check that p and r are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstTensor3View(mdata +
                          mpr.mstart + p * mpr.mstride +
                          mrr.mstart + r * mrr.mstride,
                          msr, mbr, mcr,
                          s,   b,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(const Range& s,
                                                     Index b,
                                                     const Range& p,
                                                     Index r,
                                                     const Range& c) const
{
  // Check that b and r are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );

  return ConstTensor3View(mdata +
                          mbr.mstart + b * mbr.mstride +
                          mrr.mstart + r * mrr.mstride,
                          msr, mpr, mcr,
                          s,   p,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(const Range& s,
                                                     Index b,
                                                     const Range& p,
                                                     const Range& r,
                                                     Index c       ) const
{
  // Check that b and c are valid:
  assert( 0 <= b );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( c < mcr.mextent );

  return ConstTensor3View(mdata +
                          mbr.mstart + b * mbr.mstride +
                          mcr.mstart + c * mcr.mstride,
                          msr, mpr, mrr,
                          s,   p,   r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(const Range& s,
                                                     Index b,
                                                     Index p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that b and p are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );

  return ConstTensor3View(mdata +
                          mbr.mstart + b * mbr.mstride +
                          mpr.mstart + p * mpr.mstride,
                          msr, mrr, mcr,
                          s,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(Index s,
                                                     const Range& b,
                                                     Index p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that s and p are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );

  return ConstTensor3View(mdata +
                          msr.mstart + s * msr.mstride +
                          mpr.mstart + p * mpr.mstride,
                          mbr, mrr, mcr,
                          b,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(Index s,
                                                     const Range& b,
                                                     const Range& p,
                                                     Index r,
                                                     const Range& c) const
{
  // Check that s and r are valid:
  assert( 0 <= s );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( r < mrr.mextent );

  return ConstTensor3View(mdata +
                          msr.mstart + s * msr.mstride +
                          mrr.mstart + r * mrr.mstride,
                          mbr, mpr, mcr,
                          b,   p,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(Index s,
                                                     const Range& b,
                                                     const Range& p,
                                                     const Range& r,
                                                     Index c       ) const
{
  // Check that s and c are valid:
  assert( 0 <= s );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( c < mcr.mextent );

  return ConstTensor3View(mdata +
                          msr.mstart + s * msr.mstride +
                          mcr.mstart + c * mcr.mstride,
                          mbr, mpr, mrr,
                          b,   p,   r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
inline ConstTensor3View ConstTensor5View::operator()(Index s,
                                                     Index b,
                                                     const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that s and b are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );

  return ConstTensor3View(mdata +
                          msr.mstart + s * msr.mstride +
                          mbr.mstart + b * mbr.mstride,
                          mpr, mrr, mcr,
                          p,   r,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(const Range& s,
                                                    const Range& b,
                                                    Index p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that p, r and c are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         msr, mbr,
                         s,   b);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(const Range& s,
                                                    Index b,
                                                    const Range& p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that b, r and c are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mrr.mstart + r * mpr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         msr, mpr,
                         s,   p);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(const Range& s,
                                                    Index b,
                                                    Index p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that b, p and c are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         msr, mrr,
                         s,   r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(const Range& s,
                                                    Index b,
                                                    Index p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that b, p and r are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstMatrixView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         msr, mcr,
                         s,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(Index s,
                                                    const Range& b,
                                                    Index p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that s, p and r are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstMatrixView(mdata +
                         msr.mstart + s * msr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         mbr, mcr,
                         b,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(Index s,
                                                    const Range& b,
                                                    Index p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that s, p and c are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         msr.mstart + s * msr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mbr, mrr,
                         b,   r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(Index s,
                                                    const Range& b,
                                                    const Range& p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that s, r and c are valid:
  assert( 0 <= s );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         msr.mstart + s * msr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mbr, mpr,
                         b,   p);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(Index s,
                                                    Index b,
                                                    const Range& p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that s, b and c are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         msr.mstart + s * msr.mstride +
                         mbr.mstart + b * mbr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mpr, mrr,
                         p,   r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(Index s,
                                                    Index b,
                                                    const Range& p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that s, b and r are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );

  return ConstMatrixView(mdata +
                         msr.mstart + s * msr.mstride +
                         mbr.mstart + b * mbr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         mpr, mcr,
                         p,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
inline ConstMatrixView ConstTensor5View::operator()(Index s,
                                                    Index b,
                                                    Index p,
                                                    const Range& r,
                                                    const Range& c) const
{
  // Check that s, b and p are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );

  return ConstMatrixView(mdata +
                         msr.mstart + s * msr.mstride +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride,
                         mrr, mcr,
                         r,   c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
inline ConstVectorView ConstTensor5View::operator()(const Range& s,
                                                    Index b,
                                                    Index p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that b, p, r and c are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         msr,
                         s);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
inline ConstVectorView ConstTensor5View::operator()(Index s,
                                                    const Range& b,
                                                    Index p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that s, p, r and c are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         msr.mstart + s * msr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mbr,
                         b);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
inline ConstVectorView ConstTensor5View::operator()(Index s,
                                                    Index b,
                                                    const Range& p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that s, b, r and c are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         msr.mstart + s * msr.mstride +
                         mbr.mstart + b * mbr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mpr,
                         p);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
inline ConstVectorView ConstTensor5View::operator()(Index s,
                                                    Index b,
                                                    Index p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that s, b, p and c are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         msr.mstart + s * msr.mstride +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mrr,
                         r);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
inline ConstVectorView ConstTensor5View::operator()(Index s,
                                                    Index b,
                                                    Index p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that s, b, p and r are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstVectorView(mdata +
                         msr.mstart + s * msr.mstride +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         mcr,
                         c);
}

/** Plain const index operator. */
inline Numeric ConstTensor5View::operator()(Index s,
                                            Index b,
                                            Index p,
                                            Index r,
                                            Index c) const
{
  // Check if indices are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return *( mdata +
            msr.mstart + s * msr.mstride +
            mbr.mstart + b * mbr.mstride +
            mpr.mstart + p * mpr.mstride +
            mrr.mstart + r * mrr.mstride +
            mcr.mstart + c * mcr.mstride );
}

/** Return const iterator to first shelf. */
inline ConstIterator5D ConstTensor5View::begin() const
{
  return ConstIterator5D( ConstTensor4View(mdata + msr.mstart,
                                           mbr, mpr, mrr, mcr),
                          msr.mstride );
}

/** Return const iterator behind last shelf. */
inline ConstIterator5D ConstTensor5View::end() const
{
  return ConstIterator5D( ConstTensor4View(mdata + msr.mstart +
                                           (msr.mextent) * msr.mstride,
                                           mbr, mpr, mrr, mcr),
                          msr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor5View::ConstTensor5View() :
  msr(0,0,1),
  mbr(0,0,1),
  mpr(0,0,1),
  mrr(0,0,1),
  mcr(0,0,1),
  mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor5 to initialize
    its own Tensor5View part. The book range br must have a stride to
    account for the length of one book. The shelf range sr must have a
    stride to account for the length of one shelf. */
inline ConstTensor5View::ConstTensor5View(Numeric *data,
                                          const Range& sr,
                                          const Range& br,
                                          const Range& pr,
                                          const Range& rr,
                                          const Range& cr) :
  msr(sr),
  mbr(br),
  mpr(pr),
  mrr(rr),
  mcr(cr),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
inline ConstTensor5View::ConstTensor5View(Numeric *data,
                                          const Range& ps,
                                          const Range& pb,
                                          const Range& pp,
                                          const Range& pr,
                                          const Range& pc,
                                          const Range& ns,
                                          const Range& nb,
                                          const Range& np,
                                          const Range& nr,
                                          const Range& nc) :
  msr(ps,ns),
  mbr(pb,nb),
  mpr(pp,np),
  mrr(pr,nr),
  mcr(pc,nc),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Tensor to print each book in turn. */
inline std::ostream& operator<<(std::ostream& os, const ConstTensor5View& v)
{
  // Page iterators:
  ConstIterator5D is = v.begin();
  const ConstIterator5D end_shelf = v.end();

  if ( is != end_shelf ) {
    os << *is;
    ++is;
  }

  for ( ; is != end_shelf; ++is ) {
    os << "\n\n";
    os << *is;
  }

  return os;
}


// Functions for Tensor5View:
// -------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor5. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline ConstTensor5View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor4View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                const Range& p,
                                                const Range& r,
                                                Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor4View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                const Range& p,
                                                Index r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor4View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                Index p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor4View Tensor5View::operator()(const Range& s,
                                                Index b,
                                                const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor4View Tensor5View::operator()(Index s,
                                                const Range& b,
                                                const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                const Range& p,
                                                Index r,
                                                Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                Index p,
                                                const Range& r,
                                                Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(const Range& s,
                                                const Range& b,
                                                Index p,
                                                Index r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(const Range& s,
                                                Index b,
                                                const Range& p,
                                                Index r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(const Range& s,
                                                Index b,
                                                const Range& p,
                                                const Range& r,
                                                Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(const Range& s,
                                                Index b,
                                                Index p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(Index s,
                                                const Range& b,
                                                Index p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(Index s,
                                                const Range& b,
                                                const Range& p,
                                                Index r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(Index s,
                                                const Range& b,
                                                const Range& p,
                                                const Range& r,
                                                Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor5View::operator()(Index s,
                                                Index b,
                                                const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(const Range& s,
                                               const Range& b,
                                               Index p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(const Range& s,
                                               Index b,
                                               const Range& p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(const Range& s,
                                               Index b,
                                               Index p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(const Range& s,
                                               Index b,
                                               Index p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(Index s,
                                               const Range& b,
                                               Index p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(Index s,
                                               const Range& b,
                                               Index p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(Index s,
                                               const Range& b,
                                               const Range& p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(Index s,
                                               Index b,
                                               const Range& p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(Index s,
                                               Index b,
                                               const Range& p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor5View::operator()(Index s,
                                               Index b,
                                               Index p,
                                               const Range& r,
                                               const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor5View::operator()(const Range& s,
                                               Index b,
                                               Index p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor5View::operator()(Index s,
                                               const Range& b,
                                               Index p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor5View::operator()(Index s,
                                               Index b,
                                               const Range& p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor5View::operator()(Index s,
                                               Index b,
                                               Index p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor5View::operator()(Index s,
                                               Index b,
                                               Index p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline Numeric Tensor5View::operator()(Index s,
                                       Index b,
                                       Index p,
                                       Index r,
                                       Index c) const
{
  return ConstTensor5View::operator()(s,b,p,r,c);
}

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor5. This allows
    correct recursive behavior.  */
inline Tensor5View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           const Range& p,
                                           const Range& r,
                                           const Range& c)
{
  return Tensor5View(mdata,
                     msr, mbr, mpr, mrr, mcr,
                     s,   b,   p,   r,   c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
inline Tensor4View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           const Range& p,
                                           const Range& r,
                                           Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return Tensor4View(mdata +
                     mcr.mstart + c * mcr.mstride,
                     msr, mbr, mpr, mrr,
                     s,   b,   p,   r);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
inline Tensor4View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           const Range& p,
                                           Index r,
                                           const Range& c)
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return Tensor4View(mdata +
                     mrr.mstart + r * mrr.mstride,
                     msr, mbr, mpr, mcr,
                     s,   b,   p,   c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
inline Tensor4View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           Index p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that p is valid:
  assert( 0 <= p );
  assert( p < mpr.mextent );

  return Tensor4View(mdata +
                     mpr.mstart + p * mpr.mstride,
                     msr, mbr, mrr, mcr,
                     s,   b,   r,   c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
inline Tensor4View Tensor5View::operator()(const Range& s,
                                           Index b,
                                           const Range& p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that b is valid:
  assert( 0 <= b );
  assert( b < mbr.mextent );

  return Tensor4View(mdata +
                     mbr.mstart + b * mbr.mstride,
                     msr, mpr, mrr, mcr,
                     s,   p,   r,   c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
inline Tensor4View Tensor5View::operator()(Index s,
                                           const Range& b,
                                           const Range& p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that s is valid:
  assert( 0 <= s );
  assert( s < msr.mextent );

  return Tensor4View(mdata +
                     msr.mstart + s * msr.mstride,
                     mbr, mpr, mrr, mcr,
                     b,   p,   r,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           const Range& p,
                                           Index r,
                                           Index c)
{
  // Check that r and c is valid:
  assert( 0 <= r );
  assert( 0 <= c );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return Tensor3View(mdata +
                     mrr.mstart + r * mrr.mstride +
                     mcr.mstart + c * mcr.mstride,
                     msr, mbr, mpr,
                     s,   b,   p);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           Index p,
                                           const Range& r,
                                           Index c)
{
  // Check that p and c are valid:
  assert( 0 <= p );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return Tensor3View(mdata +
                     mpr.mstart + p * mpr.mstride +
                     mcr.mstart + c * mcr.mstride,
                     msr, mbr, mrr,
                     s,   b,   r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(const Range& s,
                                           const Range& b,
                                           Index p,
                                           Index r,
                                           const Range& c)
{
  // Check that p and r are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return Tensor3View(mdata +
                     mpr.mstart + p * mpr.mstride +
                     mrr.mstart + r * mrr.mstride,
                     msr, mbr, mcr,
                     s,   b,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(const Range& s,
                                           Index b,
                                           const Range& p,
                                           Index r,
                                           const Range& c)
{
  // Check that b and r are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );

  return Tensor3View(mdata +
                     mbr.mstart + b * mbr.mstride +
                     mrr.mstart + r * mrr.mstride,
                     msr, mpr, mcr,
                     s,   p,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(const Range& s,
                                           Index b,
                                           const Range& p,
                                           const Range& r,
                                           Index c)
{
  // Check that b and c are valid:
  assert( 0 <= b );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( c < mcr.mextent );

  return Tensor3View(mdata +
                     mbr.mstart + b * mbr.mstride +
                     mcr.mstart + c * mcr.mstride,
                     msr, mpr, mrr,
                     s,   p,   r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(const Range& s,
                                           Index b,
                                           Index p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that b and p are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );

  return Tensor3View(mdata +
                     mbr.mstart + b * mbr.mstride +
                     mpr.mstart + p * mpr.mstride,
                     msr, mrr, mcr,
                     s,   r,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(Index s,
                                           const Range& b,
                                           Index p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that s and p are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );

  return Tensor3View(mdata +
                     msr.mstart + s * msr.mstride +
                     mpr.mstart + p * mpr.mstride,
                     mbr, mrr, mcr,
                     b,   r,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(Index s,
                                           const Range& b,
                                           const Range& p,
                                           Index r,
                                           const Range& c)
{
  // Check that s and r are valid:
  assert( 0 <= s );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( r < mrr.mextent );

  return Tensor3View(mdata +
                     msr.mstart + s * msr.mstride +
                     mrr.mstart + r * mrr.mstride,
                     mbr, mpr, mcr,
                     b,   p,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(Index s,
                                           const Range& b,
                                           const Range& p,
                                           const Range& r,
                                           Index c)
{
  // Check that s and c are valid:
  assert( 0 <= s );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( c < mcr.mextent );

  return Tensor3View(mdata +
                     msr.mstart + s * msr.mstride +
                     mcr.mstart + c * mcr.mstride,
                     mbr, mpr, mrr,
                     b,   p,   r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
inline Tensor3View Tensor5View::operator()(Index s,
                                           Index b,
                                           const Range& p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that s and b are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );

  return Tensor3View(mdata +
                     msr.mstart + s * msr.mstride +
                     mbr.mstart + b * mbr.mstride,
                     mpr, mrr, mcr,
                     p,   r,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(const Range& s,
                                          const Range& b,
                                          Index p,
                                          Index r,
                                          Index c)
{
  // Check that p, r and c are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    msr, mbr,
                    s,   b);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(const Range& s,
                                          Index b,
                                          const Range& p,
                                          Index r,
                                          Index c)
{
  // Check that b, r and c are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mrr.mstart + r * mpr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    msr, mpr,
                    s,   p);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(const Range& s,
                                          Index b,
                                          Index p,
                                          const Range& r,
                                          Index c)
{
  // Check that b, p and c are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    msr, mrr,
                    s,   r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(const Range& s,
                                          Index b,
                                          Index p,
                                          Index r,
                                          const Range& c)
{
  // Check that b, p and r are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return MatrixView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    msr, mcr,
                    s,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(Index s,
                                          const Range& b,
                                          Index p,
                                          Index r,
                                          const Range& c)
{
  // Check that s, p and r are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return MatrixView(mdata +
                    msr.mstart + s * msr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    mbr, mcr,
                    b,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(Index s,
                                          const Range& b,
                                          Index p,
                                          const Range& r,
                                          Index c       )
{
  // Check that s, p and c are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    msr.mstart + s * msr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mbr, mrr,
                    b,   r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(Index s,
                                          const Range& b,
                                          const Range& p,
                                          Index r,
                                          Index c       )
{
  // Check that s, r and c are valid:
  assert( 0 <= s );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    msr.mstart + s * msr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mbr, mpr,
                    b,   p);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(Index s,
                                          Index b,
                                          const Range& p,
                                          const Range& r,
                                          Index c       )
{
  // Check that s, b and c are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    msr.mstart + s * msr.mstride +
                    mbr.mstart + b * mbr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mpr, mrr,
                    p,   r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(Index s,
                                          Index b,
                                          const Range& p,
                                          Index r,
                                          const Range& c)
{
  // Check that s, b and r are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );

  return MatrixView(mdata +
                    msr.mstart + s * msr.mstride +
                    mbr.mstart + b * mbr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    mpr, mcr,
                    p,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
inline MatrixView Tensor5View::operator()(Index s,
                                          Index b,
                                          Index p,
                                          const Range& r,
                                          const Range& c)
{
  // Check that s, b and p are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );

  return MatrixView(mdata +
                    msr.mstart + s * msr.mstride +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride,
                    mrr, mcr,
                    r,   c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
inline VectorView Tensor5View::operator()(const Range& s,
                                          Index b,
                                          Index p,
                                          Index r,
                                          Index c)
{
  // Check that b, p, r and c are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    msr,
                    s);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
inline VectorView Tensor5View::operator()(Index s,
                                          const Range& b,
                                          Index p,
                                          Index r,
                                          Index c)
{
  // Check that s, p, r and c are valid:
  assert( 0 <= s );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    msr.mstart + s * msr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mbr,
                    b);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
inline VectorView Tensor5View::operator()(Index s,
                                          Index b,
                                          const Range& p,
                                          Index r,
                                          Index c)
{
  // Check that s, b, r and c are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    msr.mstart + s * msr.mstride +
                    mbr.mstart + b * mbr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mpr,
                    p);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
inline VectorView Tensor5View::operator()(Index s,
                                          Index b,
                                          Index p,
                                          const Range& r,
                                          Index c)
{
  // Check that s, b, p and c are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    msr.mstart + s * msr.mstride +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mrr,
                    r);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
inline VectorView Tensor5View::operator()(Index s,
                                          Index b,
                                          Index p,
                                          Index r,
                                          const Range& c)
{
  // Check that s, b, p and r are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return VectorView(mdata +
                    msr.mstart + s * msr.mstride +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    mcr,
                    c);
}

/** Plain const index operator. */
inline Numeric& Tensor5View::operator()(Index s,
                                        Index b,
                                        Index p,
                                        Index r,
                                        Index c)
{
  // Check if indices are valid:
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( s < msr.mextent );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return *( mdata +
            msr.mstart + s * msr.mstride +
            mbr.mstart + b * mbr.mstride +
            mpr.mstart + p * mpr.mstride +
            mrr.mstart + r * mrr.mstride +
            mcr.mstart + c * mcr.mstride );
}
/** Return const iterator to first shelf. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator5D Tensor5View::begin() const
{
  return ConstTensor5View::begin();
}

/** Return const iterator behind last shelf. */
inline ConstIterator5D Tensor5View::end() const
{
  return ConstTensor5View::end();
}

/** Return iterator to first shelf. */
inline Iterator5D Tensor5View::begin()
{
  return Iterator5D( Tensor4View(mdata + msr.mstart,
                                 mbr, mpr, mrr, mcr),
                     msr.mstride );
}

/** Return iterator behind last shelf. */
inline Iterator5D Tensor5View::end()
{
  return Iterator5D( Tensor4View(mdata + msr.mstart +
                                 (msr.mextent) * msr.mstride,
                                 mbr, mpr, mrr, mcr),
                     msr.mstride );
}

/** Assignment operator. This copies the data from another Tensor5View
    to this Tensor5View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor5View by
    setting its range. */
inline Tensor5View& Tensor5View::operator=(const ConstTensor5View& m)
{
  // Check that sizes are compatible:
  assert( msr.mextent == m.msr.mextent );
  assert( mbr.mextent == m.mbr.mextent );
  assert( mpr.mextent == m.mpr.mextent );
  assert( mrr.mextent == m.mrr.mextent );
  assert( mcr.mextent == m.mcr.mextent );

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from Tensor5View to Tensor5View. This is a tricky
    one. The problem is that since Tensor5View is derived from
    ConstTensor5View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline Tensor5View& Tensor5View::operator=(const Tensor5View& m)
{
  // Check that sizes are compatible:
  assert( msr.mextent == m.msr.mextent );
  assert( mbr.mextent == m.mbr.mextent );
  assert( mpr.mextent == m.mpr.mextent );
  assert( mrr.mextent == m.mrr.mextent );
  assert( mcr.mextent == m.mcr.mextent );

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Tensor5. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline Tensor5View& Tensor5View::operator=(const Tensor5& m)
{
  // Check that sizes are compatible:
  assert( msr.mextent == m.msr.mextent );
  assert( mbr.mextent == m.mbr.mextent );
  assert( mpr.mextent == m.mpr.mextent );
  assert( mrr.mextent == m.mrr.mextent );
  assert( mcr.mextent == m.mcr.mextent );

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor5View will set all elements to this
    value. */
inline Tensor5View& Tensor5View::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
inline Tensor5View& Tensor5View::operator*=(Numeric x)
{
  const Iterator5D es = end();
  for ( Iterator5D s = begin(); s != es ; ++s )
  {
    *s *= x;
  }
  return *this;
}

/** Division by scalar. */
inline Tensor5View& Tensor5View::operator/=(Numeric x)
{
  const Iterator5D es = end();
  for ( Iterator5D s = begin(); s != es ; ++s )
  {
    *s /= x;
  }
  return *this;
}

/** Addition of scalar. */
inline Tensor5View& Tensor5View::operator+=(Numeric x)
{
  const Iterator5D es = end();
  for ( Iterator5D s = begin(); s != es ; ++s )
  {
    *s += x;
  }
  return *this;
}

/** Subtraction of scalar. */
inline Tensor5View& Tensor5View::operator-=(Numeric x)
{
  const Iterator5D es = end();
  for ( Iterator5D s = begin(); s != es ; ++s )
  {
    *s -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor5. */
inline Tensor5View& Tensor5View::operator*=(const ConstTensor5View& x)
{
  assert( nshelves() == x.nshelves() );
  assert( nbooks()   == x.nbooks()   );
  assert( npages()   == x.npages()   );
  assert( nrows()    == x.nrows()    );
  assert( ncols()    == x.ncols()    );
  ConstIterator5D  xs = x.begin();
  Iterator5D        s = begin();
  const Iterator5D es = end();
  for ( ; s != es ; ++s, ++xs )
    {
      *s *= *xs;
    }
  return *this;
}

/** Element-vise division by another Tensor5. */
inline Tensor5View& Tensor5View::operator/=(const ConstTensor5View& x)
{
  assert( nshelves() == x.nshelves() );
  assert( nbooks() == x.nbooks()     );
  assert( npages() == x.npages()     );
  assert( nrows()  == x.nrows()      );
  assert( ncols()  == x.ncols()      );
  ConstIterator5D  xs = x.begin();
  Iterator5D        s = begin();
  const Iterator5D es = end();
  for ( ; s != es ; ++s, ++xs )
    {
      *s /= *xs;
    }
  return *this;
}

/** Element-vise addition of another Tensor5. */
inline Tensor5View& Tensor5View::operator+=(const ConstTensor5View& x)
{
  assert( nshelves() == x.nshelves() );
  assert( nbooks() == x.nbooks()     );
  assert( npages() == x.npages()     );
  assert( nrows()  == x.nrows()      );
  assert( ncols()  == x.ncols()      );
  ConstIterator5D  xs = x.begin();
  Iterator5D        s = begin();
  const Iterator5D es = end();
  for ( ; s != es ; ++s, ++xs )
    {
      *s += *xs;
    }
  return *this;
}

/** Element-vise subtraction of another Tensor5. */
inline Tensor5View& Tensor5View::operator-=(const ConstTensor5View& x)
{
  assert( nshelves() == x.nshelves() );
  assert( nbooks() == x.nbooks()     );
  assert( npages() == x.npages()     );
  assert( nrows()  == x.nrows()      );
  assert( ncols()  == x.ncols()      );
  ConstIterator5D  xs = x.begin();
  Iterator5D        s = begin();
  const Iterator5D es = end();
  for ( ; s != es ; ++s, ++xs )
    {
      *s -= *xs;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Tensor5. */
inline Tensor5View::Tensor5View() :
  ConstTensor5View()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor5 to initialize its
    own Tensor5View part. The row range rr must have a
    stride to account for the length of one row. */
inline Tensor5View::Tensor5View(Numeric *data,
                                const Range& sr,
                                const Range& br,
                                const Range& pr,
                                const Range& rr,
                                const Range& cr) :
  ConstTensor5View(data, sr, br, pr, rr, cr)
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
inline Tensor5View::Tensor5View(Numeric *data,
                                const Range& ps,
                                const Range& pb,
                                const Range& pp,
                                const Range& pr,
                                const Range& pc,
                                const Range& ns,
                                const Range& nb,
                                const Range& np,
                                const Range& nr,
                                const Range& nc) :
  ConstTensor5View(data, ps, pb, pp, pr, pc, ns, nb, np, nr, nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
inline void copy(ConstIterator5D origin,
                 const ConstIterator5D& end,
                 Iterator5D target)
{
  for ( ; origin != end ; ++origin, ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy( origin->begin(), origin->end(), target->begin() );
    }
}

/** Copy a scalar to all elements. */
inline void copy(Numeric x,
                 Iterator5D target,
                 const Iterator5D& end)
{
  for ( ; target != end ; ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy( x, target->begin(), target->end() );
    }
}


// Functions for Tensor5:
// ---------------------

/** Default constructor. */
inline Tensor5::Tensor5() :
  Tensor5View::Tensor5View()
{
  // Nothing to do here. However, note that the default constructor
  // for Tensor5View has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized.
}

/** Constructor setting size. This constructor has to set the strides
    in the shelf, book, page and row ranges correctly! */
inline Tensor5::Tensor5(Index s, Index b, Index p, Index r, Index c) :
  Tensor5View( new Numeric[s*b*p*r*c],
               Range( 0, s, b*p*r*c ),
               Range( 0, b, p*r*c ),
               Range( 0, p, r*c ),
               Range( 0, r, c ),
               Range( 0, c) )
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Tensor5::Tensor5(Index s, Index b, Index p, Index r, Index c, Numeric fill) :
  Tensor5View( new Numeric[s*b*p*r*c],
               Range( 0, s, b*p*r*c ),
               Range( 0, b, p*r*c ),
               Range( 0, p, r*c ),
               Range( 0, r, c ),
               Range( 0, c) )
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x = mdata; x < mdata + s*b*p*r*c; ++x )
    *x = fill;
}

/** Copy constructor from Tensor5View. This automatically sets the size
    and copies the data. */
inline Tensor5::Tensor5(const ConstTensor5View& m) :
  Tensor5View( new Numeric[m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols()],
               Range( 0, m.nshelves(), m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
               Range( 0, m.nbooks(), m.npages()*m.nrows()*m.ncols() ),
               Range( 0, m.npages(), m.nrows()*m.ncols() ),
               Range( 0, m.nrows(), m.ncols() ),
               Range( 0, m.ncols() ) )
{
  copy( m.begin(), m.end(), begin() );
}

/** Copy constructor from Tensor5. This automatically sets the size
    and copies the data. */
inline Tensor5::Tensor5(const Tensor5& m) :
  Tensor5View( new Numeric[m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols()],
               Range( 0, m.nshelves(), m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
               Range( 0, m.nbooks(), m.npages()*m.nrows()*m.ncols() ),
               Range( 0, m.npages(), m.nrows()*m.ncols() ),
               Range( 0, m.nrows(), m.ncols() ),
               Range( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty tensor, then it will have
  // dimensions of size 0. But these are used to initialize the stride
  // for higher dimensions! Thus, this method has to be consistent
  // with the behaviour of Range::Range. For now, Range::Range allows
  // also stride 0.
  copy( m.begin(), m.end(), begin() );
}

/** Assignment operator from another tensor. It is important that this
    operator exists. Otherwise the = operator seems to copy references
    instead of content in some cases.

    The Behavior of this one is a bit special: If the size of the
    target tensor is 0 then it will be automatically resized to match
    (this is needed to have the correct initialization for constructed
    classes that use the assignment operator to initialize their
    data).
*/
inline Tensor5& Tensor5::operator=(const Tensor5& m)
{
  //  cout << "Tensor5 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "              n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid tensor, so we just
  // have to check one.
  if ( 0 == mcr.mextent )
    {
      // Adjust if previously empty.
      resize( m.msr.mextent, m.mbr.mextent, m.mpr.mextent, m.mrr.mextent, m.mcr.mextent );
    }
  else
    {
      // Check that sizes are compatible:
      assert( msr.mextent == m.msr.mextent );
      assert( mbr.mextent == m.mbr.mextent );
      assert( mpr.mextent == m.mpr.mextent );
      assert( mrr.mextent == m.mrr.mextent );
      assert( mcr.mextent == m.mcr.mextent );
    }

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
inline Tensor5& Tensor5::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
inline void Tensor5::resize(Index s, Index b, Index p, Index r, Index c)
{
  assert( 0 <= s );
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );

  if ( msr.mextent != s ||
       mbr.mextent != b ||
       mpr.mextent != p ||
       mrr.mextent != r ||
       mcr.mextent != c )
    {
      delete mdata;
      mdata = new Numeric[s*b*p*r*c];

      msr.mstart = 0;
      msr.mextent = s;
      msr.mstride = b*p*r*c;

      mbr.mstart = 0;
      mbr.mextent = b;
      mbr.mstride = p*r*c;

      mpr.mstart = 0;
      mpr.mextent = p;
      mpr.mstride = r*c;

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
}

/** Destructor for Tensor5. This is important, since Tensor5 uses new to
    allocate storage. */
inline Tensor5::~Tensor5()
{
//   cout << "Destroying a Tensor5:\n"
//        << *this << "\n........................................\n";
  delete mdata;
}


/** A generic transform function for tensors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for tensors! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    The two views may be the same one, in which case the
    conversion happens in place.

    \retval   y   the results of the function acting on each element of x
    \param    my_func a function (e.g., sqrt)
    \param    x   a tensor */
inline void transform( Tensor5View y,
                       double (&my_func)(double),
                       ConstTensor5View x )
{
  // Check dimensions:
  assert( y.nshelves() == x.nshelves() );
  assert( y.nbooks()   == x.nbooks()   );
  assert( y.npages()   == x.npages()   );
  assert( y.nrows()    == x.nrows()    );
  assert( y.ncols()    == x.ncols()    );

  const ConstIterator5D xe = x.end();
  ConstIterator5D       xi = x.begin();
  Iterator5D            yi = y.begin();
  for ( ; xi != xe; ++xi, ++yi )
    {
      // Use the transform function of lower dimensional tensors
      // recursively:
      transform( *yi, my_func, *xi );
    }
}

/** Max function, tensor version. */
inline Numeric max(const ConstTensor5View& x)
{
  const ConstIterator5D xe = x.end();
  ConstIterator5D       xi = x.begin();

  // Initial value for max:
  Numeric themax = max(*xi);
  ++xi;

  for ( ; xi != xe ; ++xi )
    {
      // Use the max function of lower dimensional tensors
      // recursively:
      Numeric maxi = max(*xi);
      if ( maxi > themax )
        themax = maxi;
    }

  return themax;
}

/** Min function, tensor version. */
inline Numeric min(const ConstTensor5View& x)
{
  const ConstIterator5D xe = x.end();
  ConstIterator5D       xi = x.begin();

  // Initial value for min:
  Numeric themin = min(*xi);
  ++xi;

  for ( ; xi != xe ; ++xi )
    {
      // Use the min function of lower dimensional tensors
      // recursively:
      Numeric mini = min(*xi);
      if ( mini < themin )
        themin = mini;
    }

  return themin;
}

#endif    // matpackV_h
