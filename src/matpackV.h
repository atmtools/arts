/* Copyright (C) 2002-2012
   Stefan Buehler <sbuehler@ltu.se>
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

#include "matpackIV.h"

/** The outermost iterator class for rank 5 tensors. This takes into
    account the defined strided. */
class Iterator5D {
public:
  // Constructors:
  // Functions for Iterator5D
  // ------------------------

  /** Default constructor. */
  Iterator5D() : msv(), mstride(0) { /* Nothing to do here. */ }

  /** Explicit constructor. */
  Iterator5D(const Tensor4View& x, Index stride) : msv(x), mstride(stride)
      { /* Nothing to do here. */ }

  // Operators:
  /** Prefix increment operator. */
  Iterator5D& operator++() { msv.mdata += mstride; return *this; }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const Iterator5D& other) const
    { if ( msv.mdata +
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
  Tensor4View* operator->() { return &msv; }

  /** Dereferencing. */
  Tensor4View& operator*() { return msv; }

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
  /** Default constructor. */
  ConstIterator5D() : msv(), mstride(0) { /* Nothing to do here. */ }

  /** Explicit constructor. */
  ConstIterator5D(const ConstTensor4View& x, Index stride)
    : msv(x), mstride(stride)
      { /* Nothing to do here. */ }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator5D& operator++() { msv.mdata += mstride; return *this; }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const ConstIterator5D& other) const
    { if ( msv.mdata +
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
  const ConstTensor4View* operator->() const { return &msv; }

  /** Dereferencing. */
  const ConstTensor4View& operator*() const { return msv; }

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
  bool empty() const;
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

  /** Plain const index operator. */
  Numeric operator()(Index s, Index b, Index p, Index r, Index c) const
    { // Check if indices are valid:
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

      return get(s, b, p, r, c);
    }

  /** Get element implementation without assertions. */
  Numeric get(Index s, Index b, Index p, Index r, Index c) const
    {
      return *( mdata +
                msr.mstart + s * msr.mstride +
                mbr.mstart + b * mbr.mstride +
                mpr.mstart + p * mpr.mstride +
                mrr.mstart + r * mrr.mstride +
                mcr.mstart + c * mcr.mstride );
    }

  // Functions returning iterators:
  ConstIterator5D begin() const;
  ConstIterator5D end()   const;

  //! Destructor
  virtual ~ConstTensor5View() {}

  // Friends:
  friend class Tensor5View;
  friend class ConstIterator6D;
  friend class ConstTensor6View;
  friend class ConstTensor7View;

  // Special constructor to make a Tensor5 view of a Tensor4.
  ConstTensor5View(const ConstTensor4View& a);

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

  /** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
  Numeric operator()(Index s, Index b, Index p, Index r, Index c) const
    { return ConstTensor5View::operator()(s,b,p,r,c); }

  /** Get element implementation without assertions. */
  Numeric get(Index s, Index b, Index p, Index r, Index c) const
    { return ConstTensor5View::get(s,b,p,r,c); }

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

  /** Plain const index operator. */
  Numeric& operator()(Index s, Index b, Index p, Index r, Index c)
    { // Check if indices are valid:
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

      return get(s, b, p, r, c);
    }

  /** Get element implementation without assertions. */
  Numeric& get(Index s, Index b, Index p, Index r, Index c)
    {
      return *( mdata +
                msr.mstart + s * msr.mstride +
                mbr.mstart + b * mbr.mstride +
                mpr.mstart + p * mpr.mstride +
                mrr.mstart + r * mrr.mstride +
                mcr.mstart + c * mcr.mstride );
    }

  // Conversion to a plain C-array
  const Numeric *get_c_array() const;
  Numeric *get_c_array();
  
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

  //! Destructor
  virtual ~Tensor5View() {}

  // Friends:
  // friend class VectorView;
  // friend ConstTensor5View transpose(ConstTensor5View m);
  // friend Tensor5View transpose(Tensor5View m);
  friend class Iterator6D;
  friend class Tensor6View;
  friend class Tensor7View;

  // Special constructor to make a Tensor5 view of a Tensor4.
  Tensor5View(const Tensor4View& a);

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
  Tensor5& operator=(Tensor5 x);
  Tensor5& operator=(Numeric x);

  // Resize function:
  void resize(Index s, Index b, Index p, Index r, Index c);

  // Swap function:
  friend void swap(Tensor5& t1, Tensor5& t2);

  // Destructor:
  virtual ~Tensor5();
};


// Function declarations:
// ----------------------

void copy(ConstIterator5D origin,
          const ConstIterator5D& end,
          Iterator5D target);

void copy(Numeric x,
          Iterator5D target,
          const Iterator5D& end);

void transform( Tensor5View y,
                double (&my_func)(double),
                ConstTensor5View x );

Numeric max(const ConstTensor5View& x);

Numeric min(const ConstTensor5View& x);

std::ostream& operator<<(std::ostream& os, const ConstTensor5View& v);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor5view_get_elem (Tensor5View& tv, Index s, Index b, Index p,
                                    Index r, Index c);

#endif
////////////////////////////////

#endif    // matpackV_h
