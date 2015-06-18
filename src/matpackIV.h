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
  Implementation of Tensors of Rank 4.

  Based on Tensor3 by Stefan Buehler.

  The four dimensions are called: book, page, row, column.

  \author Wolfram-Andre Haas
  \date   2002-03-01
 */

#ifndef matpackIV_h
#define matpackIV_h

#include "matpackIII.h"

/** The outermost iterator class for rank 4 tensors. This takes into
    account the defined strided. */
class Iterator4D {
public:
  // Constructors:
  /** Default constructor. */
  Iterator4D() : msv(), mstride(0) { /* Nothing to do here. */ }

  /** Explicit constructor. */
  Iterator4D(const Tensor3View& x, Index stride)
    : msv(x), mstride(stride)
      { /* Nothing to do here. */ }

  // Operators:
  /** Prefix increment operator. */
  Iterator4D& operator++() { msv.mdata += mstride; return *this; }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const Iterator4D& other) const
    { if ( msv.mdata +
           msv.mpr.mstart +
           msv.mrr.mstart +
           msv.mcr.mstart
           !=
           other.msv.mdata +
           other.msv.mpr.mstart +
           other.msv.mrr.mstart +
           other.msv.mcr.mstart )
        return true;
      else
        return false;
    }

  Tensor3View* operator->();
  Tensor3View& operator*();

private:
  /** Current position. */
  Tensor3View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator4D. */
class ConstIterator4D {
public:
  // Constructors:
  // Functions for ConstIterator4D
  // -----------------------------

  /** Default constructor. */
  ConstIterator4D() : msv(), mstride(0) { /* Nothing to do here. */ }

  /** Explicit constructor. */
  ConstIterator4D(const ConstTensor3View& x, Index stride)
    : msv(x), mstride(stride)
      { /* Nothing to do here. */ }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator4D& operator++() { msv.mdata += mstride; return *this; }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const ConstIterator4D& other) const
    { if ( msv.mdata +
           msv.mpr.mstart +
           msv.mrr.mstart +
           msv.mcr.mstart
           !=
           other.msv.mdata +
           other.msv.mpr.mstart +
           other.msv.mrr.mstart +
           other.msv.mcr.mstart )
        return true;
      else
        return false;
    }

  const ConstTensor3View* operator->() const;
  const ConstTensor3View& operator*()  const;

private:
  /** Current position. */
  ConstTensor3View msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor4:
class Tensor4;


/** A constant view of a Tensor4.

    This, together with the derived class Tensor4View, contains the
    main implementation of a Tensor4. It defines the concepts of
    Tensor4View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor4View from a subrange of
    a Tensor4View.

    The four dimensions of the tensor are called: book, page, row, column.

    The class Tensor4 is just a special case of a Tensor4View
    which also allocates storage. */
class ConstTensor4View {
public:
  // Member functions:

  bool empty() const;

  Index nbooks() const;
  Index npages() const;
  Index nrows()  const;
  Index ncols()  const;

  // Const index operators:
  ConstTensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& b, Index p,        const Range& r, const Range& c ) const;
  ConstTensor3View operator()( Index b,        const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index b,        Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( const Range& b, Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        const Range& p, Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        Index r,        const Range& c ) const;

  /** Plain const index operator. */
  Numeric operator()(Index b, Index p, Index r, Index c) const
    { // Check if indices are valid:
      assert( 0 <= b );
      assert( 0 <= p );
      assert( 0 <= r );
      assert( 0 <= c );
      assert( b < mbr.mextent );
      assert( p < mpr.mextent );
      assert( r < mrr.mextent );
      assert( c < mcr.mextent );

      return get(b, p, r, c);
    }

  /** Get element implementation without assertions. */
  Numeric get(Index b, Index p, Index r, Index c) const
    {
      return *( mdata +
                mbr.mstart + b * mbr.mstride +
                mpr.mstart + p * mpr.mstride +
                mrr.mstart + r * mrr.mstride +
                mcr.mstart + c * mcr.mstride );
    }

  // Functions returning iterators:
  ConstIterator4D begin() const;
  ConstIterator4D end()   const;

  //! Destructor
  virtual ~ConstTensor4View() {}

  // Friends:
  friend class Tensor4View;
  friend class ConstIterator5D;
  friend class ConstTensor5View;
  friend class ConstTensor6View;
  friend class ConstTensor7View;

  // Special constructor to make a Tensor4 view of a Tensor3.
  ConstTensor4View(const ConstTensor3View& a);

protected:
  // Constructors:
  ConstTensor4View();
  ConstTensor4View(Numeric *data,
                   const Range& b, const Range& p, const Range& r, const Range& c);
  ConstTensor4View(Numeric *data,
                   const Range& pb, const Range& pp, const Range& pr, const Range& pc,
                   const Range& nb, const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
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

/** The Tensor4View class

    This contains the main implementation of a Tensor4. It defines
    the concepts of Tensor4View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor4View
    from a subrange of a Tensor4View.

    The class Tensor4 is just a special case of a Tensor4View
    which also allocates storage. */
class Tensor4View : public ConstTensor4View {
public:

  // Const index operators:
  ConstTensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& b, Index p,        const Range& r, const Range& c ) const;
  ConstTensor3View operator()( Index b,        const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index b,        Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( const Range& b, Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        const Range& p, Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        Index r,        const Range& c ) const;

  /** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
  Numeric operator()(Index b, Index p, Index r, Index c) const
    { return ConstTensor4View::operator()(b,p,r,c); }

  /** Get element implementation without assertions. */
  Numeric get(Index b, Index p, Index r, Index c) const
    { return ConstTensor4View::get(b,p,r,c); }

  // Non-const index operators:

  Tensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c );

  Tensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        );
  Tensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c );
  Tensor3View operator()( const Range& b, Index p,        const Range& r, const Range& c );
  Tensor3View operator()( Index b,        const Range& p, const Range& r, const Range& c );

  MatrixView  operator()( const Range& b, const Range& p, Index r,        Index c        );
  MatrixView  operator()( const Range& b, Index p,        const Range& r, Index c        );
  MatrixView  operator()( const Range& b, Index p,        Index r,        const Range& c );
  MatrixView  operator()( Index b,        const Range& p, Index r,        const Range& c );
  MatrixView  operator()( Index b,        const Range& p, const Range& r, Index c        );
  MatrixView  operator()( Index b,        Index p,        const Range& r, const Range& c );

  VectorView  operator()( const Range& b, Index p,        Index r,        Index c        );
  VectorView  operator()( Index b,        const Range& p, Index r,        Index c        );
  VectorView  operator()( Index b,        Index p,        const Range& r, Index c        );
  VectorView  operator()( Index b,        Index p,        Index r,        const Range& c );

  /** Plain non-const index operator. */
  Numeric& operator()(Index b, Index p, Index r, Index c)
    { // Check if indices are valid:
      assert( 0 <= b );
      assert( 0 <= p );
      assert( 0 <= r );
      assert( 0 <= c );
      assert( b < mbr.mextent );
      assert( p < mpr.mextent );
      assert( r < mrr.mextent );
      assert( c < mcr.mextent );

      return get(b, p, r, c);
    }

  /** Get element implementation without assertions. */
  Numeric& get(Index b, Index p, Index r, Index c)
    {
      return *( mdata +
                mbr.mstart + b * mbr.mstride +
                mpr.mstart + p * mpr.mstride +
                mrr.mstart + r * mrr.mstride +
                mcr.mstart + c * mcr.mstride );
    }

  // Conversion to a plain C-array
  const Numeric *get_c_array() const;
  Numeric *get_c_array();

  // Functions returning const iterators:
  ConstIterator4D begin() const;
  ConstIterator4D end()   const;
  // Functions returning iterators:
  Iterator4D begin();
  Iterator4D end();

  // Assignment operators:
  Tensor4View& operator=(const ConstTensor4View& v);
  Tensor4View& operator=(const Tensor4View& v);
  Tensor4View& operator=(const Tensor4& v);
  Tensor4View& operator=(Numeric x);

  // Other operators:
  Tensor4View& operator*=(Numeric x);
  Tensor4View& operator/=(Numeric x);
  Tensor4View& operator+=(Numeric x);
  Tensor4View& operator-=(Numeric x);

  Tensor4View& operator*=(const ConstTensor4View& x);
  Tensor4View& operator/=(const ConstTensor4View& x);
  Tensor4View& operator+=(const ConstTensor4View& x);
  Tensor4View& operator-=(const ConstTensor4View& x);

  //! Destructor
  virtual ~Tensor4View() {}

  // Friends:
  // friend class VectorView;
  // friend ConstTensor4View transpose(ConstTensor4View m);
  // friend Tensor4View transpose(Tensor4View m);
  friend class Iterator5D;
  friend class Tensor5View;
  friend class Tensor6View;
  friend class Tensor7View;

  // Special constructor to make a Tensor4 view of a Tensor3.
  Tensor4View(const Tensor3View& a);

protected:
  // Constructors:
  Tensor4View();
  Tensor4View(Numeric *data,
              const Range& b, const Range& p, const Range& r, const Range& c);
  Tensor4View(Numeric *data,
              const Range& pb, const Range& pp, const Range& pr, const Range& pc,
              const Range& nb, const Range& np, const Range& nr, const Range& nc);
};


/** The Tensor4 class. This is a Tensor4View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor4View. Additionally defined here
    are:

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor4 : public Tensor4View {
public:
  // Constructors:
  Tensor4();
  Tensor4(Index b, Index p, Index r, Index c);
  Tensor4(Index b, Index p, Index r, Index c, Numeric fill);
  Tensor4(const ConstTensor4View& v);
  Tensor4(const Tensor4& v);

  // Assignment operators:
  Tensor4& operator=(Tensor4 x);
  Tensor4& operator=(Numeric x);

  // Resize function:
  void resize(Index b, Index p, Index r, Index c);

  // Swap function:
  friend void swap(Tensor4& t1, Tensor4& t2);

  // Destructor:
  virtual ~Tensor4();
};


// Function declarations:
// ----------------------

void copy(ConstIterator4D origin,
          const ConstIterator4D& end,
          Iterator4D target);

void copy(Numeric x,
          Iterator4D target,
          const Iterator4D& end);

void transform( Tensor4View y,
                double (&my_func)(double),
                ConstTensor4View x );

Numeric max(const ConstTensor4View& x);

Numeric min(const ConstTensor4View& x);

std::ostream& operator<<(std::ostream& os, const ConstTensor4View& v);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor4view_get_elem (Tensor4View& tv,
                                    Index b, Index p, Index r, Index c);

#endif
////////////////////////////////

#endif    // matpackIV_h
