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
  Implementation of Tensors of Rank 3.

  The three dimensions are called: page, row, column.
  
  \author Stefan Buehler
  \date   2001-11-22
 */

#ifndef matpackIII_h
#define matpackIII_h

#include <iomanip>
#include "matpackI.h"

/** The outermost iterator class for rank 3 tensors. This takes into
    account the defined strided. */
class Iterator3D {
public:
  // Constructors:
  Iterator3D();
  Iterator3D(const Iterator3D& o);
  Iterator3D(const MatrixView& x, Index stride);

  // Operators:
  Iterator3D& operator++();
  bool operator!=(const Iterator3D& other) const;

  /** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
  MatrixView* const Iterator3D::operator->() { return &msv; }

  /** Dereferencing. */
  MatrixView& Iterator3D::operator*() { return msv; }
 
private:
  /** Current position. */
  MatrixView msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator3D. */
class ConstIterator3D {
public:
  // Constructors:
  ConstIterator3D();
  ConstIterator3D(const ConstIterator3D& o);
  ConstIterator3D(const ConstMatrixView& x, Index stride);

  // Operators:
  ConstIterator3D& operator++();
  bool operator!=(const ConstIterator3D& other) const;

  /** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
  const ConstMatrixView* ConstIterator3D::operator->() const { return &msv; }

  /** Dereferencing. */
  const ConstMatrixView& ConstIterator3D::operator*() const { return msv; }


private:
  /** Current position. */
  ConstMatrixView msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor3:
class Tensor3;


/** A constant view of a Tensor3.

    This, together with the derived class Tensor3View, contains the
    main implementation of a Tensor3. It defines the concepts of
    Tensor3View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor3View from a subrange of
    a Tensor3View.

    The three dimensions of the tensor are called: page, row, column.

    The class Tensor3 is just a special case of a Tensor3View
    which also allocates storage. */
class ConstTensor3View {
public:
  // Member functions:
  /** Returns the number of pages. */
  Index npages() const { return mpr.mextent; }

  /** Returns the number of rows. */
  Index nrows() const { return mrr.mextent; }

  /** Returns the number of columns. */
  Index ncols() const { return mcr.mextent; }

  // Const index operators:
  ConstTensor3View operator()( const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( Index p,        Index r,        const Range& c ) const;
  ConstVectorView  operator()( Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( const Range& p, Index r,        Index c        ) const;

  Numeric          operator()( Index p,        Index r,        Index c) const;

  // Functions returning iterators:
  ConstIterator3D begin() const;
  ConstIterator3D end() const;
  
  // Friends:
  friend class Tensor3View;
  friend class ConstIterator4D;
  friend class ConstTensor4View;
  friend class ConstTensor5View;
  friend class ConstTensor6View;
  friend class ConstTensor7View;

  // Special constructor to make a Tensor3 view of a matrix.
  ConstTensor3View(const ConstMatrixView& a);

protected:
  // Constructors:
  ConstTensor3View();
  ConstTensor3View(Numeric *data,
                   const Range& p, const Range& r, const Range& c);
  ConstTensor3View(Numeric *data,
                   const Range& pp, const Range& pr, const Range& pc,
                   const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The page range of mdata that is actually used. */
  Range mpr;
  /** The row range of mdata that is actually used. */
  Range mrr;
  /** The column range of mdata that is actually used. */
  Range mcr;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The Tensor3View class

    This contains the main implementation of a Tensor3. It defines
    the concepts of Tensor3View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor3View
    from a subrange of a Tensor3View. 

    The class Tensor3 is just a special case of a Tensor3View
    which also allocates storage. */
class Tensor3View : public ConstTensor3View {
public:

  // Const index operators:
  ConstTensor3View operator()( const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( Index p,        Index r,        const Range& c ) const;
  ConstVectorView  operator()( Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( const Range& p, Index r,        Index c        ) const;

  Numeric          operator()( Index p,        Index r,        Index c        ) const;

  // Non-const index operators:

  Tensor3View operator()( const Range& p, const Range& r, const Range& c );

  MatrixView  operator()( const Range& p, const Range& r, Index c        );
  MatrixView  operator()( const Range& p, Index r,        const Range& c );
  MatrixView  operator()( Index p,        const Range& r, const Range& c );

  VectorView  operator()( Index p,        Index r,        const Range& c );
  VectorView  operator()( Index p,        const Range& r, Index c        );
  VectorView  operator()( const Range& p, Index r,        Index c        );

  /** Plain non-const index operator. */
  Numeric&    operator()( Index p,        Index r,        Index c)
    {
      // Check if indices are valid:
      assert( 0<=p );
      assert( 0<=r );
      assert( 0<=c );
      assert( p<mpr.mextent );
      assert( r<mrr.mextent );
      assert( c<mcr.mextent );

      return *( mdata +
                mpr.mstart + p*mpr.mstride +
                mrr.mstart + r*mrr.mstride +
                mcr.mstart + c*mcr.mstride );
    }

  // Functions returning const iterators:
  ConstIterator3D begin() const;
  ConstIterator3D end() const;
  // Functions returning iterators:
  Iterator3D begin();
  Iterator3D end();
  
  // Assignment operators:
  Tensor3View& operator=(const ConstTensor3View& v);
  Tensor3View& operator=(const Tensor3View& v);
  Tensor3View& operator=(const Tensor3& v);
  Tensor3View& operator=(Numeric x);

  // Other operators:
  Tensor3View& operator*=(Numeric x);
  Tensor3View& operator/=(Numeric x);
  Tensor3View& operator+=(Numeric x);
  Tensor3View& operator-=(Numeric x);

  Tensor3View& operator*=(const ConstTensor3View& x);
  Tensor3View& operator/=(const ConstTensor3View& x);
  Tensor3View& operator+=(const ConstTensor3View& x);
  Tensor3View& operator-=(const ConstTensor3View& x);

  // Friends:
  friend class Iterator4D;
  friend class Tensor4View;
  friend class Tensor5View;
  friend class Tensor6View;
  friend class Tensor7View;

  // Special constructor to make a Tensor3 view of a matrix.
  Tensor3View(const MatrixView& a);

protected:
  // Constructors:
  Tensor3View();
  Tensor3View(Numeric *data, const Range& p, const Range& r, const Range& c);
  Tensor3View(Numeric *data,
	      const Range& pp, const Range& pr, const Range& pc,
	      const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor3 class. This is a Tensor3View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor3View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor3 : public Tensor3View {
public:
  // Constructors:
  Tensor3();
  Tensor3(Index p, Index r, Index c);
  Tensor3(Index p, Index r, Index c, Numeric fill);
  Tensor3(const ConstTensor3View& v);
  Tensor3(const Tensor3& v);

  // Assignment operators:
  Tensor3& operator=(const Tensor3& x);
  Tensor3& operator=(Numeric x);

  // Resize function:
  void resize(Index p, Index r, Index c);

  // Destructor:
  ~Tensor3();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator3D origin,
                 const ConstIterator3D& end,
                 Iterator3D target);

inline void copy(Numeric x,
                 Iterator3D target,
                 const Iterator3D& end);

void transform( Tensor3View y,
                double (&my_func)(double),
                ConstTensor3View x );

Numeric max(const ConstTensor3View& x);

Numeric min(const ConstTensor3View& x);

std::ostream& operator<<(std::ostream& os, const ConstTensor3View& v);


// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor3. */
typedef Array<Tensor3> ArrayOfTensor3;

typedef Array<ArrayOfTensor3> ArrayOfArrayOfTensor3;


#endif    // matpackIII_h
