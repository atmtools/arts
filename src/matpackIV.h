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
  Implementation of Tensors of Rank 4.

  Based on Tensor3 by Stefan Buehler.

  The four dimensions are called: book, page, row, column.

  \author Wolfram-Andre Haas
  \date   2002-03-01
 */

#ifndef matpackIV_h
#define matpackIV_h

#include <iomanip>
#include "matpackI.h"
#include "matpackIII.h"

/** The outermost iterator class for rank 4 tensors. This takes into
    account the defined strided. */
class Iterator4D {
public:
  // Constructors:
  Iterator4D();
  Iterator4D(const Iterator4D& o);
  Iterator4D(const Tensor3View& x, Index stride);

  // Operators:
  Iterator4D& operator++();
  bool operator!=(const Iterator4D& other) const;
  Tensor3View* const operator->();
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
  ConstIterator4D();
  ConstIterator4D(const ConstIterator4D& o);
  ConstIterator4D(const ConstTensor3View& x, Index stride);

  // Operators:
  ConstIterator4D& operator++();
  bool operator!=(const ConstIterator4D& other) const;
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
  Index nbooks() const;
  Index npages() const;
  Index nrows()  const;
  Index ncols()  const;

  // Const index operators:
  ConstTensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& b, Index p,        const Range &r, const Range& c ) const;
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

  Numeric          operator()( Index b,        Index p,        Index r,        Index c        ) const;

  // Functions returning iterators:
  ConstIterator4D begin() const;
  ConstIterator4D end()   const;

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
  ConstTensor3View operator()( const Range& b, Index p,        const Range &r, const Range& c ) const;
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

  Numeric          operator()( Index b,        Index p,        Index r,        Index c        ) const;

  // Non-const index operators:

  Tensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c );

  Tensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        );
  Tensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c );
  Tensor3View operator()( const Range& b, Index p,        const Range &r, const Range& c );
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

  Numeric&    operator()( Index b,        Index p,        Index r,        Index c        );

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
  Tensor4& operator=(const Tensor4& x);
  Tensor4& operator=(Numeric x);

  // Resize function:
  void resize(Index b, Index p, Index r, Index c);

  // Destructor:
  ~Tensor4();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator4D origin,
                 const ConstIterator4D& end,
                 Iterator4D target);

inline void copy(Numeric x,
                 Iterator4D target,
                 const Iterator4D& end);

void transform( Tensor4View y,
                double (&my_func)(double),
                ConstTensor4View x );

Numeric max(const ConstTensor4View& x);

Numeric min(const ConstTensor4View& x);

std::ostream& operator<<(std::ostream& os, const ConstTensor4View& v);


// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor4. */
typedef Array<Tensor4> ArrayOfTensor4;


#endif    // matpackIV_h
