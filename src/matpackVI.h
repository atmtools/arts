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
  Implementation of Tensors of Rank 6.

  Dimensions are called: vitrine, shelf, book, page, row, column.
  or short:              v,       s,     b,    p,    r,   c
  
  \author Stefan Buehler
  \date   2001-11-22
 */

#ifndef matpackVI_h
#define matpackVI_h

#include <iomanip>
#include "matpackV.h"

/** The outermost iterator class for rank 6 tensors. This takes into
    account the defined strided. */
class Iterator6D {
public:
  // Constructors:
  Iterator6D();
  Iterator6D(const Iterator6D& o);
  Iterator6D(const Tensor5View& x, Index stride);

  // Operators:
  Iterator6D& operator++();
  bool operator!=(const Iterator6D& other) const;
  Tensor5View* const operator->();
  Tensor5View& operator*();
  
private:
  /** Current position. */
  Tensor5View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator6D. */
class ConstIterator6D {
public:
  // Constructors:
  ConstIterator6D();
  ConstIterator6D(const ConstIterator6D& o);
  ConstIterator6D(const ConstTensor5View& x, Index stride);

  // Operators:
  ConstIterator6D& operator++();
  bool operator!=(const ConstIterator6D& other) const;
  const ConstTensor5View* operator->() const;
  const ConstTensor5View& operator*() const;

private:
  /** Current position. */
  ConstTensor5View msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor6:
class Tensor6;


/** A constant view of a Tensor6.

    This, together with the derived class Tensor6View, contains the
    main implementation of a Tensor6. It defines the concepts of
    Tensor6View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor6View from a subrange of
    a Tensor6View.

    Dimensions are called: vitrine, shelf, book, page, row, column.
    or short:              v,       s,     b,    p,    r,   c

    The class Tensor6 is just a special case of a Tensor6View
    which also allocates storage. */
class ConstTensor6View {
public:
  // Member functions:
  Index nvitrines() const;
  Index nshelves()  const;
  Index nbooks()    const;
  Index npages()    const;
  Index nrows()     const;
  Index ncols()     const;

  // Const index operators:

  // Result 6D (1 combination)
  // ------
  ConstTensor6View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 5D (6 combinations)
  // -----|
  ConstTensor5View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ----|-
  ConstTensor5View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ---|--
  ConstTensor5View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // --|---
  ConstTensor5View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|----
  ConstTensor5View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-----
  ConstTensor5View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 4D (5+4+3+2+1 = 15 combinations)
  // ----||
  ConstTensor4View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ---|-|
  ConstTensor4View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // --|--|
  ConstTensor4View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|---|
  ConstTensor4View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |----|
  ConstTensor4View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ---||-
  ConstTensor4View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // --|-|-
  ConstTensor4View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|--|-
  ConstTensor4View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |---|-
  ConstTensor4View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // --||--
  ConstTensor4View operator()( const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -|-|--
  ConstTensor4View operator()( const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |--|--
  ConstTensor4View operator()( Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -||---
  ConstTensor4View operator()( const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-|---
  ConstTensor4View operator()( Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||----
  ConstTensor4View operator()( Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
  // ---|||
  ConstTensor3View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --|-||
  ConstTensor3View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|--||
  ConstTensor3View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |---||
  ConstTensor3View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // --||-|
  ConstTensor3View operator()( const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -|-|-|
  ConstTensor3View operator()( const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |--|-|
  ConstTensor3View operator()( Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -||--|
  ConstTensor3View operator()( const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |-|--|
  ConstTensor3View operator()( Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // ||---|
  ConstTensor3View operator()( Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // --|||-
  ConstTensor3View operator()( const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|-||-
  ConstTensor3View operator()( const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |--||-
  ConstTensor3View operator()( Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -||-|-
  ConstTensor3View operator()( const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |-|-|-
  ConstTensor3View operator()( Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // ||--|-
  ConstTensor3View operator()( Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // -|||--
  ConstTensor3View operator()( const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |-||--
  ConstTensor3View operator()( Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // ||-|--
  ConstTensor3View operator()( Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |||---
  ConstTensor3View operator()( Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 2D (15 combinations)
  // IIII--
  ConstMatrixView  operator()( Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // III-I-
  ConstMatrixView  operator()( Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // II-II-
  ConstMatrixView  operator()( Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // I-III-
  ConstMatrixView  operator()( Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -IIII-
  ConstMatrixView  operator()( const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // III--I
  ConstMatrixView  operator()( Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // II-I-I
  ConstMatrixView  operator()( Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // I-II-I
  ConstMatrixView  operator()( Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -III-I
  ConstMatrixView  operator()( const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // II--II
  ConstMatrixView  operator()( Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // I-I-II
  ConstMatrixView  operator()( Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -II-II
  ConstMatrixView  operator()( const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // I--III
  ConstMatrixView  operator()( Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -I-III
  ConstMatrixView  operator()( const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --IIII
  ConstMatrixView  operator()( const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 1D (6 combinations)
  // IIIII-
  ConstVectorView  operator()( Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // IIII-I
  ConstVectorView  operator()( Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // III-II
  ConstVectorView  operator()( Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // II-III
  ConstVectorView  operator()( Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // I-IIII
  ConstVectorView  operator()( Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -IIIII
  ConstVectorView  operator()( const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result scalar (1 combination)
  // IIIIII
  Numeric          operator()( Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Functions returning iterators:
  ConstIterator6D begin() const;
  ConstIterator6D end() const;
  
  // Friends:
  friend class ConstIterator7D;
  friend class Tensor6View;
  friend class ConstTensor7View;

protected:
  // Constructors:
  ConstTensor6View();
  ConstTensor6View(Numeric *data,
		   const Range& v, const Range& s, const Range& b,
		   const Range& p, const Range& r, const Range& c);
  ConstTensor6View(Numeric *data,
		   const Range& pv, const Range& ps, const Range& pb,
		   const Range& pp, const Range& pr, const Range& pc,
		   const Range& nv, const Range& ns, const Range& nb,
		   const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The vitrine range of mdata that is actually used. */
  Range mvr;
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

/** The Tensor6View class

    This contains the main implementation of a Tensor6. It defines
    the concepts of Tensor6View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor6View
    from a subrange of a Tensor6View. 

    The class Tensor6 is just a special case of a Tensor6View
    which also allocates storage. */
class Tensor6View : public ConstTensor6View {
public:

  // Const index operators:

  // Result 6D (1 combination)
  // ------
  ConstTensor6View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 5D (6 combinations)
  // -----|
  ConstTensor5View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ----|-
  ConstTensor5View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ---|--
  ConstTensor5View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // --|---
  ConstTensor5View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|----
  ConstTensor5View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-----
  ConstTensor5View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 4D (5+4+3+2+1 = 15 combinations)
  // ----||
  ConstTensor4View operator()( const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ---|-|
  ConstTensor4View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // --|--|
  ConstTensor4View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|---|
  ConstTensor4View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |----|
  ConstTensor4View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ---||-
  ConstTensor4View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // --|-|-
  ConstTensor4View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|--|-
  ConstTensor4View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |---|-
  ConstTensor4View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // --||--
  ConstTensor4View operator()( const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -|-|--
  ConstTensor4View operator()( const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |--|--
  ConstTensor4View operator()( Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -||---
  ConstTensor4View operator()( const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-|---
  ConstTensor4View operator()( Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||----
  ConstTensor4View operator()( Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
  // ---|||
  ConstTensor3View operator()( const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --|-||
  ConstTensor3View operator()( const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|--||
  ConstTensor3View operator()( const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |---||
  ConstTensor3View operator()( Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // --||-|
  ConstTensor3View operator()( const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -|-|-|
  ConstTensor3View operator()( const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |--|-|
  ConstTensor3View operator()( Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -||--|
  ConstTensor3View operator()( const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |-|--|
  ConstTensor3View operator()( Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // ||---|
  ConstTensor3View operator()( Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // --|||-
  ConstTensor3View operator()( const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|-||-
  ConstTensor3View operator()( const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |--||-
  ConstTensor3View operator()( Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -||-|-
  ConstTensor3View operator()( const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |-|-|-
  ConstTensor3View operator()( Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // ||--|-
  ConstTensor3View operator()( Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // -|||--
  ConstTensor3View operator()( const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |-||--
  ConstTensor3View operator()( Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // ||-|--
  ConstTensor3View operator()( Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |||---
  ConstTensor3View operator()( Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 2D (15 combinations)
  // IIII--
  ConstMatrixView  operator()( Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // III-I-
  ConstMatrixView  operator()( Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // II-II-
  ConstMatrixView  operator()( Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // I-III-
  ConstMatrixView  operator()( Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -IIII-
  ConstMatrixView  operator()( const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // III--I
  ConstMatrixView  operator()( Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // II-I-I
  ConstMatrixView  operator()( Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // I-II-I
  ConstMatrixView  operator()( Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -III-I
  ConstMatrixView  operator()( const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // II--II
  ConstMatrixView  operator()( Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // I-I-II
  ConstMatrixView  operator()( Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -II-II
  ConstMatrixView  operator()( const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // I--III
  ConstMatrixView  operator()( Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -I-III
  ConstMatrixView  operator()( const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --IIII
  ConstMatrixView  operator()( const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 1D (6 combinations)
  // IIIII-
  ConstVectorView  operator()( Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // IIII-I
  ConstVectorView  operator()( Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // III-II
  ConstVectorView  operator()( Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // II-III
  ConstVectorView  operator()( Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // I-IIII
  ConstVectorView  operator()( Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -IIIII
  ConstVectorView  operator()( const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result scalar (1 combination)
  // IIIIII
  Numeric          operator()( Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Non-const index operators:

  // Result 6D (1 combination)
  // ------
  Tensor6View operator()( const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 5D (6 combinations)
  // -----|
  Tensor5View operator()( const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ----|-
  Tensor5View operator()( const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ---|--
  Tensor5View operator()( const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // --|---
  Tensor5View operator()( const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // -|----
  Tensor5View operator()( const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // |-----
  Tensor5View operator()( Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 4D (5+4+3+2+1 = 15 combinations)
  // ----||
  Tensor4View operator()( const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ---|-|
  Tensor4View operator()( const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // --|--|
  Tensor4View operator()( const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -|---|
  Tensor4View operator()( const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // |----|
  Tensor4View operator()( Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ---||-
  Tensor4View operator()( const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // --|-|-
  Tensor4View operator()( const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -|--|-
  Tensor4View operator()( const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // |---|-
  Tensor4View operator()( Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // --||--
  Tensor4View operator()( const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -|-|--
  Tensor4View operator()( const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // |--|--
  Tensor4View operator()( Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // -||---
  Tensor4View operator()( const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |-|---
  Tensor4View operator()( Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // ||----
  Tensor4View operator()( Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
  // ---|||
  Tensor3View operator()( const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // --|-||
  Tensor3View operator()( const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // -|--||
  Tensor3View operator()( const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |---||
  Tensor3View operator()( Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // --||-|
  Tensor3View operator()( const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // -|-|-|
  Tensor3View operator()( const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |--|-|
  Tensor3View operator()( Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // -||--|
  Tensor3View operator()( const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |-|--|
  Tensor3View operator()( Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // ||---|
  Tensor3View operator()( Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // --|||-
  Tensor3View operator()( const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // -|-||-
  Tensor3View operator()( const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |--||-
  Tensor3View operator()( Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // -||-|-
  Tensor3View operator()( const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |-|-|-
  Tensor3View operator()( Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // ||--|-
  Tensor3View operator()( Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // -|||--
  Tensor3View operator()( const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |-||--
  Tensor3View operator()( Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // ||-|--
  Tensor3View operator()( Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // |||---
  Tensor3View operator()( Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);

  // Result 2D (15 combinations)
  // IIII--
  MatrixView  operator()( Index        v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // III-I-
  MatrixView  operator()( Index        v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // II-II-
  MatrixView  operator()( Index        v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // I-III-
  MatrixView  operator()( Index        v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // -IIII-
  MatrixView  operator()( const Range& v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // III--I
  MatrixView  operator()( Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // II-I-I
  MatrixView  operator()( Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // I-II-I
  MatrixView  operator()( Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // -III-I
  MatrixView  operator()( const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // II--II
  MatrixView  operator()( Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // I-I-II
  MatrixView  operator()( Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // -II-II
  MatrixView  operator()( const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // I--III
  MatrixView  operator()( Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // -I-III
  MatrixView  operator()( const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // --IIII
  MatrixView  operator()( const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result 1D (6 combinations)
  // IIIII-
  VectorView  operator()( Index        v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // IIII-I
  VectorView  operator()( Index        v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // III-II
  VectorView  operator()( Index        v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // II-III
  VectorView  operator()( Index        v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // I-IIII
  VectorView  operator()( Index        v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // -IIIII
  VectorView  operator()( const Range& v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result scalar (1 combination)
  // IIIIII
  Numeric&          operator()( Index        v, Index        s, Index        b,
				Index        p, Index        r, Index        c);


  // Functions returning const iterators:
  ConstIterator6D begin() const;
  ConstIterator6D end() const;
  // Functions returning iterators:
  Iterator6D begin();
  Iterator6D end();
  
  // Assignment operators:
  Tensor6View& operator=(const ConstTensor6View& v);
  Tensor6View& operator=(const Tensor6View& v);
  Tensor6View& operator=(const Tensor6& v);
  Tensor6View& operator=(Numeric x);

  // Other operators:
  Tensor6View& operator*=(Numeric x);
  Tensor6View& operator/=(Numeric x);
  Tensor6View& operator+=(Numeric x);
  Tensor6View& operator-=(Numeric x);

  Tensor6View& operator*=(const ConstTensor6View& x);
  Tensor6View& operator/=(const ConstTensor6View& x);
  Tensor6View& operator+=(const ConstTensor6View& x);
  Tensor6View& operator-=(const ConstTensor6View& x);

  // Friends:
  friend class Iterator7D;
  friend class Tensor7View;

protected:
  // Constructors:
  Tensor6View();
  Tensor6View(Numeric *data,
	      const Range& v, const Range& s, const Range& b,
	      const Range& p, const Range& r, const Range& c);
  Tensor6View(Numeric *data,
	      const Range& pv, const Range& ps, const Range& pb,
	      const Range& pp, const Range& pr, const Range& pc,
	      const Range& nv, const Range& ns, const Range& nb,
	      const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor6 class. This is a Tensor6View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor6View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor6 : public Tensor6View {
public:
  // Constructors:
  Tensor6();
  Tensor6(Index        v, Index        s, Index        b,
	  Index        p, Index        r, Index        c);
  Tensor6(Index        v, Index        s, Index        b,
	  Index        p, Index        r, Index        c,
	  Numeric fill);
  Tensor6(const ConstTensor6View& v);
  Tensor6(const Tensor6& v);

  // Assignment operators:
  Tensor6& operator=(const Tensor6& x);
  Tensor6& operator=(Numeric x);

  // Resize function:
  void resize(Index        v, Index        s, Index        b,
	      Index        p, Index        r, Index        c);

  // Destructor:
  ~Tensor6();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator6D origin,
                 const ConstIterator6D& end,
                 Iterator6D target);

inline void copy(Numeric x,
                 Iterator6D target,
                 const Iterator6D& end);



// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor6. */
typedef Array<Tensor6> ArrayOfTensor6;



// Functions for Iterator6D
// ------------------------

/** Default constructor. */
inline Iterator6D::Iterator6D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator6D::Iterator6D(const Iterator6D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator6D::Iterator6D(const Tensor5View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline Iterator6D& Iterator6D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy.
    FIXME: Is it really necessary to have such a complicated check
    here? It could be sufficient to just test
    msv.mdata!=other.msv.mdata. */
inline bool Iterator6D::operator!=(const Iterator6D& other) const
{
  if ( msv.mdata +
       msv.msr.mstart +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart 
       !=
       other.msv.mdata +
       other.msv.msr.mstart +
       other.msv.mbr.mstart +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline Tensor5View* const Iterator6D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline Tensor5View& Iterator6D::operator*()
{
  return msv;
}

// Functions for ConstIterator6D
// -----------------------------

/** Default constructor. */
inline ConstIterator6D::ConstIterator6D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator6D::ConstIterator6D(const ConstIterator6D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator6D::ConstIterator6D(const ConstTensor5View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline ConstIterator6D& ConstIterator6D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy.
    FIXME: Is it really necessaary to have such a complicated check
    here? It could be sufficient to just test
    msv.mdata!=other.msv.mdata. */
inline bool ConstIterator6D::operator!=(const ConstIterator6D& other) const
{
  if ( msv.mdata +
       msv.msr.mstart +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.msr.mstart +
       other.msv.mbr.mstart +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
inline const ConstTensor5View* ConstIterator6D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstTensor5View& ConstIterator6D::operator*() const
{
  return msv;
}



// Functions for ConstTensor6View:
// ------------------------------

/** Returns the number of vitrines. */
inline Index ConstTensor6View::nvitrines() const
{
  return mvr.mextent;
}

/** Returns the number of shelves. */
inline Index ConstTensor6View::nshelves() const
{
  return msr.mextent;
}

/** Returns the number of books. */
inline Index ConstTensor6View::nbooks() const
{
  return mbr.mextent;
}

/** Returns the number of pages. */
inline Index ConstTensor6View::npages() const
{
  return mpr.mextent;
}

/** Returns the number of rows. */
inline Index ConstTensor6View::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor6View::ncols() const
{
  return mcr.mextent;
}

// Const index operators:

// Result 6D (1 combination)
// ------
inline ConstTensor6View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View(mdata,
			  mvr, msr, mbr, mpr, mrr, mcr,
			  v,   s,   b,   p,   r,   c);
}

#define CHECK(x)   assert( 0 <= x ); assert( x < m ## x ## r.mextent )
#define OFFSET(x)  m ## x ## r.mstart + x * m ## x ## r.mstride

// Result 5D (6 combinations)
// -----|
inline ConstTensor5View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(c),
                          mvr, msr, mbr, mpr, mrr,
                          v,   s,   b,   p,   r);
}

// ----|-
inline ConstTensor5View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(r);
  return ConstTensor5View(mdata + OFFSET(r),
                          mvr, msr, mbr, mpr, mcr,
                          v,   s,   b,   p,   c);
}

// ---|--
inline ConstTensor5View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(p);
  return ConstTensor5View(mdata + OFFSET(p),
                          mvr, msr, mbr, mrr, mcr,
                          v,   s,   b,   r,   c);
}

// --|---
inline ConstTensor5View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(b);
  return ConstTensor5View(mdata + OFFSET(b),
                          mvr, msr, mpr, mrr, mcr,
                          v,   s,   p,   r,   c);
}

// -|----
inline ConstTensor5View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(s);
  return ConstTensor5View(mdata + OFFSET(s),
                          mvr, mbr, mpr, mrr, mcr,
                          v,   b,   p,   r,   c);
}

// |-----
inline ConstTensor5View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  return ConstTensor5View(mdata + OFFSET(v),
                          msr, mbr, mpr, mrr, mcr,
                          s,   b,   p,   r,   c);
}


// Result 4D (5+4+3+2+1 = 15 combinations)
// ----||
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(r) + OFFSET(c),
                          mvr, msr, mbr, mpr,
                          v,   s,   b,   p   );
}

// ---|-|
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(p) + OFFSET(c),
                          mvr, msr, mbr, mrr,
                          v,   s,   b,   r    );
}

// --|--|
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(c),
                          mvr, msr, mpr, mrr,
                          v,   s,   p,   r    );
}

// -|---|
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(c),
                          mvr, mbr, mpr, mrr,
                          v,   b,   p,   r    );
}

// |----|
inline ConstTensor4View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(c),
                          msr, mbr, mpr, mrr,
                          s,   b,   p,   r    );
}

// ---||-
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(p) + OFFSET(r),
                          mvr, msr, mbr, mcr,
                          v,   s,   b,   c    );
}

// --|-|-
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(r),
                          mvr, msr, mpr, mcr,
                          v,   s,   p,   c    );
}

// -|--|-
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(r),
                          mvr, mbr, mpr, mcr,
                          v,   b,   p,   c    );
}

// |---|-
inline ConstTensor4View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(r),
                          msr, mbr, mpr, mcr,
                          s,   b,   p,   c    );
}

// --||--
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(p),
                          mvr, msr, mrr, mcr,
                          v,   s,   r,   c    );
}

// -|-|--
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(s);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(p),
                          mvr, mbr, mrr, mcr,
                          v,   b,   r,   c    );
}

// |--|--
inline ConstTensor4View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(p),
                          msr, mbr, mrr, mcr,
                          s,   b,   r,   c    );
}

// -||---
inline ConstTensor4View ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b),
                          mvr, mpr, mrr, mcr,
                          v,   p,   r,   c    );
}

// |-|---
inline ConstTensor4View ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b),
                          msr, mpr, mrr, mcr,
                          s,   p,   r,   c    );
}

// ||----
inline ConstTensor4View ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s),
                          mbr, mpr, mrr, mcr,
                          b,   p,   r,   c    );
}


// Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
// ---|||
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mbr, 
                          v,   s,   b     );
}

// --|-||
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mpr, 
                          v,   s,   p     );
}

// -|--||
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mvr, mbr, mpr, 
                          v,   b,   p     );
}

// |---||
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(r) + OFFSET(c),
                          msr, mbr, mpr, 
                          s,   b,   p     );
}

// --||-|
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mvr, msr, mrr, 
                          v,   s,   r     );
}

// -|-|-|
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mvr, mbr, mrr, 
                          v,   b,   r     );
}

// |--|-|
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(c),
                          msr, mbr, mrr, 
                          s,   b,   r     );
}

// -||--|
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mvr, mpr, mrr, 
                          v,   p,   r     );
}

// |-|--|
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(c),
                          msr, mpr, mrr, 
                          s,   p,   r     );
}

// ||---|
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(c),
                          mbr, mpr, mrr, 
                          b,   p,   r     );
}

// --|||-
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mvr, msr, mcr, 
                          v,   s,   c     );
}

// -|-||-
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mvr, mbr, mcr, 
                          v,   b,   c     );
}

// |--||-
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r),
                          msr, mbr, mcr, 
                          s,   b,   c     );
}

// -||-|-
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mvr, mpr, mcr, 
                          v,   p,   c     );
}

// |-|-|-
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r),
                          msr, mpr, mcr, 
                          s,   p,   c     );
}

// ||--|-
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r),
                          mbr, mpr, mcr, 
                          b,   p,   c     );
}

// -|||--
inline ConstTensor3View ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mvr, mrr, mcr, 
                          v,   r,   c     );
}

// |-||--
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p),
                          msr, mrr, mcr, 
                          s,   r,   c     );
}

// ||-|--
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p),
                          mbr, mrr, mcr, 
                          b,   r,   c     );
}

// |||---
inline ConstTensor3View ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b),
                          mpr, mrr, mcr, 
                          p,   r,   c     );
}


// Result 2D (15 combinations)
// IIII--
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mrr, mcr, 
                          r,   c     );
}

// III-I-
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mpr, mcr, 
                          p,   c     );
}

// II-II-
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mbr, mcr, 
                          b,   c     );
}

// I-III-
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          msr, mcr, 
                          s,   c     );
}

// -IIII-
inline ConstMatrixView  ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mvr, mcr, 
                          v,   c     );
}

// III--I
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mpr, mrr, 
                          p,   r     );
}

// II-I-I
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mbr, mrr, 
                          b,   r     );
}

// I-II-I
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          msr, mrr, 
                          s,   r     );
}

// -III-I
inline ConstMatrixView  ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mvr, mrr, 
                          v,   r     );
}

// II--II
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mbr, mpr, 
                          b,   p     );
}

// I-I-II
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          msr, mpr, 
                          s,   p     );
}

// -II-II
inline ConstMatrixView  ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mvr, mpr, 
                          v,   p     );
}

// I--III
inline ConstMatrixView  ConstTensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          msr, mbr, 
                          s,   b     );
}

// -I-III
inline ConstMatrixView  ConstTensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, mbr, 
                          v,   b     );
}

// --IIII
inline ConstMatrixView  ConstTensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, msr, 
                          v,   s     );
}


// Result 1D (6 combinations)
// IIIII-
inline ConstVectorView  ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstVectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mcr,  
                          c      );
}

// IIII-I
inline ConstVectorView  ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mrr,  
                          r      );
}

// III-II
inline ConstVectorView  ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mpr,  
                          p      );
}

// II-III
inline ConstVectorView  ConstTensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mbr,  
                          b      );
}

// I-IIII
inline ConstVectorView  ConstTensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          msr,  
                          s      );
}

// -IIIII
inline ConstVectorView  ConstTensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr,  
                          v      );
}


// Result scalar (1 combination)
// IIIIII
inline Numeric          ConstTensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return                *(mdata +
			  OFFSET(v) + OFFSET(s) + OFFSET(b) +
			  OFFSET(p) + OFFSET(r) + OFFSET(c)    );
}



/** Return const iterator to first sub-tensor. */
inline ConstIterator6D ConstTensor6View::begin() const
{
  return ConstIterator6D( ConstTensor5View(mdata+mvr.mstart,
					  msr,
					  mbr,
					  mpr,
					  mrr,
					  mcr),
			  mvr.mstride);
}

/** Return const iterator behind last sub-tensor. */
inline ConstIterator6D ConstTensor6View::end() const
{
  return ConstIterator6D( ConstTensor5View(mdata + mvr.mstart +
                                          (mvr.mextent)*mvr.mstride,
					  msr,
					  mbr,
					  mpr,
                                          mrr,
                                          mcr),
                          mvr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor6View::ConstTensor6View() :
  mvr(0,0,1), msr(0,0,1), mbr(0,0,1), 
  mpr(0,0,1), mrr(0,0,1), mcr(0,0,1),
  mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor6 to initialize
    its own Tensor6View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
inline ConstTensor6View::ConstTensor6View(Numeric *data,
                                          const Range& v,
                                          const Range& s,
                                          const Range& b,
                                          const Range& p,
                                          const Range& r,
                                          const Range& c) :
  mvr(v),
  msr(s),
  mbr(b),
  mpr(p),
  mrr(r),
  mcr(c),
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
inline ConstTensor6View::ConstTensor6View(Numeric *data,
                                          const Range& pv,
                                          const Range& ps,
                                          const Range& pb,
                                          const Range& pp,
                                          const Range& pr,
                                          const Range& pc,
                                          const Range& nv,
                                          const Range& ns,
                                          const Range& nb,
                                          const Range& np,
                                          const Range& nr,
                                          const Range& nc) :
  mvr(pv,nv),
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
    Tensor5 to print each page in turn. */
inline std::ostream& operator<<(std::ostream& os, const ConstTensor6View& v)
{
  // Page iterators:
  ConstIterator6D ip=v.begin();
  const ConstIterator6D end_page=v.end();

  if ( ip!=end_page )
    {
      os << *ip;
      ++ip;
    }

  for ( ; ip!=end_page; ++ip )
    {
      os << "\n\n";
      os << *ip;
    }

  return os;
}


// Functions for Tensor6View:
// -------------------------

// Const index operators:

// Result 6D (1 combination)
// ------
inline ConstTensor6View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Result 5D (6 combinations)
// -----|
inline ConstTensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ----|-
inline ConstTensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ---|--
inline ConstTensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --|---
inline ConstTensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|----
inline ConstTensor5View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |-----
inline ConstTensor5View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Result 4D (5+4+3+2+1 = 15 combinations)
// ----||
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ---|-|
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --|--|
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|---|
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |----|
inline ConstTensor4View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ---||-
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --|-|-
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|--|-
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |---|-
inline ConstTensor4View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --||--
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|-|--
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |--|--
inline ConstTensor4View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -||---
inline ConstTensor4View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |-|---
inline ConstTensor4View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ||----
inline ConstTensor4View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
// ---|||
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --|-||
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|--||
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |---||
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --||-|
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|-|-|
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |--|-|
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -||--|
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |-|--|
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ||---|
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --|||-
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|-||-
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |--||-
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -||-|-
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |-|-|-
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ||--|-
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -|||--
inline ConstTensor3View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |-||--
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// ||-|--
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// |||---
inline ConstTensor3View Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Result 2D (15 combinations)
// IIII--
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// III-I-
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// II-II-
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// I-III-
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -IIII-
inline ConstMatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// III--I
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// II-I-I
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// I-II-I
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -III-I
inline ConstMatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// II--II
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// I-I-II
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -II-II
inline ConstMatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// I--III
inline ConstMatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -I-III
inline ConstMatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// --IIII
inline ConstMatrixView  Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Result 1D (6 combinations)
// IIIII-
inline ConstVectorView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// IIII-I
inline ConstVectorView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// III-II
inline ConstVectorView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// II-III
inline ConstVectorView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// I-IIII
inline ConstVectorView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}

// -IIIII
inline ConstVectorView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Result scalar (1 combination)
// IIIIII
inline Numeric          Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor6View::operator()(v,s,b,p,r,c);  
}


// Non-const index operators:

// Result 6D (1 combination)
// ------
inline Tensor6View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  return Tensor6View(mdata,
			  mvr, msr, mbr, mpr, mrr, mcr,
			  v,   s,   b,   p,   r,   c);
}


// Result 5D (6 combinations)
// -----|
inline Tensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(c);
  return Tensor5View(mdata + OFFSET(c),
                          mvr, msr, mbr, mpr, mrr,
                          v,   s,   b,   p,   r);
}

// ----|-
inline Tensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(r);
  return Tensor5View(mdata + OFFSET(r),
                          mvr, msr, mbr, mpr, mcr,
                          v,   s,   b,   p,   c);
}

// ---|--
inline Tensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(p);
  return Tensor5View(mdata + OFFSET(p),
                          mvr, msr, mbr, mrr, mcr,
                          v,   s,   b,   r,   c);
}

// --|---
inline Tensor5View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(b);
  return Tensor5View(mdata + OFFSET(b),
                          mvr, msr, mpr, mrr, mcr,
                          v,   s,   p,   r,   c);
}

// -|----
inline Tensor5View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(s);
  return Tensor5View(mdata + OFFSET(s),
                          mvr, mbr, mpr, mrr, mcr,
                          v,   b,   p,   r,   c);
}

// |-----
inline Tensor5View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  return Tensor5View(mdata + OFFSET(v),
                          msr, mbr, mpr, mrr, mcr,
                          s,   b,   p,   r,   c);
}


// Result 4D (5+4+3+2+1 = 15 combinations)
// ----||
inline Tensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(r) + OFFSET(c),
                          mvr, msr, mbr, mpr,
                          v,   s,   b,   p    );
}

// ---|-|
inline Tensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(p) + OFFSET(c),
                          mvr, msr, mbr, mrr,
                          v,   s,   b,   r    );
}

// --|--|
inline Tensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(c),
                          mvr, msr, mpr, mrr,
                          v,   s,   p,   r    );
}

// -|---|
inline Tensor4View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(c),
                          mvr, mbr, mpr, mrr,
                          v,   b,   p,   r    );
}

// |----|
inline Tensor4View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(c),
                          msr, mbr, mpr, mrr,
                          s,   b,   p,   r    );
}

// ---||-
inline Tensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(p) + OFFSET(r),
                          mvr, msr, mbr, mcr,
                          v,   s,   b,   c    );
}

// --|-|-
inline Tensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(r),
                          mvr, msr, mpr, mcr,
                          v,   s,   p,   c    );
}

// -|--|-
inline Tensor4View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(r),
                          mvr, mbr, mpr, mcr,
                          v,   b,   p,   c    );
}

// |---|-
inline Tensor4View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(r),
                          msr, mbr, mpr, mcr,
                          s,   b,   p,   c    );
}

// --||--
inline Tensor4View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(p),
                          mvr, msr, mrr, mcr,
                          v,   s,   r,   c    );
}

// -|-|--
inline Tensor4View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(s);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(p),
                          mvr, mbr, mrr, mcr,
                          v,   b,   r,   c    );
}

// |--|--
inline Tensor4View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(p),
                          msr, mbr, mrr, mcr,
                          s,   b,   r,   c    );
}

// -||---
inline Tensor4View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b),
                          mvr, mpr, mrr, mcr,
                          v,   p,   r,   c    );
}

// |-|---
inline Tensor4View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b),
                          msr, mpr, mrr, mcr,
                          s,   p,   r,   c    );
}

// ||----
inline Tensor4View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s),
                          mbr, mpr, mrr, mcr,
                          b,   p,   r,   c    );
}


// Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
// ---|||
inline Tensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mbr, 
                          v,   s,   b     );
}

// --|-||
inline Tensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mpr, 
                          v,   s,   p     );
}

// -|--||
inline Tensor3View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mvr, mbr, mpr, 
                          v,   b,   p     );
}

// |---||
inline Tensor3View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(r) + OFFSET(c),
                          msr, mbr, mpr, 
                          s,   b,   p     );
}

// --||-|
inline Tensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mvr, msr, mrr, 
                          v,   s,   r     );
}

// -|-|-|
inline Tensor3View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mvr, mbr, mrr, 
                          v,   b,   r     );
}

// |--|-|
inline Tensor3View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(c),
                          msr, mbr, mrr, 
                          s,   b,   r     );
}

// -||--|
inline Tensor3View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mvr, mpr, mrr, 
                          v,   p,   r     );
}

// |-|--|
inline Tensor3View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(c),
                          msr, mpr, mrr, 
                          s,   p,   r     );
}

// ||---|
inline Tensor3View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(c),
                          mbr, mpr, mrr, 
                          b,   p,   r     );
}

// --|||-
inline Tensor3View Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mvr, msr, mcr, 
                          v,   s,   c     );
}

// -|-||-
inline Tensor3View Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mvr, mbr, mcr, 
                          v,   b,   c     );
}

// |--||-
inline Tensor3View Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r),
                          msr, mbr, mcr, 
                          s,   b,   c     );
}

// -||-|-
inline Tensor3View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mvr, mpr, mcr, 
                          v,   p,   c     );
}

// |-|-|-
inline Tensor3View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r),
                          msr, mpr, mcr, 
                          s,   p,   c     );
}

// ||--|-
inline Tensor3View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r),
                          mbr, mpr, mcr, 
                          b,   p,   c     );
}

// -|||--
inline Tensor3View Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mvr, mrr, mcr, 
                          v,   r,   c     );
}

// |-||--
inline Tensor3View Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p),
                          msr, mrr, mcr, 
                          s,   r,   c     );
}

// ||-|--
inline Tensor3View Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p),
                          mbr, mrr, mcr, 
                          b,   r,   c     );
}

// |||---
inline Tensor3View Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b),
                          mpr, mrr, mcr, 
                          p,   r,   c     );
}


// Result 2D (15 combinations)
// IIII--
inline MatrixView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mrr, mcr, 
                          r,   c     );
}

// III-I-
inline MatrixView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mpr, mcr, 
                          p,   c     );
}

// II-II-
inline MatrixView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mbr, mcr, 
                          b,   c     );
}

// I-III-
inline MatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          msr, mcr, 
                          s,   c     );
}

// -IIII-
inline MatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mvr, mcr, 
                          v,   c     );
}

// III--I
inline MatrixView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mpr, mrr, 
                          p,   r     );
}

// II-I-I
inline MatrixView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mbr, mrr, 
                          b,   r     );
}

// I-II-I
inline MatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          msr, mrr, 
                          s,   r     );
}

// -III-I
inline MatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mvr, mrr, 
                          v,   r     );
}

// II--II
inline MatrixView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mbr, mpr, 
                          b,   p     );
}

// I-I-II
inline MatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          msr, mpr, 
                          s,   p     );
}

// -II-II
inline MatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mvr, mpr, 
                          v,   p     );
}

// I--III
inline MatrixView  Tensor6View::operator()
  ( Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          msr, mbr, 
                          s,   b     );
}

// -I-III
inline MatrixView  Tensor6View::operator()
  ( const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, mbr, 
                          v,   b     );
}

// --IIII
inline MatrixView  Tensor6View::operator()
  ( const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, msr, 
                          v,   s     );
}


// Result 1D (6 combinations)
// IIIII-
inline VectorView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  VectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mcr,  
                          c      );
}

// IIII-I
inline VectorView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  VectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mrr,  
                          r      );
}

// III-II
inline VectorView  Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mpr,  
                          p      );
}

// II-III
inline VectorView  Tensor6View::operator()
  ( Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mbr,  
                          b      );
}

// I-IIII
inline VectorView  Tensor6View::operator()
  ( Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          msr,  
                          s      );
}

// -IIIII
inline VectorView  Tensor6View::operator()
  ( const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr,  
                          v      );
}


// Result scalar (1 combination)
// IIIIII
inline Numeric&         Tensor6View::operator()
  ( Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return                *(mdata +
			  OFFSET(v) + OFFSET(s) + OFFSET(b) +
			  OFFSET(p) + OFFSET(r) + OFFSET(c)    );
}



/** Return const iterator to sub-tensor. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator6D Tensor6View::begin() const
{
  return ConstTensor6View::begin();
}

/** Return const iterator behind last sub-tensor. */
inline ConstIterator6D Tensor6View::end() const
{
  return ConstTensor6View::end();
}

/** Return iterator to first sub-tensor. */
inline Iterator6D Tensor6View::begin()
{
  return Iterator6D( Tensor5View(mdata+mvr.mstart,
				msr,
				mbr,
				mpr,
				mrr,
				mcr),
		     mvr.mstride);
}

/** Return iterator behind last sub-tensor. */
inline Iterator6D Tensor6View::end()
{
  return Iterator6D( Tensor5View(mdata + mvr.mstart +
				(mvr.mextent)*mvr.mstride,
				msr,
				mbr,
				mpr,
				mrr,
				mcr),
		     mvr.mstride );
}

/** Assignment operator. This copies the data from another Tensor6View
    to this Tensor6View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor6View by
    setting its range. */
inline Tensor6View& Tensor6View::operator=(const ConstTensor6View& m)
{
  // Check that sizes are compatible:
  assert(mvr.mextent==m.mvr.mextent);
  assert(msr.mextent==m.msr.mextent);
  assert(mbr.mextent==m.mbr.mextent);
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from Tensor6View to Tensor6View. This is a tricky
    one. The problem is that since Tensor6View is derived from
    ConstTensor6View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline Tensor6View& Tensor6View::operator=(const Tensor6View& m)
{
  // Check that sizes are compatible:
  assert(mvr.mextent==m.mvr.mextent);
  assert(msr.mextent==m.msr.mextent);
  assert(mbr.mextent==m.mbr.mextent);
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Tensor6. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline Tensor6View& Tensor6View::operator=(const Tensor6& m)
{
  // Check that sizes are compatible:
  assert(mvr.mextent==m.mvr.mextent);
  assert(msr.mextent==m.msr.mextent);
  assert(mbr.mextent==m.mbr.mextent);
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor6View will set all elements to this
    value. */
inline Tensor6View& Tensor6View::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
inline Tensor6View& Tensor6View::operator*=(Numeric x)
{
  const Iterator6D ep=end();  
  for ( Iterator6D p=begin(); p!=ep ; ++p )
  {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
inline Tensor6View& Tensor6View::operator/=(Numeric x)
{
  const Iterator6D ep=end();  
  for ( Iterator6D p=begin(); p!=ep ; ++p )
  {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
inline Tensor6View& Tensor6View::operator+=(Numeric x)
{
  const Iterator6D ep=end();  
  for ( Iterator6D p=begin(); p!=ep ; ++p )
  {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
inline Tensor6View& Tensor6View::operator-=(Numeric x)
{
  const Iterator6D ep=end();  
  for ( Iterator6D p=begin(); p!=ep ; ++p )
  {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor6. */
inline Tensor6View& Tensor6View::operator*=(const ConstTensor6View& x)
{
  assert( nvitrines() == x.nvitrines() );
  assert( nshelves()  == x.nshelves()  );
  assert( nbooks()    == x.nbooks()    );
  assert( npages()    == x.npages()    );
  assert( nrows()     == x.nrows()     );
  assert( ncols()     == x.ncols()     );
  ConstIterator6D  xp = x.begin();
  Iterator6D        p = begin();
  const Iterator6D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p *= *xp;
    }
  return *this;
}

/** Element-vise division by another Tensor6. */
inline Tensor6View& Tensor6View::operator/=(const ConstTensor6View& x)
{
  assert( nvitrines() == x.nvitrines() );
  assert( nshelves()  == x.nshelves()  );
  assert( nbooks()    == x.nbooks()    );
  assert( npages()    == x.npages()    );
  assert( nrows()     == x.nrows()     );
  assert( ncols()     == x.ncols()     );
  ConstIterator6D  xp = x.begin();
  Iterator6D        p = begin();
  const Iterator6D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p /= *xp;
    }
  return *this;
}

/** Element-vise addition of another Tensor6. */
inline Tensor6View& Tensor6View::operator+=(const ConstTensor6View& x)
{
  assert( nvitrines() == x.nvitrines() );
  assert( nshelves()  == x.nshelves()  );
  assert( nbooks()    == x.nbooks()    );
  assert( npages()    == x.npages()    );
  assert( nrows()     == x.nrows()     );
  assert( ncols()     == x.ncols()     );
  ConstIterator6D  xp = x.begin();
  Iterator6D        p = begin();
  const Iterator6D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p += *xp;
    }
  return *this;
}

/** Element-vise subtraction of another Tensor6. */
inline Tensor6View& Tensor6View::operator-=(const ConstTensor6View& x)
{
  assert( nvitrines() == x.nvitrines() );
  assert( nshelves()  == x.nshelves()  );
  assert( nbooks()    == x.nbooks()    );
  assert( npages()    == x.npages()    );
  assert( nrows()     == x.nrows()     );
  assert( ncols()     == x.ncols()     );
  ConstIterator6D  xp = x.begin();
  Iterator6D        p = begin();
  const Iterator6D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p -= *xp;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Tensor6. */
inline Tensor6View::Tensor6View() :
  ConstTensor6View()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor6 to initialize its
    own Tensor6View part. The row range rr must have a
    stride to account for the length of one row. */
inline Tensor6View::Tensor6View(Numeric *data,
				const Range& v,
				const Range& s,
				const Range& b,
				const Range& p,
				const Range& r,
				const Range& c) :
  ConstTensor6View(data, v, s, b, p, r, c)
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
inline Tensor6View::Tensor6View(Numeric *data,
				const Range& pv,
				const Range& ps,
				const Range& pb,
				const Range& pp,
				const Range& pr,
				const Range& pc,
				const Range& nv,
				const Range& ns,
				const Range& nb,
				const Range& np,
				const Range& nr,
				const Range& nc) :
  ConstTensor6View(data,pv,ps,pb,pp,pr,pc,nv,ns,nb,np,nr,nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
inline void copy(ConstIterator6D origin,
                 const ConstIterator6D& end,
                 Iterator6D target)
{
  for ( ; origin!=end ; ++origin,++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy(origin->begin(),
	   origin->end(),
	   target->begin());
    }
}

/** Copy a scalar to all elements. */
inline void copy(Numeric x,
                 Iterator6D target,
                 const Iterator6D& end)
{
  for ( ; target!=end ; ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy(x,target->begin(),target->end());
    }
}


// Functions for Tensor6:
// ---------------------

/** Default constructor. */
inline Tensor6::Tensor6() :
  Tensor6View::Tensor6View()
{
  // Nothing to do here. However, note that the default constructor
  // for Tensor6View has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized. 
}

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
inline Tensor6::Tensor6(Index v, Index s, Index b,
                        Index p, Index r, Index c) :
  Tensor6View( new Numeric[v*s*b*p*r*c],
	       Range(0,v,s*b*p*r*c),
	       Range(0,s,b*p*r*c),
	       Range(0,b,p*r*c),
	       Range(0,p,r*c),
	       Range(0,r,c),
	       Range(0,c))
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Tensor6::Tensor6(Index v, Index s, Index b,
                        Index p, Index r, Index c, Numeric fill) :
  Tensor6View( new Numeric[v*s*b*p*r*c],
	       Range(0,v,s*b*p*r*c),
	       Range(0,s,b*p*r*c),
	       Range(0,b,p*r*c),
	       Range(0,p,r*c),
	       Range(0,r,c),
	       Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  const Numeric *stop = mdata+v*s*b*p*r*c;
  for ( Numeric *x=mdata; x<stop; ++x )
    *x = fill;
}

/** Copy constructor from Tensor6View. This automatically sets the size
    and copies the data. */
inline Tensor6::Tensor6(const ConstTensor6View& m) :
  Tensor6View( new Numeric[m.nvitrines()*m.nshelves()*m.nbooks(),
			   m.npages()*m.nrows()*m.ncols()],
	       Range( 0, m.nvitrines(), m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.nshelves(), m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.nbooks(), m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.npages(), m.nrows()*m.ncols() ),
	       Range( 0, m.nrows(), m.ncols() ),
	       Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Tensor6. This automatically sets the size
    and copies the data. */
inline Tensor6::Tensor6(const Tensor6& m) :
  Tensor6View( new Numeric[m.nvitrines()*m.nshelves()*m.nbooks(),
			   m.npages()*m.nrows()*m.ncols()],
	       Range( 0, m.nvitrines(), m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
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
  copy(m.begin(),m.end(),begin());
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
inline Tensor6& Tensor6::operator=(const Tensor6& m)
{
  //  cout << "Tensor6 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "             n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid tensor, so we just
  // have to check one.
  if ( 0 == mcr.mextent )
    {
      // Adjust if previously empty.
      resize( m.mvr.mextent, m.msr.mextent, m.mbr.mextent,
	      m.mpr.mextent, m.mrr.mextent, m.mcr.mextent ); 
    }
  else
    {
      // Check that sizes are compatible:
      assert( mvr.mextent==m.mvr.mextent );
      assert( msr.mextent==m.msr.mextent );
      assert( mbr.mextent==m.mbr.mextent );
      assert( mpr.mextent==m.mpr.mextent );
      assert( mrr.mextent==m.mrr.mextent );
      assert( mcr.mextent==m.mcr.mextent );
    }

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
inline Tensor6& Tensor6::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values. */
inline void Tensor6::resize(Index v, Index s, Index b,
			    Index p, Index r, Index c)
{
  assert( 0<=v );
  assert( 0<=s );
  assert( 0<=b );
  assert( 0<=p );
  assert( 0<=r );
  assert( 0<=c );

  if ( mvr.mextent!=v ||
       msr.mextent!=s ||
       mbr.mextent!=b ||
       mpr.mextent!=p ||
       mrr.mextent!=r ||
       mcr.mextent!=c )
    {
      delete mdata;
      mdata = new Numeric[v*s*b*p*r*c];

      mvr.mstart = 0;
      mvr.mextent = v;
      mvr.mstride = s*b*p*r*c;

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

/** Destructor for Tensor6. This is important, since Tensor6 uses new to
    allocate storage. */
inline Tensor6::~Tensor6()
{
//   cout << "Destroying a Tensor6:\n"
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

    \param   y Output:   The results of the function acting on each element of x.
    \param    my_func A function (e.g., sqrt).
    \param    x   A tensor. */
inline void transform( Tensor6View y,
                       double (&my_func)(double),
                       ConstTensor6View x )
{
  // Check dimensions:
  assert( y.nvitrines() == x.nvitrines() );
  assert( y.nshelves()  == x.nshelves()  );
  assert( y.nbooks()    == x.nbooks()    );
  assert( y.npages() 	== x.npages()    );
  assert( y.nrows()  	== x.nrows()     );
  assert( y.ncols()  	== x.ncols()     );

  const ConstIterator6D xe = x.end();
  ConstIterator6D       xi = x.begin();
  Iterator6D            yi = y.begin();
  for ( ; xi!=xe; ++xi, ++yi )
    {
      // Use the transform function of lower dimensional tensors
      // recursively:
      transform(*yi,my_func,*xi);
    }
}

/** Max function, tensor version. */
inline Numeric max(const ConstTensor6View& x)
{
  const ConstIterator6D xe = x.end();
  ConstIterator6D       xi = x.begin();

  // Initial value for max:
  Numeric themax = max(*xi);
  ++xi;

  for ( ; xi!=xe ; ++xi )
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
inline Numeric min(const ConstTensor6View& x)
{
  const ConstIterator6D xe = x.end();
  ConstIterator6D       xi = x.begin();

  // Initial value for min:
  Numeric themin = min(*xi);
  ++xi;

  for ( ; xi!=xe ; ++xi )
    {
      // Use the min function of lower dimensional tensors
      // recursively:
      Numeric mini = min(*xi);
      if ( mini < themin )
        themin = mini;
    }

  return themin;
}



#endif    // matpackVI_h
