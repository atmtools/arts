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
   Implementation of Tensors of Rank 7.

   Dimensions are called: library, vitrine, shelf, book, page, row, column.
   or short:              l,       v,       s,     b,    p,    r,   c
  
   \author Stefan Buehler
   \date   2001-11-22
*/

#ifndef matpackVII_h
#define matpackVII_h

#include <iomanip>
#include "matpackVI.h"

/** The outermost iterator class for rank 7 tensors. This takes into
    account the defined strided. */
class Iterator7D {
public:
  // Constructors:
  Iterator7D();
  Iterator7D(const Iterator7D& o);
  Iterator7D(const Tensor6View& x, Index stride);

  // Operators:
  Iterator7D& operator++();
  bool operator!=(const Iterator7D& other) const;
  Tensor6View* const operator->();
  Tensor6View& operator*();
  
private:
  /** Current position. */
  Tensor6View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator7D. */
class ConstIterator7D {
public:
  // Constructors:
  ConstIterator7D();
  ConstIterator7D(const ConstIterator7D& o);
  ConstIterator7D(const ConstTensor6View& x, Index stride);

  // Operators:
  ConstIterator7D& operator++();
  bool operator!=(const ConstIterator7D& other) const;
  const ConstTensor6View* operator->() const;
  const ConstTensor6View& operator*() const;

private:
  /** Current position. */
  ConstTensor6View msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor7:
class Tensor7;


/** A constant view of a Tensor7.

This, together with the derived class Tensor7View, contains the
main implementation of a Tensor7. It defines the concepts of
Tensor7View. Plus additionally the recursive subrange operator,
which makes it possible to create a Tensor7View from a subrange of
a Tensor7View.

Dimensions are called: library, vitrine, shelf, book, page, row, column.
or short:              l,       v,       s,     b,    p,    r,   c

The class Tensor7 is just a special case of a Tensor7View
which also allocates storage. */
class ConstTensor7View {
public:
  // Member functions:
  Index nlibraries() const;
  Index nvitrines()  const;
  Index nshelves()   const;
  Index nbooks()     const;
  Index npages()     const;
  Index nrows()      const;
  Index ncols()      const;

  // Const index operators:

  // Result 7D (1 combination)
  // -------
  ConstTensor7View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 6D (7 combinations)
  // ------|
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -----|-
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ----|--
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ---|---
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // --|----
  ConstTensor6View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-----
  ConstTensor6View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |------
  ConstTensor6View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ----|-|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ---|--|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // --|---|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -|----|
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-----|
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ----||-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ---|-|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // --|--|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // -|---|-
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |----|-
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ---||--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // --|-|--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|--|--
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |---|--
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // --||---
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-|---
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |--|---
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -||----
  ConstTensor5View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-|----
  ConstTensor5View operator()( Index l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-----
  ConstTensor5View operator()( Index l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---|-||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --|--||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -|---||
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |----||
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ---||-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|-|-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -|--|-|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |---|-|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // --||--|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|-|--|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |--|--|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -||---|
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-|---|
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||----|
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ---|||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --|-||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -|--||-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |---||-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // --||-|-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|-|-|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |--|-|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -||--|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |-|--|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||---|-
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // --|||--
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -|-||--
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |--||--
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||-|--
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|-|--
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||--|--
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|||---
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-||---
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-|---
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||----
  ConstTensor4View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||-|--
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||-||--
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|||--
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||||--
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |||--|-
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||-|-|-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |-||-|-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|||-|-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // ||--||-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |-|-||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -||-||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |--|||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|-|||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --||||-
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||---|
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||-|--|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |-||--|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|||--|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // ||--|-|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |-|-|-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -||-|-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |--||-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -|-||-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|||-|
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||---||
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |-|--||
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -||--||
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |--|-||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|-|-||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --||-||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |---|||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -|--|||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --|-|||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---||||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // ||||-|-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |||-||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ||-|||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |-||||-
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|||||-
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // ||||--|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |||-|-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ||-||-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |-|||-|
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -||||-|
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |||--||
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ||-|-||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |-||-||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|||-||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // ||--|||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |-|-|||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -||-|||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |--||||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -|-||||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // --|||||
  ConstMatrixView  operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 1D (7 combinations)
  // ||||||-
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||||-|
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||||-||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |||-|||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ||-||||
  ConstVectorView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // |-|||||
  ConstVectorView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -||||||
  ConstVectorView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result scalar (1 combination)
  // |||||||
  Numeric          operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;


  // Functions returning iterators:
  ConstIterator7D begin() const;
  ConstIterator7D end() const;
  
  // Friends:
  friend class Tensor7View;


protected:
  // Constructors:
  ConstTensor7View();
  ConstTensor7View(Numeric *data,
		   const Range& l,
		   const Range& v, const Range& s, const Range& b,
		   const Range& p, const Range& r, const Range& c);
  ConstTensor7View(Numeric *data,
		   const Range& pl,
		   const Range& pv, const Range& ps, const Range& pb,
		   const Range& pp, const Range& pr, const Range& pc,
		   const Range& nl,
		   const Range& nv, const Range& ns, const Range& nb,
		   const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The library range of mdata that is actually used. */
  Range mlr;
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

/** The Tensor7View class

This contains the main implementation of a Tensor7. It defines
the concepts of Tensor7View. Plus additionally the recursive
subrange operator, which makes it possible to create a Tensor7View
from a subrange of a Tensor7View. 

The class Tensor7 is just a special case of a Tensor7View
which also allocates storage. */
class Tensor7View : public ConstTensor7View {
public:

  // Const index operators:

  // Result 7D (1 combination)
  // -------
  ConstTensor7View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 6D (7 combinations)
  // ------|
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -----|-
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ----|--
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ---|---
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // --|----
  ConstTensor6View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-----
  ConstTensor6View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |------
  ConstTensor6View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ----|-|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ---|--|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // --|---|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -|----|
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-----|
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ----||-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ---|-|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // --|--|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // -|---|-
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |----|-
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ---||--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // --|-|--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|--|--
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |---|--
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // --||---
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-|---
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |--|---
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -||----
  ConstTensor5View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-|----
  ConstTensor5View operator()( Index l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-----
  ConstTensor5View operator()( Index l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---|-||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --|--||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -|---||
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |----||
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ---||-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|-|-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -|--|-|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |---|-|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // --||--|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|-|--|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |--|--|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -||---|
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-|---|
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||----|
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ---|||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --|-||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -|--||-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |---||-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // --||-|-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|-|-|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |--|-|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -||--|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |-|--|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||---|-
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // --|||--
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -|-||--
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |--||--
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||-|--
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|-|--
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||--|--
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|||---
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-||---
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-|---
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||----
  ConstTensor4View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||-|--
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||-||--
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|||--
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||||--
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |||--|-
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||-|-|-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |-||-|-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|||-|-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // ||--||-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |-|-||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -||-||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |--|||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|-|||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --||||-
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||---|
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||-|--|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |-||--|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|||--|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // ||--|-|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |-|-|-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -||-|-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |--||-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -|-||-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|||-|
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||---||
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |-|--||
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -||--||
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |--|-||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|-|-||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --||-||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |---|||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -|--|||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --|-|||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---||||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // ||||-|-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |||-||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ||-|||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |-||||-
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|||||-
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // ||||--|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |||-|-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ||-||-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |-|||-|
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -||||-|
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |||--||
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ||-|-||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |-||-||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|||-||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // ||--|||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |-|-|||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -||-|||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |--||||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -|-||||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // --|||||
  ConstMatrixView  operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 1D (7 combinations)
  // ||||||-
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||||-|
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||||-||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |||-|||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ||-||||
  ConstVectorView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // |-|||||
  ConstVectorView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -||||||
  ConstVectorView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result scalar (1 combination)
  // |||||||
  Numeric          operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;


  // Non-const index operators:

  // Result 7D (1 combination)
  // -------
  Tensor7View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 6D (7 combinations)
  // ------|
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // -----|-
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ----|--
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // ---|---
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // --|----
  Tensor6View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // -|-----
  Tensor6View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // |------
  Tensor6View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ----|-|
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // ---|--|
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // --|---|
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // -|----|
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // |-----|
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ----||-
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // ---|-|-
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // --|--|-
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // -|---|-
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // |----|-
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ---||--
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // --|-|--
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // -|--|--
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // |---|--
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // --||---
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // -|-|---
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |--|---
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // -||----
  Tensor5View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // |-|----
  Tensor5View operator()( Index l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // ||-----
  Tensor5View operator()( Index l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // ---|-||
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // --|--||
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // -|---||
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |----||
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ---||-|
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // --|-|-|
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // -|--|-|
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |---|-|
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // --||--|
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -|-|--|
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |--|--|
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -||---|
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // |-|---|
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ||----|
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ---|||-
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // --|-||-
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // -|--||-
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |---||-
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // --||-|-
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -|-|-|-
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |--|-|-
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -||--|-
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // |-|--|-
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ||---|-
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // --|||--
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -|-||--
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |--||--
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -||-|--
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // |-|-|--
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // ||--|--
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // -|||---
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |-||---
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // ||-|---
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |||----
  Tensor4View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |||-|--
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // ||-||--
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |-|||--
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -||||--
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |||--|-
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ||-|-|-
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |-||-|-
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -|||-|-
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // ||--||-
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |-|-||-
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // -||-||-
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |--|||-
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // -|-|||-
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // --||||-
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // |||---|
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ||-|--|
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |-||--|
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -|||--|
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // ||--|-|
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |-|-|-|
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // -||-|-|
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |--||-|
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // -|-||-|
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // --|||-|
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // ||---||
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |-|--||
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // -||--||
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |--|-||
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // -|-|-||
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // --||-||
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // |---|||
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // -|--|||
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // --|-|||
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // ---||||
  Tensor3View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // ||||-|-
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |||-||-
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // ||-|||-
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // |-||||-
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // -|||||-
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // ||||--|
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |||-|-|
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // ||-||-|
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // |-|||-|
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // -||||-|
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // |||--||
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ||-|-||
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // |-||-||
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // -|||-||
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // ||--|||
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // |-|-|||
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // -||-|||
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // |--||||
  MatrixView  operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // -|-||||
  MatrixView  operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // --|||||
  MatrixView  operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result 1D (7 combinations)
  // ||||||-
  VectorView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // |||||-|
  VectorView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // ||||-||
  VectorView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // |||-|||
  VectorView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // ||-||||
  VectorView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // |-|||||
  VectorView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);
  // -||||||
  VectorView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result scalar (1 combination)
  // |||||||
  Numeric&         operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c);


  // Functions returning const iterators:
  ConstIterator7D begin() const;
  ConstIterator7D end() const;
  // Functions returning iterators:
  Iterator7D begin();
  Iterator7D end();
  
  // Assignment operators:
  Tensor7View& operator=(const ConstTensor7View& v);
  Tensor7View& operator=(const Tensor7View& v);
  Tensor7View& operator=(const Tensor7& v);
  Tensor7View& operator=(Numeric x);

  // Other operators:
  Tensor7View& operator*=(Numeric x);
  Tensor7View& operator/=(Numeric x);
  Tensor7View& operator+=(Numeric x);
  Tensor7View& operator-=(Numeric x);

  Tensor7View& operator*=(const ConstTensor7View& x);
  Tensor7View& operator/=(const ConstTensor7View& x);
  Tensor7View& operator+=(const ConstTensor7View& x);
  Tensor7View& operator-=(const ConstTensor7View& x);

  // Friends:

protected:
  // Constructors:
  Tensor7View();
  Tensor7View(Numeric *data,
	      const Range& l,
	      const Range& v, const Range& s, const Range& b,
	      const Range& p, const Range& r, const Range& c);
  Tensor7View(Numeric *data,
	      const Range& pl,
	      const Range& pv, const Range& ps, const Range& pb,
	      const Range& pp, const Range& pr, const Range& pc,
	      const Range& nl,
	      const Range& nv, const Range& ns, const Range& nb,
	      const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor7 class. This is a Tensor7View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor7View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor7 : public Tensor7View {
public:
  // Constructors:
  Tensor7();
  Tensor7(Index        l,
	  Index        v, Index        s, Index        b,
	  Index        p, Index        r, Index        c);
  Tensor7(Index        l,
	  Index        v, Index        s, Index        b,
	  Index        p, Index        r, Index        c,
	  Numeric fill);
  Tensor7(const ConstTensor7View& v);
  Tensor7(const Tensor7& v);

  // Assignment operators:
  Tensor7& operator=(const Tensor7& x);
  Tensor7& operator=(Numeric x);

  // Resize function:
  void resize(Index        l,
	      Index        v, Index        s, Index        b,
	      Index        p, Index        r, Index        c);

  // Destructor:
  ~Tensor7();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator7D origin,
                 const ConstIterator7D& end,
                 Iterator7D target);

inline void copy(Numeric x,
                 Iterator7D target,
                 const Iterator7D& end);



// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor7. */
typedef Array<Tensor7> ArrayOfTensor7;



// Functions for Iterator7D
// ------------------------

/** Default constructor. */
inline Iterator7D::Iterator7D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator7D::Iterator7D(const Iterator7D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator7D::Iterator7D(const Tensor6View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline Iterator7D& Iterator7D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy.
    FIXME: Is it really necessary to have such a complicated check
    here? It could be sufficient to just test
    msv.mdata!=other.msv.mdata. */
inline bool Iterator7D::operator!=(const Iterator7D& other) const
{
  if ( msv.mdata +
       msv.mvr.mstart +
       msv.msr.mstart +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart 
       !=
       other.msv.mdata +
       other.msv.mvr.mstart +
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
inline Tensor6View* const Iterator7D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline Tensor6View& Iterator7D::operator*()
{
  return msv;
}

// Functions for ConstIterator7D
// -----------------------------

/** Default constructor. */
inline ConstIterator7D::ConstIterator7D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator7D::ConstIterator7D(const ConstIterator7D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator7D::ConstIterator7D(const ConstTensor6View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here. 
}

/** Prefix increment operator. */
inline ConstIterator7D& ConstIterator7D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. 
    FIXME: Is it really necessary to have such a complicated check
    here? It could be sufficient to just test
    msv.mdata!=other.msv.mdata. */
inline bool ConstIterator7D::operator!=(const ConstIterator7D& other) const
{
  if ( msv.mdata +
       msv.mvr.mstart +
       msv.msr.mstart +
       msv.mbr.mstart +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mvr.mstart +
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
inline const ConstTensor6View* ConstIterator7D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstTensor6View& ConstIterator7D::operator*() const
{
  return msv;
}



// Functions for ConstTensor7View:
// ------------------------------

/** Returns the number of libraries. */
inline Index ConstTensor7View::nlibraries() const
{
  return mlr.mextent;
}

/** Returns the number of vitrines. */
inline Index ConstTensor7View::nvitrines() const
{
  return mvr.mextent;
}

/** Returns the number of shelves. */
inline Index ConstTensor7View::nshelves() const
{
  return msr.mextent;
}

/** Returns the number of books. */
inline Index ConstTensor7View::nbooks() const
{
  return mbr.mextent;
}

/** Returns the number of pages. */
inline Index ConstTensor7View::npages() const
{
  return mpr.mextent;
}

/** Returns the number of rows. */
inline Index ConstTensor7View::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor7View::ncols() const
{
  return mcr.mextent;
}

// Const index operators:

// -------
inline ConstTensor7View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View(mdata,
			  mlr, mvr, msr, mbr, mpr, mrr, mcr,
			  l,   v,   s,   b,   p,   r,   c);
}
// |------
inline ConstTensor6View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  return ConstTensor6View(mdata + OFFSET(l),
			  mvr, msr, mbr, mpr, mrr, mcr,
			  v,   s,   b,   p,   r,   c);
}

// ------|
inline ConstTensor6View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(c);
  return ConstTensor6View(mdata + OFFSET(c),
                          mlr, mvr, msr, mbr, mpr, mrr,
                          l,   v,   s,   b,   p,   r);
}
// |-----|
inline ConstTensor5View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(l) + OFFSET(c),
                          mvr, msr, mbr, mpr, mrr,
                          v,   s,   b,   p,   r);
}

// -----|-
inline ConstTensor6View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(r);
  return ConstTensor6View(mdata + OFFSET(r),
                          mlr, mvr, msr, mbr, mpr, mcr,
                          l,   v,   s,   b,   p,   c);
}
// |----|-
inline ConstTensor5View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(r);
  return ConstTensor5View(mdata + OFFSET(l) + OFFSET(r),
                          mvr, msr, mbr, mpr, mcr,
                          v,   s,   b,   p,   c);
}

// ----|--
inline ConstTensor6View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(p);
  return ConstTensor6View(mdata + OFFSET(p),
                          mlr, mvr, msr, mbr, mrr, mcr,
                          l,   v,   s,   b,   r,   c);
}
// |---|--
inline ConstTensor5View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(p);
  return ConstTensor5View(mdata + OFFSET(l) + OFFSET(p),
                          mvr, msr, mbr, mrr, mcr,
                          v,   s,   b,   r,   c);
}

// ---|---
inline ConstTensor6View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(b);
  return ConstTensor6View(mdata + OFFSET(b),
                          mlr, mvr, msr, mpr, mrr, mcr,
                          l,   v,   s,   p,   r,   c);
}
// |--|---
inline ConstTensor5View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(b);
  return ConstTensor5View(mdata + OFFSET(l) + OFFSET(b),
                          mvr, msr, mpr, mrr, mcr,
                          v,   s,   p,   r,   c);
}

// --|----
inline ConstTensor6View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(s);
  return ConstTensor6View(mdata + OFFSET(s),
                          mlr, mvr, mbr, mpr, mrr, mcr,
                          l,   v,   b,   p,   r,   c);
}
// |-|----
inline ConstTensor5View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  return ConstTensor5View(mdata + OFFSET(l) + OFFSET(s),
                          mvr, mbr, mpr, mrr, mcr,
                          v,   b,   p,   r,   c);
}

// -|-----
inline ConstTensor6View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  return ConstTensor6View(mdata + OFFSET(v),
                          mlr, msr, mbr, mpr, mrr, mcr,
                          l,   s,   b,   p,   r,   c);
}
// ||-----
inline ConstTensor5View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  return ConstTensor5View(mdata + OFFSET(l) + OFFSET(v),
                          msr, mbr, mpr, mrr, mcr,
                          s,   b,   p,   r,   c);
}

// -----||
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(r);
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(r) + OFFSET(c),
                          mlr, mvr, msr, mbr, mpr,
                          l,   v,   s,   b,   p   );
}
// |----||
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mbr, mpr,
                          v,   s,   b,   p   );
}

// ----|-|
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(p);
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(p) + OFFSET(c),
                          mlr, mvr, msr, mbr, mrr,
                          l,   v,   s,   b,   r    );
}
// |---|-|
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(c),
                          mvr, msr, mbr, mrr,
                          v,   s,   b,   r    );
}

// ---|--|
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(b);
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(b) + OFFSET(c),
                          mlr, mvr, msr, mpr, mrr,
                          l,   v,   s,   p,   r    );
}
// |--|--|
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(c),
                          mvr, msr, mpr, mrr,
                          v,   s,   p,   r    );
}

// --|---|
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(s) + OFFSET(c),
                          mlr, mvr, mbr, mpr, mrr,
                          l,   v,   b,   p,   r    );
}
// |-|---|
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(c),
                          mvr, mbr, mpr, mrr,
                          v,   b,   p,   r    );
}

// -|----|
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(c);
  return ConstTensor5View(mdata + OFFSET(v) + OFFSET(c),
                          mlr, msr, mbr, mpr, mrr,
                          l,   s,   b,   p,   r    );
}
// ||----|
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(c),
                          msr, mbr, mpr, mrr,
                          s,   b,   p,   r    );
}

// ----||-
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(p);
  CHECK(r);
  return ConstTensor5View(mdata + OFFSET(p) + OFFSET(r),
                          mlr, mvr, msr, mbr, mcr,
                          l,   v,   s,   b,   c    );
}
// |---||-
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r),
                          mvr, msr, mbr, mcr,
                          v,   s,   b,   c    );
}

// ---|-|-
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(b);
  CHECK(r);
  return ConstTensor5View(mdata + OFFSET(b) + OFFSET(r),
                          mlr, mvr, msr, mpr, mcr,
                          l,   v,   s,   p,   c    );
}
// |--|-|-
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r),
                          mvr, msr, mpr, mcr,
                          v,   s,   p,   c    );
}

// --|--|-
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(r);
  return ConstTensor5View(mdata + OFFSET(s) + OFFSET(r),
                          mlr, mvr, mbr, mpr, mcr,
                          l,   v,   b,   p,   c    );
}
// |-|--|-
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r),
                          mvr, mbr, mpr, mcr,
                          v,   b,   p,   c    );
}

// -|---|-
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(r);
  return ConstTensor5View(mdata + OFFSET(v) + OFFSET(r),
                          mlr, msr, mbr, mpr, mcr,
                          l,   s,   b,   p,   c    );
}
// ||---|-
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r),
                          msr, mbr, mpr, mcr,
                          s,   b,   p,   c    );
}

// ---||--
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(b);
  CHECK(p);
  return ConstTensor5View(mdata + OFFSET(b) + OFFSET(p),
                          mlr, mvr, msr, mrr, mcr,
                          l,   v,   s,   r,   c    );
}
// |--||--
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p),
                          mvr, msr, mrr, mcr,
                          v,   s,   r,   c    );
}

// --|-|--
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(s);
  CHECK(p);
  return ConstTensor5View(mdata + OFFSET(s) + OFFSET(p),
                          mlr, mvr, mbr, mrr, mcr,
                          l,   v,   b,   r,   c    );
}
// |-|-|--
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p),
                          mvr, mbr, mrr, mcr,
                          v,   b,   r,   c    );
}

// -|--|--
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(p);
  return ConstTensor5View(mdata + OFFSET(v) + OFFSET(p),
                          mlr, msr, mbr, mrr, mcr,
                          l,   s,   b,   r,   c    );
}
// ||--|--
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p),
                          msr, mbr, mrr, mcr,
                          s,   b,   r,   c    );
}

// --||---
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  return ConstTensor5View(mdata + OFFSET(s) + OFFSET(b),
                          mlr, mvr, mpr, mrr, mcr,
                          l,   v,   p,   r,   c    );
}
// |-||---
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b),
                          mvr, mpr, mrr, mcr,
                          v,   p,   r,   c    );
}

// -|-|---
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  return ConstTensor5View(mdata + OFFSET(v) + OFFSET(b),
                          mlr, msr, mpr, mrr, mcr,
                          l,   s,   p,   r,   c    );
}
// ||-|---
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b),
                          msr, mpr, mrr, mcr,
                          s,   p,   r,   c    );
}

// -||----
inline ConstTensor5View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  return ConstTensor5View(mdata + OFFSET(v) + OFFSET(s),
                          mlr, mbr, mpr, mrr, mcr,
                          l,   b,   p,   r,   c    );
}
// |||----
inline ConstTensor4View ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s),
                          mbr, mpr, mrr, mcr,
                          b,   p,   r,   c    );
}

// ----|||
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr, mvr, msr, mbr, 
                          l,   v,   s,   b     );
}
// |---|||
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mbr, 
                          v,   s,   b     );
}

// ---|-||
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mlr, mvr, msr, mpr, 
                          l,   v,   s,   p     );
}
// |--|-||
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mvr, msr, mpr, 
                          v,   s,   p     );
}

// --|--||
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mlr, mvr, mbr, mpr, 
                          l,   v,   b,   p     );
}
// |-|--||
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mvr, mbr, mpr, 
                          v,   b,   p     );
}

// -|---||
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(r) + OFFSET(c),
                          mlr, msr, mbr, mpr, 
                          l,   s,   b,   p     );
}
// ||---||
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r) + OFFSET(c),
                          msr, mbr, mpr, 
                          s,   b,   p     );
}

// ---||-|
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mlr, mvr, msr, mrr, 
                          l,   v,   s,   r     );
}
// |--||-|
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mvr, msr, mrr, 
                          v,   s,   r     );
}

// --|-|-|
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mlr, mvr, mbr, mrr, 
                          l,   v,   b,   r     );
}
// |-|-|-|
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mvr, mbr, mrr, 
                          v,   b,   r     );
}

// -|--|-|
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(c),
                          mlr, msr, mbr, mrr, 
                          l,   s,   b,   r     );
}
// ||--|-|
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(c),
                          msr, mbr, mrr, 
                          s,   b,   r     );
}

// --||--|
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mlr, mvr, mpr, mrr, 
                          l,   v,   p,   r     );
}
// |-||--|
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mvr, mpr, mrr, 
                          v,   p,   r     );
}

// -|-|--|
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(c),
                          mlr, msr, mpr, mrr, 
                          l,   s,   p,   r     );
}
// ||-|--|
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(c),
                          msr, mpr, mrr, 
                          s,   p,   r     );
}

// -||---|
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(c),
                          mlr, mbr, mpr, mrr, 
                          l,   b,   p,   r     );
}
// |||---|
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(c),
                          mbr, mpr, mrr, 
                          b,   p,   r     );
}

// ---|||-
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mlr, mvr, msr, mcr, 
                          l,   v,   s,   c     );
}
// |--|||-
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mvr, msr, mcr, 
                          v,   s,   c     );
}

// --|-||-
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mlr, mvr, mbr, mcr, 
                          l,   v,   b,   c     );
}
// |-|-||-
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mvr, mbr, mcr, 
                          v,   b,   c     );
}

// -|--||-
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r),
                          mlr, msr, mbr, mcr, 
                          l,   s,   b,   c     );
}
// ||--||-
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r),
                          msr, mbr, mcr, 
                          s,   b,   c     );
}

// --||-|-
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mlr, mvr, mpr, mcr, 
                          l,   v,   p,   c     );
}
// |-||-|-
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mvr, mpr, mcr, 
                          v,   p,   c     );
}

// -|-|-|-
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r),
                          mlr, msr, mpr, mcr, 
                          l,   s,   p,   c     );
}
// ||-|-|-
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r),
                          msr, mpr, mcr, 
                          s,   p,   c     );
}

// -||--|-
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r),
                          mlr, mbr, mpr, mcr, 
                          l,   b,   p,   c     );
}
// |||--|-
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r),
                          mbr, mpr, mcr, 
                          b,   p,   c     );
}

// --|||--
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mlr, mvr, mrr, mcr, 
                          l,   v,   r,   c     );
}
// |-|||--
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mvr, mrr, mcr, 
                          v,   r,   c     );
}

// -|-||--
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p),
                          mlr, msr, mrr, mcr, 
                          l,   s,   r,   c     );
}
// ||-||--
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p),
                          msr, mrr, mcr, 
                          s,   r,   c     );
}

// -||-|--
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p),
                          mlr, mbr, mrr, mcr, 
                          l,   b,   r,   c     );
}
// |||-|--
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p),
                          mbr, mrr, mcr, 
                          b,   r,   c     );
}

// -|||---
inline ConstTensor4View ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b),
                          mlr, mpr, mrr, mcr, 
                          l,   p,   r,   c     );
}
// ||||---
inline ConstTensor3View ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b),
                          mpr, mrr, mcr, 
                          p,   r,   c     );
}

// -||||--
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
			   mlr, mrr, mcr, 
			   l,   r,   c     );
}
// |||||--
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
			  mrr, mcr, 
			  r,   c     );
}

// -|||-|-
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
			   mlr, mpr, mcr, 
			   l,   p,   c     );
}
// ||||-|-
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
			  mpr, mcr, 
			  p,   c     );
}

// -||-||-
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
			   mlr, mbr, mcr, 
			   l,   b,   c     );
}
// |||-||-
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
			  mbr, mcr, 
			  b,   c     );
}

// -|-|||-
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
			   mlr, msr, mcr, 
			   l,   s,   c     );
}
// ||-|||-
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
			  msr, mcr, 
			  s,   c     );
}

// --||||-
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
			   mlr, mvr, mcr, 
			   l,   v,   c     );
}
// |-||||-
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
			  mvr, mcr, 
			  v,   c     );
}

// -|||--|
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
			   mlr, mpr, mrr, 
			   l,   p,   r     );
}
// ||||--|
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
			  mpr, mrr, 
			  p,   r     );
}

// -||-|-|
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
			   mlr, mbr, mrr, 
			   l,   b,   r     );
}
// |||-|-|
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
			  mbr, mrr, 
			  b,   r     );
}

// -|-||-|
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
			   mlr, msr, mrr, 
			   l,   s,   r     );
}
// ||-||-|
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
			  msr, mrr, 
			  s,   r     );
}

// --|||-|
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
			   mlr, mvr, mrr, 
			   l,   v,   r     );
}
// |-|||-|
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
			  mvr, mrr, 
			  v,   r     );
}

// -||--||
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
			   mlr, mbr, mpr, 
			   l,   b,   p     );
}
// |||--||
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
			  mbr, mpr, 
			  b,   p     );
}

// -|-|-||
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
			   mlr, msr, mpr, 
			   l,   s,   p     );
}
// ||-|-||
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
			  msr, mpr, 
			  s,   p     );
}

// --||-||
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
			   mlr, mvr, mpr, 
			   l,   v,   p     );
}
// |-||-||
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
			  mvr, mpr, 
			  v,   p     );
}

// -|--|||
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
			   mlr, msr, mbr, 
			   l,   s,   b     );
}
// ||--|||
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
			  msr, mbr, 
			  s,   b     );
}

// --|-|||
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
			   mlr, mvr, mbr, 
			   l,   v,   b     );
}
// |-|-|||
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
			  mvr, mbr, 
			  v,   b     );
}

// ---||||
inline ConstTensor3View  ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstTensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
			   mlr, mvr, msr, 
			   l,   v,   s     );
}
// |--||||
inline ConstMatrixView  ConstTensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
			  mvr, msr, 
			  v,   s     );
}

// -|||||-
inline ConstMatrixView  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mlr, mcr,  
                          l,   c      );
}
// ||||||-
inline ConstVectorView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mcr,  
                          c      );
}

// -||||-|
inline ConstMatrixView  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mlr, mrr,  
                          l,   r      );
}
// |||||-|
inline ConstVectorView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mrr,  
                          r      );
}

// -|||-||
inline ConstMatrixView  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mlr, mpr,  
                          l,   p      );
}
// ||||-||
inline ConstVectorView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mpr,  
                          p      );
}

// -||-|||
inline ConstMatrixView  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr, mbr,  
                          l,   b      );
}
// |||-|||
inline ConstVectorView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mbr,  
                          b      );
}

// -|-||||
inline ConstMatrixView  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr, msr,  
                          l,   s      );
}
// ||-||||
inline ConstVectorView  ConstTensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          msr,  
                          s      );
}

// --|||||
inline ConstMatrixView  ConstTensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstMatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr, mvr,  
                          l,   v      );
}
// |-|||||
inline ConstVectorView  ConstTensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr,  
                          v      );
}

// -||||||
inline ConstVectorView  ConstTensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  ConstVectorView( mdata +
			   OFFSET(v) + OFFSET(s) + OFFSET(b) +
			   OFFSET(p) + OFFSET(r) + OFFSET(c),
			   mlr,
			   l    );
}
// |||||||
inline Numeric          ConstTensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return                *(mdata + OFFSET(l) +
			  OFFSET(v) + OFFSET(s) + OFFSET(b) +
			  OFFSET(p) + OFFSET(r) + OFFSET(c)    );
}


/** Return const iterator to first sub-tensor. */
inline ConstIterator7D ConstTensor7View::begin() const
{
  return ConstIterator7D( ConstTensor6View(mdata+mlr.mstart,
					   mvr,
					   msr,
					   mbr,
					   mpr,
					   mrr,
					   mcr),
			  mlr.mstride);
}

/** Return const iterator behind last sub-tensor. */
inline ConstIterator7D ConstTensor7View::end() const
{
  return ConstIterator7D( ConstTensor6View(mdata + mlr.mstart +
                                          (mlr.mextent)*mlr.mstride,
					   mvr,
					   msr,
					   mbr,
					   mpr,
					   mrr,
					   mcr),
                          mlr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor7View::ConstTensor7View() :
  mlr(0,0,1),
  mvr(0,0,1), msr(0,0,1), mbr(0,0,1), 
  mpr(0,0,1), mrr(0,0,1), mcr(0,0,1),
  mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor7 to initialize
    its own Tensor7View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
inline ConstTensor7View::ConstTensor7View(Numeric *data,
					  const Range& l,
                                          const Range& v,
                                          const Range& s,
                                          const Range& b,
                                          const Range& p,
                                          const Range& r,
                                          const Range& c) :
  mlr(l),
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
inline ConstTensor7View::ConstTensor7View(Numeric *data,
					  const Range& pl,
                                          const Range& pv,
                                          const Range& ps,
                                          const Range& pb,
                                          const Range& pp,
                                          const Range& pr,
                                          const Range& pc,
					  const Range& nl,
                                          const Range& nv,
                                          const Range& ns,
                                          const Range& nb,
                                          const Range& np,
                                          const Range& nr,
                                          const Range& nc) :
  mlr(pl,nl),
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
    Tensor6 to print each page in turn. */
inline std::ostream& operator<<(std::ostream& os, const ConstTensor7View& v)
{
  // Page iterators:
  ConstIterator7D ip=v.begin();
  const ConstIterator7D end_page=v.end();

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


// Functions for Tensor7View:
// -------------------------

// Const index operators:

// -------
inline ConstTensor7View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |------
inline ConstTensor6View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ------|
inline ConstTensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-----|
inline ConstTensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -----|-
inline ConstTensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |----|-
inline ConstTensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ----|--
inline ConstTensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |---|--
inline ConstTensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---|---
inline ConstTensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--|---
inline ConstTensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|----
inline ConstTensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|----
inline ConstTensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-----
inline ConstTensor6View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-----
inline ConstTensor5View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -----||
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |----||
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ----|-|
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |---|-|
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---|--|
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--|--|
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|---|
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|---|
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|----|
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||----|
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ----||-
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |---||-
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---|-|-
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--|-|-
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|--|-
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|--|-
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|---|-
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||---|-
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---||--
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--||--
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|-|--
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|-|--
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|--|--
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||--|--
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --||---
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-||---
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-|---
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-|---
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||----
inline ConstTensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||----
inline ConstTensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ----|||
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |---|||
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---|-||
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--|-||
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|--||
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|--||
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|---||
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||---||
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---||-|
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--||-|
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|-|-|
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|-|-|
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|--|-|
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||--|-|
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --||--|
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-||--|
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-|--|
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-|--|
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||---|
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||---|
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---|||-
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--|||-
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|-||-
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|-||-
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|--||-
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||--||-
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --||-|-
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-||-|-
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-|-|-
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-|-|-
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||--|-
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||--|-
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|||--
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|||--
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-||--
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-||--
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||-|--
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||-|--
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|||---
inline ConstTensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||||---
inline ConstTensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||||--
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||||--
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|||-|-
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||||-|-
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||-||-
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||-||-
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-|||-
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-|||-
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --||||-
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-||||-
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|||--|
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||||--|
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||-|-|
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||-|-|
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-||-|
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-||-|
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|||-|
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|||-|
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||--||
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||--||
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-|-||
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-|-||
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --||-||
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-||-||
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|--|||
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||--|||
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|-|||
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|-|||
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// ---||||
inline ConstTensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |--||||
inline ConstMatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|||||-
inline ConstMatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||||||-
inline ConstVectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||||-|
inline ConstMatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||||-|
inline ConstVectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|||-||
inline ConstMatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||||-||
inline ConstVectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||-|||
inline ConstMatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||-|||
inline ConstVectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -|-||||
inline ConstMatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// ||-||||
inline ConstVectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// --|||||
inline ConstMatrixView  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |-|||||
inline ConstVectorView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}

// -||||||
inline ConstVectorView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}
// |||||||
inline Numeric          Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) const
{
  return ConstTensor7View::operator()(l,v,s,b,p,r,c);    
}


// Non-const index operators:

// -------
inline Tensor7View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  return Tensor7View(mdata,
		     mlr, mvr, msr, mbr, mpr, mrr, mcr,
		     l,   v,   s,   b,   p,   r,   c);
}
// |------
inline Tensor6View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  return Tensor6View(mdata + OFFSET(l),
		     mvr, msr, mbr, mpr, mrr, mcr,
		     v,   s,   b,   p,   r,   c);
}

// ------|
inline Tensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(c);
  return Tensor6View(mdata + OFFSET(c),
		     mlr, mvr, msr, mbr, mpr, mrr,
		     l,   v,   s,   b,   p,   r);
}
// |-----|
inline Tensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(c);
  return Tensor5View(mdata + OFFSET(l) + OFFSET(c),
		     mvr, msr, mbr, mpr, mrr,
		     v,   s,   b,   p,   r);
}

// -----|-
inline Tensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(r);
  return Tensor6View(mdata + OFFSET(r),
		     mlr, mvr, msr, mbr, mpr, mcr,
		     l,   v,   s,   b,   p,   c);
}
// |----|-
inline Tensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(r);
  return Tensor5View(mdata + OFFSET(l) + OFFSET(r),
		     mvr, msr, mbr, mpr, mcr,
		     v,   s,   b,   p,   c);
}

// ----|--
inline Tensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(p);
  return Tensor6View(mdata + OFFSET(p),
		     mlr, mvr, msr, mbr, mrr, mcr,
		     l,   v,   s,   b,   r,   c);
}
// |---|--
inline Tensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(p);
  return Tensor5View(mdata + OFFSET(l) + OFFSET(p),
		     mvr, msr, mbr, mrr, mcr,
		     v,   s,   b,   r,   c);
}

// ---|---
inline Tensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(b);
  return Tensor6View(mdata + OFFSET(b),
		     mlr, mvr, msr, mpr, mrr, mcr,
		     l,   v,   s,   p,   r,   c);
}
// |--|---
inline Tensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(b);
  return Tensor5View(mdata + OFFSET(l) + OFFSET(b),
		     mvr, msr, mpr, mrr, mcr,
		     v,   s,   p,   r,   c);
}

// --|----
inline Tensor6View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(s);
  return Tensor6View(mdata + OFFSET(s),
		     mlr, mvr, mbr, mpr, mrr, mcr,
		     l,   v,   b,   p,   r,   c);
}
// |-|----
inline Tensor5View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  return Tensor5View(mdata + OFFSET(l) + OFFSET(s),
		     mvr, mbr, mpr, mrr, mcr,
		     v,   b,   p,   r,   c);
}

// -|-----
inline Tensor6View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  return Tensor6View(mdata + OFFSET(v),
		     mlr, msr, mbr, mpr, mrr, mcr,
		     l,   s,   b,   p,   r,   c);
}
// ||-----
inline Tensor5View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  return Tensor5View(mdata + OFFSET(l) + OFFSET(v),
		     msr, mbr, mpr, mrr, mcr,
		     s,   b,   p,   r,   c);
}

// -----||
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(r);
  CHECK(c);
  return Tensor5View(mdata + OFFSET(r) + OFFSET(c),
		     mlr, mvr, msr, mbr, mpr,
		     l,   v,   s,   b,   p   );
}
// |----||
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(r) + OFFSET(c),
		     mvr, msr, mbr, mpr,
		     v,   s,   b,   p   );
}

// ----|-|
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(p);
  CHECK(c);
  return Tensor5View(mdata + OFFSET(p) + OFFSET(c),
		     mlr, mvr, msr, mbr, mrr,
		     l,   v,   s,   b,   r    );
}
// |---|-|
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(c),
		     mvr, msr, mbr, mrr,
		     v,   s,   b,   r    );
}

// ---|--|
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(b);
  CHECK(c);
  return Tensor5View(mdata + OFFSET(b) + OFFSET(c),
		     mlr, mvr, msr, mpr, mrr,
		     l,   v,   s,   p,   r    );
}
// |--|--|
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(c),
		     mvr, msr, mpr, mrr,
		     v,   s,   p,   r    );
}

// --|---|
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(c);
  return Tensor5View(mdata + OFFSET(s) + OFFSET(c),
		     mlr, mvr, mbr, mpr, mrr,
		     l,   v,   b,   p,   r    );
}
// |-|---|
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(c),
		     mvr, mbr, mpr, mrr,
		     v,   b,   p,   r    );
}

// -|----|
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(c);
  return Tensor5View(mdata + OFFSET(v) + OFFSET(c),
		     mlr, msr, mbr, mpr, mrr,
		     l,   s,   b,   p,   r    );
}
// ||----|
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(c),
		     msr, mbr, mpr, mrr,
		     s,   b,   p,   r    );
}

// ----||-
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(p);
  CHECK(r);
  return Tensor5View(mdata + OFFSET(p) + OFFSET(r),
		     mlr, mvr, msr, mbr, mcr,
		     l,   v,   s,   b,   c    );
}
// |---||-
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r),
		     mvr, msr, mbr, mcr,
		     v,   s,   b,   c    );
}

// ---|-|-
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(b);
  CHECK(r);
  return Tensor5View(mdata + OFFSET(b) + OFFSET(r),
		     mlr, mvr, msr, mpr, mcr,
		     l,   v,   s,   p,   c    );
}
// |--|-|-
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r),
		     mvr, msr, mpr, mcr,
		     v,   s,   p,   c    );
}

// --|--|-
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(r);
  return Tensor5View(mdata + OFFSET(s) + OFFSET(r),
		     mlr, mvr, mbr, mpr, mcr,
		     l,   v,   b,   p,   c    );
}
// |-|--|-
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r),
		     mvr, mbr, mpr, mcr,
		     v,   b,   p,   c    );
}

// -|---|-
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(r);
  return Tensor5View(mdata + OFFSET(v) + OFFSET(r),
		     mlr, msr, mbr, mpr, mcr,
		     l,   s,   b,   p,   c    );
}
// ||---|-
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r),
		     msr, mbr, mpr, mcr,
		     s,   b,   p,   c    );
}

// ---||--
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(b);
  CHECK(p);
  return Tensor5View(mdata + OFFSET(b) + OFFSET(p),
		     mlr, mvr, msr, mrr, mcr,
		     l,   v,   s,   r,   c    );
}
// |--||--
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p),
		     mvr, msr, mrr, mcr,
		     v,   s,   r,   c    );
}

// --|-|--
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(s);
  CHECK(p);
  return Tensor5View(mdata + OFFSET(s) + OFFSET(p),
		     mlr, mvr, mbr, mrr, mcr,
		     l,   v,   b,   r,   c    );
}
// |-|-|--
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p),
		     mvr, mbr, mrr, mcr,
		     v,   b,   r,   c    );
}

// -|--|--
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(p);
  return Tensor5View(mdata + OFFSET(v) + OFFSET(p),
		     mlr, msr, mbr, mrr, mcr,
		     l,   s,   b,   r,   c    );
}
// ||--|--
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p),
		     msr, mbr, mrr, mcr,
		     s,   b,   r,   c    );
}

// --||---
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  return Tensor5View(mdata + OFFSET(s) + OFFSET(b),
		     mlr, mvr, mpr, mrr, mcr,
		     l,   v,   p,   r,   c    );
}
// |-||---
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b),
		     mvr, mpr, mrr, mcr,
		     v,   p,   r,   c    );
}

// -|-|---
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  return Tensor5View(mdata + OFFSET(v) + OFFSET(b),
		     mlr, msr, mpr, mrr, mcr,
		     l,   s,   p,   r,   c    );
}
// ||-|---
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b),
		     msr, mpr, mrr, mcr,
		     s,   p,   r,   c    );
}

// -||----
inline Tensor5View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  return Tensor5View(mdata + OFFSET(v) + OFFSET(s),
		     mlr, mbr, mpr, mrr, mcr,
		     l,   b,   p,   r,   c    );
}
// |||----
inline Tensor4View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s),
		     mbr, mpr, mrr, mcr,
		     b,   p,   r,   c    );
}

// ----|||
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mlr, mvr, msr, mbr, 
		     l,   v,   s,   b     );
}
// |---|||
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mvr, msr, mbr, 
		     v,   s,   b     );
}

// ---|-||
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(r) + OFFSET(c),
		     mlr, mvr, msr, mpr, 
		     l,   v,   s,   p     );
}
// |--|-||
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		     mvr, msr, mpr, 
		     v,   s,   p     );
}

// --|--||
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(r) + OFFSET(c),
		     mlr, mvr, mbr, mpr, 
		     l,   v,   b,   p     );
}
// |-|--||
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r) + OFFSET(c),
		     mvr, mbr, mpr, 
		     v,   b,   p     );
}

// -|---||
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(r) + OFFSET(c),
		     mlr, msr, mbr, mpr, 
		     l,   s,   b,   p     );
}
// ||---||
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r) + OFFSET(c),
		     msr, mbr, mpr, 
		     s,   b,   p     );
}

// ---||-|
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(c),
		     mlr, mvr, msr, mrr, 
		     l,   v,   s,   r     );
}
// |--||-|
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		     mvr, msr, mrr, 
		     v,   s,   r     );
}

// --|-|-|
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(c),
		     mlr, mvr, mbr, mrr, 
		     l,   v,   b,   r     );
}
// |-|-|-|
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(c),
		     mvr, mbr, mrr, 
		     v,   b,   r     );
}

// -|--|-|
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(c),
		     mlr, msr, mbr, mrr, 
		     l,   s,   b,   r     );
}
// ||--|-|
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(c),
		     msr, mbr, mrr, 
		     s,   b,   r     );
}

// --||--|
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(c),
		     mlr, mvr, mpr, mrr, 
		     l,   v,   p,   r     );
}
// |-||--|
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(c),
		     mvr, mpr, mrr, 
		     v,   p,   r     );
}

// -|-|--|
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(c),
		     mlr, msr, mpr, mrr, 
		     l,   s,   p,   r     );
}
// ||-|--|
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(c),
		     msr, mpr, mrr, 
		     s,   p,   r     );
}

// -||---|
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(c),
		     mlr, mbr, mpr, mrr, 
		     l,   b,   p,   r     );
}
// |||---|
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(c),
		     mbr, mpr, mrr, 
		     b,   p,   r     );
}

// ---|||-
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r),
		     mlr, mvr, msr, mcr, 
		     l,   v,   s,   c     );
}
// |--|||-
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		     mvr, msr, mcr, 
		     v,   s,   c     );
}

// --|-||-
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r),
		     mlr, mvr, mbr, mcr, 
		     l,   v,   b,   c     );
}
// |-|-||-
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r),
		     mvr, mbr, mcr, 
		     v,   b,   c     );
}

// -|--||-
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r),
		     mlr, msr, mbr, mcr, 
		     l,   s,   b,   c     );
}
// ||--||-
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r),
		     msr, mbr, mcr, 
		     s,   b,   c     );
}

// --||-|-
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r),
		     mlr, mvr, mpr, mcr, 
		     l,   v,   p,   c     );
}
// |-||-|-
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r),
		     mvr, mpr, mcr, 
		     v,   p,   c     );
}

// -|-|-|-
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r),
		     mlr, msr, mpr, mcr, 
		     l,   s,   p,   c     );
}
// ||-|-|-
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r),
		     msr, mpr, mcr, 
		     s,   p,   c     );
}

// -||--|-
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r),
		     mlr, mbr, mpr, mcr, 
		     l,   b,   p,   c     );
}
// |||--|-
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r),
		     mbr, mpr, mcr, 
		     b,   p,   c     );
}

// --|||--
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p),
		     mlr, mvr, mrr, mcr, 
		     l,   v,   r,   c     );
}
// |-|||--
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p),
		     mvr, mrr, mcr, 
		     v,   r,   c     );
}

// -|-||--
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p),
		     mlr, msr, mrr, mcr, 
		     l,   s,   r,   c     );
}
// ||-||--
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p),
		     msr, mrr, mcr, 
		     s,   r,   c     );
}

// -||-|--
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p),
		     mlr, mbr, mrr, mcr, 
		     l,   b,   r,   c     );
}
// |||-|--
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p),
		     mbr, mrr, mcr, 
		     b,   r,   c     );
}

// -|||---
inline Tensor4View Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b),
		     mlr, mpr, mrr, mcr, 
		     l,   p,   r,   c     );
}
// ||||---
inline Tensor3View Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b),
		     mpr, mrr, mcr, 
		     p,   r,   c     );
}

// -||||--
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
		      mlr, mrr, mcr, 
		      l,   r,   c     );
}
// |||||--
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
		     mrr, mcr, 
		     r,   c     );
}

// -|||-|-
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
		      mlr, mpr, mcr, 
		      l,   p,   c     );
}
// ||||-|-
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
		     mpr, mcr, 
		     p,   c     );
}

// -||-||-
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
		      mlr, mbr, mcr, 
		      l,   b,   c     );
}
// |||-||-
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
		     mbr, mcr, 
		     b,   c     );
}

// -|-|||-
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		      mlr, msr, mcr, 
		      l,   s,   c     );
}
// ||-|||-
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		     msr, mcr, 
		     s,   c     );
}

// --||||-
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		      mlr, mvr, mcr, 
		      l,   v,   c     );
}
// |-||||-
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		     mvr, mcr, 
		     v,   c     );
}

// -|||--|
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
		      mlr, mpr, mrr, 
		      l,   p,   r     );
}
// ||||--|
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
		     mpr, mrr, 
		     p,   r     );
}

// -||-|-|
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
		      mlr, mbr, mrr, 
		      l,   b,   r     );
}
// |||-|-|
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
		     mbr, mrr, 
		     b,   r     );
}

// -|-||-|
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		      mlr, msr, mrr, 
		      l,   s,   r     );
}
// ||-||-|
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		     msr, mrr, 
		     s,   r     );
}

// --|||-|
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		      mlr, mvr, mrr, 
		      l,   v,   r     );
}
// |-|||-|
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		     mvr, mrr, 
		     v,   r     );
}

// -||--||
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
		      mlr, mbr, mpr, 
		      l,   b,   p     );
}
// |||--||
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
		     mbr, mpr, 
		     b,   p     );
}

// -|-|-||
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		      mlr, msr, mpr, 
		      l,   s,   p     );
}
// ||-|-||
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		     msr, mpr, 
		     s,   p     );
}

// --||-||
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		      mlr, mvr, mpr, 
		      l,   v,   p     );
}
// |-||-||
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		     mvr, mpr, 
		     v,   p     );
}

// -|--|||
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		      mlr, msr, mbr, 
		      l,   s,   b     );
}
// ||--|||
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     msr, mbr, 
		     s,   b     );
}

// --|-|||
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		      mlr, mvr, mbr, 
		      l,   v,   b     );
}
// |-|-|||
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mvr, mbr, 
		     v,   b     );
}

// ---||||
inline Tensor3View  Tensor7View::operator()
  ( const Range& l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  Tensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		      mlr, mvr, msr, 
		      l,   v,   s     );
}
// |--||||
inline MatrixView  Tensor7View::operator()
  ( Index        l,
    const Range& v, const Range& s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mvr, msr, 
		     v,   s     );
}

// -|||||-
inline MatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		     mlr, mcr,  
		     l,   c      );
}
// ||||||-
inline VectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, const Range& c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return  VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
		     mcr,  
		     c      );
}

// -||||-|
inline MatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		     mlr, mrr,  
		     l,   r      );
}
// |||||-|
inline VectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, const Range& r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return  VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
		     mrr,  
		     r      );
}

// -|||-||
inline MatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		     mlr, mpr,  
		     l,   p      );
}
// ||||-||
inline VectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    const Range& p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
		     mpr,  
		     p      );
}

// -||-|||
inline MatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mlr, mbr,  
		     l,   b      );
}
// |||-|||
inline VectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, const Range& b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mbr,  
		     b      );
}

// -|-||||
inline MatrixView  Tensor7View::operator()
  ( const Range& l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mlr, msr,  
		     l,   s      );
}
// ||-||||
inline VectorView  Tensor7View::operator()
  ( Index        l,
    Index        v, const Range& s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     msr,  
		     s      );
}

// --|||||
inline MatrixView  Tensor7View::operator()
  ( const Range& l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  MatrixView(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mlr, mvr,  
		     l,   v      );
}
// |-|||||
inline VectorView  Tensor7View::operator()
  ( Index        l,
    const Range& v, Index        s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
		     mvr,  
		     v      );
}

// -||||||
inline VectorView  Tensor7View::operator()
  ( const Range& l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c)
{
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return  VectorView( mdata +
		      OFFSET(v) + OFFSET(s) + OFFSET(b) +
		      OFFSET(p) + OFFSET(r) + OFFSET(c),
		      mlr,
		      l    );
}
// |||||||
inline Numeric&         Tensor7View::operator()
  ( Index        l,
    Index        v, Index        s, Index        b,
    Index        p, Index        r, Index        c) 
{
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return                *(mdata + OFFSET(l) +
			  OFFSET(v) + OFFSET(s) + OFFSET(b) +
			  OFFSET(p) + OFFSET(r) + OFFSET(c)    );
}


/** Return  iterator to sub-tensor. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator7D Tensor7View::begin() const
{
  return ConstTensor7View::begin();
}

/** Return const iterator behind last sub-tensor. */
inline ConstIterator7D Tensor7View::end() const
{
  return ConstTensor7View::end();
}

/** Return iterator to first sub-tensor. */
inline Iterator7D Tensor7View::begin()
{
  return Iterator7D( Tensor6View(mdata+mlr.mstart,
				 mvr,
				 msr,
				 mbr,
				 mpr,
				 mrr,
				 mcr),
		     mlr.mstride);
}

/** Return iterator behind last sub-tensor. */
inline Iterator7D Tensor7View::end()
{
  return Iterator7D( Tensor6View(mdata + mlr.mstart +
				(mlr.mextent)*mlr.mstride,
				 mvr,
				 msr,
				 mbr,
				 mpr,
				 mrr,
				 mcr),
		     mlr.mstride );
}

/** Assignment operator. This copies the data from another Tensor7View
    to this Tensor7View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor7View by
    setting its range. */
inline Tensor7View& Tensor7View::operator=(const ConstTensor7View& m)
{
  // Check that sizes are compatible:
  assert(mlr.mextent==m.mlr.mextent);
  assert(mvr.mextent==m.mvr.mextent);
  assert(msr.mextent==m.msr.mextent);
  assert(mbr.mextent==m.mbr.mextent);
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from Tensor7View to Tensor7View. This is a tricky
    one. The problem is that since Tensor7View is derived from
    ConstTensor7View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline Tensor7View& Tensor7View::operator=(const Tensor7View& m)
{
  // Check that sizes are compatible:
  assert(mlr.mextent==m.mlr.mextent);
  assert(mvr.mextent==m.mvr.mextent);
  assert(msr.mextent==m.msr.mextent);
  assert(mbr.mextent==m.mbr.mextent);
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Tensor7. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline Tensor7View& Tensor7View::operator=(const Tensor7& m)
{
  // Check that sizes are compatible:
  assert(mlr.mextent==m.mlr.mextent);
  assert(mvr.mextent==m.mvr.mextent);
  assert(msr.mextent==m.msr.mextent);
  assert(mbr.mextent==m.mbr.mextent);
  assert(mpr.mextent==m.mpr.mextent);
  assert(mrr.mextent==m.mrr.mextent);
  assert(mcr.mextent==m.mcr.mextent);

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor7View will set all elements to this
    value. */
inline Tensor7View& Tensor7View::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
inline Tensor7View& Tensor7View::operator*=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
inline Tensor7View& Tensor7View::operator/=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
inline Tensor7View& Tensor7View::operator+=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
inline Tensor7View& Tensor7View::operator-=(Numeric x)
{
  const Iterator7D ep=end();  
  for ( Iterator7D p=begin(); p!=ep ; ++p )
  {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor7. */
inline Tensor7View& Tensor7View::operator*=(const ConstTensor7View& x)
{
  assert( nlibraries() == x.nlibraries() );
  assert( nvitrines()  == x.nvitrines()  );
  assert( nshelves()   == x.nshelves()   );
  assert( nbooks()     == x.nbooks()     );
  assert( npages()     == x.npages()     );
  assert( nrows()      == x.nrows()      );
  assert( ncols()      == x.ncols()      );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p *= *xp;
    }
  return *this;
}

/** Element-vise division by another Tensor7. */
inline Tensor7View& Tensor7View::operator/=(const ConstTensor7View& x)
{
  assert( nlibraries() == x.nlibraries() );
  assert( nvitrines()  == x.nvitrines()  );
  assert( nshelves()   == x.nshelves()   );
  assert( nbooks()     == x.nbooks()     );
  assert( npages()     == x.npages()     );
  assert( nrows()      == x.nrows()      );
  assert( ncols()      == x.ncols()      );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p /= *xp;
    }
  return *this;
}

/** Element-vise addition of another Tensor7. */
inline Tensor7View& Tensor7View::operator+=(const ConstTensor7View& x)
{
  assert( nlibraries() == x.nlibraries() );
  assert( nvitrines()  == x.nvitrines()  );
  assert( nshelves()   == x.nshelves()   );
  assert( nbooks()     == x.nbooks()     );
  assert( npages()     == x.npages()     );
  assert( nrows()      == x.nrows()      );
  assert( ncols()      == x.ncols()      );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p += *xp;
    }
  return *this;
}

/** Element-vise subtraction of another Tensor7. */
inline Tensor7View& Tensor7View::operator-=(const ConstTensor7View& x)
{
  assert( nlibraries() == x.nlibraries() );
  assert( nvitrines()  == x.nvitrines()  );
  assert( nshelves()   == x.nshelves()   );
  assert( nbooks()     == x.nbooks()     );
  assert( npages()     == x.npages()     );
  assert( nrows()      == x.nrows()      );
  assert( ncols()      == x.ncols()      );
  ConstIterator7D  xp = x.begin();
  Iterator7D        p = begin();
  const Iterator7D ep = end();
  for ( ; p!=ep ; ++p,++xp )
    {
      *p -= *xp;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Tensor7. */
inline Tensor7View::Tensor7View() :
  ConstTensor7View()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor7 to initialize its
    own Tensor7View part. */
inline Tensor7View::Tensor7View(Numeric *data,
				const Range& l,
				const Range& v,
				const Range& s,
				const Range& b,
				const Range& p,
				const Range& r,
				const Range& c) :
  ConstTensor7View(data, l, v, s, b, p, r, c)
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
inline Tensor7View::Tensor7View(Numeric *data,
				const Range& pl,
				const Range& pv,
				const Range& ps,
				const Range& pb,
				const Range& pp,
				const Range& pr,
				const Range& pc,
				const Range& nl,
				const Range& nv,
				const Range& ns,
				const Range& nb,
				const Range& np,
				const Range& nr,
				const Range& nc) :
  ConstTensor7View(data,pl,pv,ps,pb,pp,pr,pc,nl,nv,ns,nb,np,nr,nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
inline void copy(ConstIterator7D origin,
                 const ConstIterator7D& end,
                 Iterator7D target)
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
                 Iterator7D target,
                 const Iterator7D& end)
{
  for ( ; target!=end ; ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy(x,target->begin(),target->end());
    }
}


// Functions for Tensor7:
// ---------------------

/** Default constructor. */
inline Tensor7::Tensor7() :
  Tensor7View::Tensor7View()
{
  // Nothing to do here. However, note that the default constructor
  // for Tensor7View has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized. 
}

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
inline Tensor7::Tensor7(Index l,
			Index v, Index s, Index b,
                        Index p, Index r, Index c) :
  Tensor7View( new Numeric[l*v*s*b*p*r*c],
	       Range(0,l,v*s*b*p*r*c),
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
inline Tensor7::Tensor7(Index l,
			Index v, Index s, Index b,
                        Index p, Index r, Index c, Numeric fill) :
  Tensor7View( new Numeric[l*v*s*b*p*r*c],
	       Range(0,l,v*s*b*p*r*c),
	       Range(0,v,s*b*p*r*c),
	       Range(0,s,b*p*r*c),
	       Range(0,b,p*r*c),
	       Range(0,p,r*c),
	       Range(0,r,c),
	       Range(0,c))
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  const Numeric *stop = mdata+l*v*s*b*p*r*c;
  for ( Numeric *x=mdata; x<stop; ++x )
    *x = fill;
}

/** Copy constructor from Tensor7View. This automatically sets the size
    and copies the data. */
inline Tensor7::Tensor7(const ConstTensor7View& m) :
  Tensor7View( new Numeric[m.nlibraries()*m.nvitrines()*m.nshelves()*m.nbooks(),
			   m.npages()*m.nrows()*m.ncols()],
	       Range( 0, m.nlibraries(), m.nvitrines()*m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.nvitrines(), m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.nshelves(), m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.nbooks(), m.npages()*m.nrows()*m.ncols() ),
	       Range( 0, m.npages(), m.nrows()*m.ncols() ),
	       Range( 0, m.nrows(), m.ncols() ),
	       Range( 0, m.ncols() ) )
{
  copy(m.begin(),m.end(),begin());
}

/** Copy constructor from Tensor7. This automatically sets the size
    and copies the data. */
inline Tensor7::Tensor7(const Tensor7& m) :
  Tensor7View( new Numeric[m.nlibraries()*m.nvitrines()*m.nshelves()*m.nbooks(),
			   m.npages()*m.nrows()*m.ncols()],
	       Range( 0, m.nlibraries(), m.nvitrines()*m.nshelves()*m.nbooks()*m.npages()*m.nrows()*m.ncols() ),
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
inline Tensor7& Tensor7::operator=(const Tensor7& m)
{
  //  cout << "Tensor7 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "             n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid tensor, so we just
  // have to check one.
  if ( 0 == mcr.mextent )
    {
      // Adjust if previously empty.
      resize( m.mlr.mextent,
	      m.mvr.mextent, m.msr.mextent, m.mbr.mextent,
	      m.mpr.mextent, m.mrr.mextent, m.mcr.mextent ); 
    }
  else
    {
      // Check that sizes are compatible:
      assert( mlr.mextent==m.mlr.mextent );
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
inline Tensor7& Tensor7::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values. */
inline void Tensor7::resize( Index l,
			     Index v, Index s, Index b,
			     Index p, Index r, Index c)
{
  assert( 0<=l );
  assert( 0<=v );
  assert( 0<=s );
  assert( 0<=b );
  assert( 0<=p );
  assert( 0<=r );
  assert( 0<=c );

  if ( mlr.mextent!=l ||
       mvr.mextent!=v ||
       msr.mextent!=s ||
       mbr.mextent!=b ||
       mpr.mextent!=p ||
       mrr.mextent!=r ||
       mcr.mextent!=c )
    {
      delete mdata;
      mdata = new Numeric[l*v*s*b*p*r*c];

      mlr.mstart = 0;
      mlr.mextent = l;
      mlr.mstride = v*s*b*p*r*c;

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

/** Destructor for Tensor7. This is important, since Tensor7 uses new to
    allocate storage. */
inline Tensor7::~Tensor7()
{
//   cout << "Destroying a Tensor7:\n"
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
inline void transform( Tensor7View y,
                       double (&my_func)(double),
                       ConstTensor7View x )
{
  // Check dimensions:
  assert( y.nlibraries() == x.nlibraries() );
  assert( y.nvitrines()  == x.nvitrines()  );
  assert( y.nshelves()   == x.nshelves()   );
  assert( y.nbooks()     == x.nbooks()     );
  assert( y.npages() 	 == x.npages()     );
  assert( y.nrows()  	 == x.nrows()      );
  assert( y.ncols()  	 == x.ncols()      );

  const ConstIterator7D xe = x.end();
  ConstIterator7D       xi = x.begin();
  Iterator7D            yi = y.begin();
  for ( ; xi!=xe; ++xi, ++yi )
    {
      // Use the transform function of lower dimensional tensors
      // recursively:
      transform(*yi,my_func,*xi);
    }
}

/** Max function, tensor version. */
inline Numeric max(const ConstTensor7View& x)
{
  const ConstIterator7D xe = x.end();
  ConstIterator7D       xi = x.begin();

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
inline Numeric min(const ConstTensor7View& x)
{
  const ConstIterator7D xe = x.end();
  ConstIterator7D       xi = x.begin();

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



#endif    // matpackVII_h
