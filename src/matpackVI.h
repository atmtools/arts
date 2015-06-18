/* Copyright (C) 2001-2012 Stefan Buehler <sbuehler@ltu.se>

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

#include "matpackV.h"

#define CHECK(x)   assert( 0 <= x ); assert( x < m ## x ## r.mextent )
#define OFFSET(x)  m ## x ## r.mstart + x * m ## x ## r.mstride

/** The outermost iterator class for rank 6 tensors. This takes into
    account the defined strided. */
class Iterator6D {
public:
  // Constructors:
  /** Default constructor. */
  Iterator6D() : msv(), mstride(0) { /* Nothing to do here. */ }

  /** Explicit constructor. */
  Iterator6D(const Tensor5View& x, Index stride) : msv(x), mstride(stride)
      { /* Nothing to do here. */ }

  // Operators:
  /** Prefix increment operator. */
  Iterator6D& operator++() { msv.mdata += mstride; return *this; }

  /** Not equal operator, needed for algorithms like copy.
      FIXME: Is it really necessary to have such a complicated check
      here? It could be sufficient to just test
      msv.mdata!=other.msv.mdata. */
  bool operator!=(const Iterator6D& other) const
    { if ( msv.mdata +
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
  Tensor5View* operator->() { return &msv; }

  /** Dereferencing. */
  Tensor5View& operator*() { return msv; }

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
  /** Default constructor. */
  ConstIterator6D() : msv(), mstride(0) { /* Nothing to do here. */ }

  /** Explicit constructor. */
  ConstIterator6D(const ConstTensor5View& x, Index stride) : msv(x), mstride(stride)
      { /* Nothing to do here. */ }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator6D& operator++() { msv.mdata += mstride; return *this; }

  /** Not equal operator, needed for algorithms like copy.
      FIXME: Is it really necessaary to have such a complicated check
      here? It could be sufficient to just test
      msv.mdata!=other.msv.mdata. */
  bool operator!=(const ConstIterator6D& other) const
    { if ( msv.mdata +
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
  const ConstTensor5View* operator->() const { return &msv; }

  /** Dereferencing. */
  const ConstTensor5View& operator*() const { return msv; }

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
  bool empty() const;
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
  Numeric operator() ( Index        v, Index        s, Index        b,
                       Index        p, Index        r, Index        c) const
      { CHECK(v);
        CHECK(s);
        CHECK(b);
        CHECK(p);
        CHECK(r);
        CHECK(c);
        return get(v, s, b, p, r, c);
      }

  /** Get element implementation without assertions. */
  Numeric get( Index        v, Index        s, Index        b,
               Index        p, Index        r, Index        c) const
      {
        return                *(mdata +
                                OFFSET(v) + OFFSET(s) + OFFSET(b) +
                                OFFSET(p) + OFFSET(r) + OFFSET(c)    );
      }


  // Functions returning iterators:
  ConstIterator6D begin() const;
  ConstIterator6D end() const;

  // Destructor:
  virtual ~ConstTensor6View() {}
  
  // Friends:
  friend class ConstIterator7D;
  friend class Tensor6View;
  friend class ConstTensor7View;

  // Special constructor to make a Tensor6 view of a Tensor5.
  ConstTensor6View(const ConstTensor5View& a);

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
  Numeric operator() ( Index        v, Index        s, Index        b,
                       Index        p, Index        r, Index        c) const
      { return ConstTensor6View::operator()(v,s,b,p,r,c);  }

  /** Get element implementation without assertions. */
  Numeric get( Index        v, Index        s, Index        b,
               Index        p, Index        r, Index        c) const
      { return ConstTensor6View::get(v,s,b,p,r,c);  }

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
  Numeric& operator() ( Index        v, Index        s, Index        b,
                        Index        p, Index        r, Index        c)
    {   CHECK(v);
        CHECK(s);
        CHECK(b);
        CHECK(p);
        CHECK(r);
        CHECK(c);
        return get(v, s, b, p, r, c);
      }

  /** Get element implementation without assertions. */
  Numeric& get( Index        v, Index        s, Index        b,
                Index        p, Index        r, Index        c)
    {
        return                *(mdata +
                                OFFSET(v) + OFFSET(s) + OFFSET(b) +
                                OFFSET(p) + OFFSET(r) + OFFSET(c)    );
    }

  // Conversion to a plain C-array
  const Numeric *get_c_array() const;
  Numeric *get_c_array();

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

  // Destructor:
  virtual ~Tensor6View() {}
  
  // Friends:
  friend class Iterator7D;
  friend class Tensor7View;

  // Special constructor to make a Tensor6 view of a Tensor5.
  Tensor6View(const Tensor5View& a);

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
  Tensor6& operator=(Tensor6 x);
  Tensor6& operator=(Numeric x);

  // Resize function:
  void resize(Index        v, Index        s, Index        b,
              Index        p, Index        r, Index        c);

  // Swap function:
  friend void swap(Tensor6& t1, Tensor6& t2);

  // Destructor:
  virtual ~Tensor6();
};


// Function declarations:
// ----------------------

void copy(ConstIterator6D origin,
          const ConstIterator6D& end,
          Iterator6D target);

void copy(Numeric x,
          Iterator6D target,
          const Iterator6D& end);

void transform( Tensor6View y,
                double (&my_func)(double),
                ConstTensor6View x );

Numeric max(const ConstTensor6View& x);

Numeric min(const ConstTensor6View& x);

std::ostream& operator<<(std::ostream& os, const ConstTensor6View& v);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor6view_get_elem (Tensor6View& tv, Index v, Index s, Index b,
                                    Index p, Index r, Index c);

#endif
////////////////////////////////

#endif    // matpackVI_h
