
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
   Implementation of Tensors of Rank 7.

   Dimensions are called: library, vitrine, shelf, book, page, row, column.
   or short:              l,       v,       s,     b,    p,    r,   c
  
   \author Stefan Buehler
   \date   2001-11-22
*/

#ifndef matpackVII_h
#define matpackVII_h

#include <utility>

#include "matpackVI.h"
#include "matpack_concepts.h"

/** The outermost iterator class for rank 7 tensors. This takes into
    account the defined strided. */
class Iterator7D {
 public:
  // Constructors:
  /** Default constructor. */
  Iterator7D() = default;

  /** Explicit constructor. */
  Iterator7D(const Tensor6View& x, Index stride)
      : msv(x), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  Iterator7D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy.
      FIXME: Is it really necessary to have such a complicated check
      here? It could be sufficient to just test
      msv.mdata!=other.msv.mdata. */
  bool operator!=(const Iterator7D& other) const {
    if (msv.mdata + msv.mvr.mstart + msv.msr.mstart + msv.mbr.mstart +
            msv.mpr.mstart + msv.mrr.mstart + msv.mcr.mstart !=
        other.msv.mdata + other.msv.mvr.mstart + other.msv.msr.mstart +
            other.msv.mbr.mstart + other.msv.mpr.mstart + other.msv.mrr.mstart +
            other.msv.mcr.mstart)
      return true;
    return false;
  }

  /** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
  Tensor6View* operator->() { return &msv; }

  /** Dereferencing. */
  Tensor6View& operator*() { return msv; }

 private:
  /** Current position. */
  Tensor6View msv;
  /** Stride. */
  Index mstride{0};
};

/** Const version of Iterator7D. */
class ConstIterator7D {
 public:
  // Constructors:
  /** Default constructor. */
  ConstIterator7D() = default;

  /** Explicit constructor. */
  ConstIterator7D(ConstTensor6View x, Index stride)
      : msv(std::move(x)), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator7D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy. 
      FIXME: Is it really necessary to have such a complicated check
      here? It could be sufficient to just test
      msv.mdata!=other.msv.mdata. */
  bool operator!=(const ConstIterator7D& other) const {
    if (msv.mdata + msv.mvr.mstart + msv.msr.mstart + msv.mbr.mstart +
            msv.mpr.mstart + msv.mrr.mstart + msv.mcr.mstart !=
        other.msv.mdata + other.msv.mvr.mstart + other.msv.msr.mstart +
            other.msv.mbr.mstart + other.msv.mpr.mstart + other.msv.mrr.mstart +
            other.msv.mcr.mstart)
      return true;
    return false;
  }

  /** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
  const ConstTensor6View* operator->() const { return &msv; }

  /** Dereferencing. */
  const ConstTensor6View& operator*() const { return msv; }

 private:
  /** Current position. */
  ConstTensor6View msv;
  /** Stride. */
  Index mstride{0};
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
  static constexpr bool matpack_type{true};
  
  constexpr ConstTensor7View(const ConstTensor7View&) = default;
  constexpr ConstTensor7View(ConstTensor7View&&) = default;
  ConstTensor7View& operator=(const ConstTensor7View&) = default;
  ConstTensor7View& operator=(ConstTensor7View&&) = default;

  // Member functions:
  [[nodiscard]] Index nlibraries() const noexcept { return mlr.mextent; }
  [[nodiscard]] Index nvitrines() const noexcept { return mvr.mextent; }
  [[nodiscard]] Index nshelves() const noexcept { return msr.mextent; }
  [[nodiscard]] Index nbooks() const noexcept { return mbr.mextent; }
  [[nodiscard]] Index npages() const noexcept { return mpr.mextent; }
  [[nodiscard]] Index nrows() const noexcept { return mrr.mextent; }
  [[nodiscard]] Index ncols() const noexcept { return mcr.mextent; }

  // Total size
  [[nodiscard]] Index size() const noexcept {
    return nlibraries() * nvitrines() * nshelves() * nbooks() * npages() *
           nrows() * ncols();
  }
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /*! Returns the shape as an array (to allow templates to just look for shape on different matpack objects) */
  [[nodiscard]] Shape<7> shape() const {
    return {nlibraries(),
            nvitrines(),
            nshelves(),
            nbooks(),
            npages(),
            nrows(),
            ncols()};
  }

  // Const index operators:

  // Result 7D (1 combination)
  // -------
  ConstTensor7View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 6D (7 combinations)
  // ------|
  ConstTensor6View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // -----|-
  ConstTensor6View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ----|--
  ConstTensor6View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // ---|---
  ConstTensor6View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // --|----
  ConstTensor6View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // -|-----
  ConstTensor6View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |------
  ConstTensor6View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // ----|-|
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // ---|--|
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // --|---|
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // -|----|
  ConstTensor5View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // |-----|
  ConstTensor5View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ----||-
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // ---|-|-
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // --|--|-
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // -|---|-
  ConstTensor5View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // |----|-
  ConstTensor5View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ---||--
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // --|-|--
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -|--|--
  ConstTensor5View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |---|--
  ConstTensor5View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // --||---
  ConstTensor5View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // -|-|---
  ConstTensor5View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |--|---
  ConstTensor5View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // -||----
  ConstTensor5View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |-|----
  ConstTensor5View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // ||-----
  ConstTensor5View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              Index c) const;
  // ---|-||
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // --|--||
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // -|---||
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // |----||
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // ---||-|
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // --|-|-|
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // -|--|-|
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // |---|-|
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // --||--|
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // -|-|--|
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // |--|--|
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // -||---|
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // |-|---|
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ||----|
  ConstTensor4View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ---|||-
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // --|-||-
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // -|--||-
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // |---||-
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // --||-|-
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // -|-|-|-
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // |--|-|-
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // -||--|-
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // |-|--|-
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ||---|-
  ConstTensor4View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // --|||--
  ConstTensor4View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -|-||--
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |--||--
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -||-|--
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |-|-|--
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // ||--|--
  ConstTensor4View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -|||---
  ConstTensor4View operator()(const Range& l,
                              Index v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |-||---
  ConstTensor4View operator()(Index l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // ||-|---
  ConstTensor4View operator()(Index l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |||----
  ConstTensor4View operator()(Index l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  ConstTensor3View operator()(Index l,
                              Index v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |||-|--
  ConstTensor3View operator()(Index l,
                              Index v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // ||-||--
  ConstTensor3View operator()(Index l,
                              Index v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |-|||--
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              Index s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -||||--
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              Index s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |||--|-
  ConstTensor3View operator()(Index l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ||-|-|-
  ConstTensor3View operator()(Index l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // |-||-|-
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // -|||-|-
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              Index s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ||--||-
  ConstTensor3View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // |-|-||-
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // -||-||-
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // |--|||-
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // -|-|||-
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // --||||-
  ConstTensor3View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // |||---|
  ConstTensor3View operator()(Index l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ||-|--|
  ConstTensor3View operator()(Index l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // |-||--|
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // -|||--|
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ||--|-|
  ConstTensor3View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // |-|-|-|
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // -||-|-|
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // |--||-|
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // -|-||-|
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // --|||-|
  ConstTensor3View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // ||---||
  ConstTensor3View operator()(Index l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // |-|--||
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // -||--||
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // |--|-||
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // -|-|-||
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // --||-||
  ConstTensor3View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // |---|||
  ConstTensor3View operator()(Index l,
                              const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              Index c) const;
  // -|--|||
  ConstTensor3View operator()(const Range& l,
                              Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              Index c) const;
  // --|-|||
  ConstTensor3View operator()(const Range& l,
                              const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              Index r,
                              Index c) const;
  // ---||||
  ConstTensor3View operator()(const Range& l,
                              const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              Index r,
                              Index c) const;

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  ConstMatrixView operator()(Index l,
                             Index v,
                             Index s,
                             Index b,
                             Index p,
                             const Range& r,
                             const Range& c) const;
  // ||||-|-
  ConstMatrixView operator()(Index l,
                             Index v,
                             Index s,
                             Index b,
                             const Range& p,
                             Index r,
                             const Range& c) const;
  // |||-||-
  ConstMatrixView operator()(Index l,
                             Index v,
                             Index s,
                             const Range& b,
                             Index p,
                             Index r,
                             const Range& c) const;
  // ||-|||-
  ConstMatrixView operator()(Index l,
                             Index v,
                             const Range& s,
                             Index b,
                             Index p,
                             Index r,
                             const Range& c) const;
  // |-||||-
  ConstMatrixView operator()(Index l,
                             const Range& v,
                             Index s,
                             Index b,
                             Index p,
                             Index r,
                             const Range& c) const;
  // -|||||-
  ConstMatrixView operator()(const Range& l,
                             Index v,
                             Index s,
                             Index b,
                             Index p,
                             Index r,
                             const Range& c) const;
  // ||||--|
  ConstMatrixView operator()(Index l,
                             Index v,
                             Index s,
                             Index b,
                             const Range& p,
                             const Range& r,
                             Index c) const;
  // |||-|-|
  ConstMatrixView operator()(Index l,
                             Index v,
                             Index s,
                             const Range& b,
                             Index p,
                             const Range& r,
                             Index c) const;
  // ||-||-|
  ConstMatrixView operator()(Index l,
                             Index v,
                             const Range& s,
                             Index b,
                             Index p,
                             const Range& r,
                             Index c) const;
  // |-|||-|
  ConstMatrixView operator()(Index l,
                             const Range& v,
                             Index s,
                             Index b,
                             Index p,
                             const Range& r,
                             Index c) const;
  // -||||-|
  ConstMatrixView operator()(const Range& l,
                             Index v,
                             Index s,
                             Index b,
                             Index p,
                             const Range& r,
                             Index c) const;
  // |||--||
  ConstMatrixView operator()(Index l,
                             Index v,
                             Index s,
                             const Range& b,
                             const Range& p,
                             Index r,
                             Index c) const;
  // ||-|-||
  ConstMatrixView operator()(Index l,
                             Index v,
                             const Range& s,
                             Index b,
                             const Range& p,
                             Index r,
                             Index c) const;
  // |-||-||
  ConstMatrixView operator()(Index l,
                             const Range& v,
                             Index s,
                             Index b,
                             const Range& p,
                             Index r,
                             Index c) const;
  // -|||-||
  ConstMatrixView operator()(const Range& l,
                             Index v,
                             Index s,
                             Index b,
                             const Range& p,
                             Index r,
                             Index c) const;
  // ||--|||
  ConstMatrixView operator()(Index l,
                             Index v,
                             const Range& s,
                             const Range& b,
                             Index p,
                             Index r,
                             Index c) const;
  // |-|-|||
  ConstMatrixView operator()(Index l,
                             const Range& v,
                             Index s,
                             const Range& b,
                             Index p,
                             Index r,
                             Index c) const;
  // -||-|||
  ConstMatrixView operator()(const Range& l,
                             Index v,
                             Index s,
                             const Range& b,
                             Index p,
                             Index r,
                             Index c) const;
  // |--||||
  ConstMatrixView operator()(Index l,
                             const Range& v,
                             const Range& s,
                             Index b,
                             Index p,
                             Index r,
                             Index c) const;
  // -|-||||
  ConstMatrixView operator()(const Range& l,
                             Index v,
                             const Range& s,
                             Index b,
                             Index p,
                             Index r,
                             Index c) const;
  // --|||||
  ConstMatrixView operator()(const Range& l,
                             const Range& v,
                             Index s,
                             Index b,
                             Index p,
                             Index r,
                             Index c) const;

  // Result 1D (7 combinations)
  // ||||||-
  ConstVectorView operator()(Index l,
                             Index v,
                             Index s,
                             Index b,
                             Index p,
                             Index r,
                             const Range& c) const;
  // |||||-|
  ConstVectorView operator()(Index l,
                             Index v,
                             Index s,
                             Index b,
                             Index p,
                             const Range& r,
                             Index c) const;
  // ||||-||
  ConstVectorView operator()(Index l,
                             Index v,
                             Index s,
                             Index b,
                             const Range& p,
                             Index r,
                             Index c) const;
  // |||-|||
  ConstVectorView operator()(Index l,
                             Index v,
                             Index s,
                             const Range& b,
                             Index p,
                             Index r,
                             Index c) const;
  // ||-||||
  ConstVectorView operator()(Index l,
                             Index v,
                             const Range& s,
                             Index b,
                             Index p,
                             Index r,
                             Index c) const;
  // |-|||||
  ConstVectorView operator()(Index l,
                             const Range& v,
                             Index s,
                             Index b,
                             Index p,
                             Index r,
                             Index c) const;
  // -||||||
  ConstVectorView operator()(const Range& l,
                             Index v,
                             Index s,
                             Index b,
                             Index p,
                             Index r,
                             Index c) const;

  // Result scalar (1 combination)
  // |||||||
  Numeric operator()(
      Index l, Index v, Index s, Index b, Index p, Index r, Index c) const {
    CHECK(l);
    CHECK(v);
    CHECK(s);
    CHECK(b);
    CHECK(p);
    CHECK(r);
    CHECK(c);
    return get(l, v, s, b, p, r, c);
  }

  /** Get element implementation without assertions. */
  [[nodiscard]] Numeric get(
      Index l, Index v, Index s, Index b, Index p, Index r, Index c) const {
    return *(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) +
             OFFSET(r) + OFFSET(c));
  }

  // Functions returning iterators:
  [[nodiscard]] ConstIterator7D begin() const;
  [[nodiscard]] ConstIterator7D end() const;

  //! Destructor.
  virtual ~ConstTensor7View() = default;

  // Friends:
  friend class Tensor7View;

  friend std::ostream& operator<<(std::ostream& os, const ConstTensor7View& v);

  // Special constructor to make a Tensor7 view of a Tensor6.
  ConstTensor7View(const ConstTensor6View& a);

 protected:
  // Constructors:
  ConstTensor7View() = default;
  ConstTensor7View(Numeric* data,
                   const Range& l,
                   const Range& v,
                   const Range& s,
                   const Range& b,
                   const Range& p,
                   const Range& r,
                   const Range& c);
  ConstTensor7View(Numeric* data,
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
                   const Range& nc);

  // Data members:
  // -------------
  /** The library range of mdata that is actually used. */
  Range mlr{0, 0, 1};
  /** The vitrine range of mdata that is actually used. */
  Range mvr{0, 0, 1};
  /** The shelf range of mdata that is actually used. */
  Range msr{0, 0, 1};
  /** The book range of mdata that is actually used. */
  Range mbr{0, 0, 1};
  /** The page range of mdata that is actually used. */
  Range mpr{0, 0, 1};
  /** The row range of mdata that is actually used. */
  Range mrr{0, 0, 1};
  /** The column range of mdata that is actually used. */
  Range mcr{0, 0, 1};
  /** Pointer to the plain C array that holds the data */
  Numeric* mdata{nullptr};
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
  // Make const methods visible from base class
  using ConstTensor7View::begin;
  using ConstTensor7View::end;
  using ConstTensor7View::operator();
  using ConstTensor7View::get;

  constexpr Tensor7View(const Tensor7View&) = default;

  // Non-const index operators:

  // Result 7D (1 combination)
  // -------
  Tensor7View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 6D (7 combinations)
  // ------|
  Tensor6View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // -----|-
  Tensor6View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ----|--
  Tensor6View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // ---|---
  Tensor6View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // --|----
  Tensor6View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // -|-----
  Tensor6View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |------
  Tensor6View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // ----|-|
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // ---|--|
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // --|---|
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // -|----|
  Tensor5View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // |-----|
  Tensor5View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ----||-
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // ---|-|-
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // --|--|-
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // -|---|-
  Tensor5View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // |----|-
  Tensor5View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ---||--
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // --|-|--
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -|--|--
  Tensor5View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |---|--
  Tensor5View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // --||---
  Tensor5View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // -|-|---
  Tensor5View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |--|---
  Tensor5View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // -||----
  Tensor5View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |-|----
  Tensor5View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // ||-----
  Tensor5View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         Index c);
  // ---|-||
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         Index c);
  // --|--||
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // -|---||
  Tensor4View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // |----||
  Tensor4View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // ---||-|
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         Index c);
  // --|-|-|
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // -|--|-|
  Tensor4View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // |---|-|
  Tensor4View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // --||--|
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // -|-|--|
  Tensor4View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // |--|--|
  Tensor4View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // -||---|
  Tensor4View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // |-|---|
  Tensor4View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ||----|
  Tensor4View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ---|||-
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         Index r,
                         const Range& c);
  // --|-||-
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // -|--||-
  Tensor4View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // |---||-
  Tensor4View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // --||-|-
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // -|-|-|-
  Tensor4View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // |--|-|-
  Tensor4View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // -||--|-
  Tensor4View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // |-|--|-
  Tensor4View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ||---|-
  Tensor4View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // --|||--
  Tensor4View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -|-||--
  Tensor4View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |--||--
  Tensor4View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -||-|--
  Tensor4View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |-|-|--
  Tensor4View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // ||--|--
  Tensor4View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -|||---
  Tensor4View operator()(const Range& l,
                         Index v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |-||---
  Tensor4View operator()(Index l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // ||-|---
  Tensor4View operator()(Index l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |||----
  Tensor4View operator()(Index l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  Tensor3View operator()(Index l,
                         Index v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |||-|--
  Tensor3View operator()(Index l,
                         Index v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // ||-||--
  Tensor3View operator()(Index l,
                         Index v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |-|||--
  Tensor3View operator()(Index l,
                         const Range& v,
                         Index s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -||||--
  Tensor3View operator()(const Range& l,
                         Index v,
                         Index s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |||--|-
  Tensor3View operator()(Index l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ||-|-|-
  Tensor3View operator()(Index l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // |-||-|-
  Tensor3View operator()(Index l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // -|||-|-
  Tensor3View operator()(const Range& l,
                         Index v,
                         Index s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ||--||-
  Tensor3View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // |-|-||-
  Tensor3View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // -||-||-
  Tensor3View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // |--|||-
  Tensor3View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         Index r,
                         const Range& c);
  // -|-|||-
  Tensor3View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         Index p,
                         Index r,
                         const Range& c);
  // --||||-
  Tensor3View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         Index p,
                         Index r,
                         const Range& c);
  // |||---|
  Tensor3View operator()(Index l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ||-|--|
  Tensor3View operator()(Index l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // |-||--|
  Tensor3View operator()(Index l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // -|||--|
  Tensor3View operator()(const Range& l,
                         Index v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ||--|-|
  Tensor3View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // |-|-|-|
  Tensor3View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // -||-|-|
  Tensor3View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // |--||-|
  Tensor3View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         Index c);
  // -|-||-|
  Tensor3View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         Index c);
  // --|||-|
  Tensor3View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         Index p,
                         const Range& r,
                         Index c);
  // ||---||
  Tensor3View operator()(Index l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // |-|--||
  Tensor3View operator()(Index l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // -||--||
  Tensor3View operator()(const Range& l,
                         Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // |--|-||
  Tensor3View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         Index c);
  // -|-|-||
  Tensor3View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         Index c);
  // --||-||
  Tensor3View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         Index r,
                         Index c);
  // |---|||
  Tensor3View operator()(Index l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         Index c);
  // -|--|||
  Tensor3View operator()(const Range& l,
                         Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         Index c);
  // --|-|||
  Tensor3View operator()(const Range& l,
                         const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         Index r,
                         Index c);
  // ---||||
  Tensor3View operator()(const Range& l,
                         const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         Index r,
                         Index c);

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  MatrixView operator()(Index l,
                        Index v,
                        Index s,
                        Index b,
                        Index p,
                        const Range& r,
                        const Range& c);
  // ||||-|-
  MatrixView operator()(Index l,
                        Index v,
                        Index s,
                        Index b,
                        const Range& p,
                        Index r,
                        const Range& c);
  // |||-||-
  MatrixView operator()(Index l,
                        Index v,
                        Index s,
                        const Range& b,
                        Index p,
                        Index r,
                        const Range& c);
  // ||-|||-
  MatrixView operator()(Index l,
                        Index v,
                        const Range& s,
                        Index b,
                        Index p,
                        Index r,
                        const Range& c);
  // |-||||-
  MatrixView operator()(Index l,
                        const Range& v,
                        Index s,
                        Index b,
                        Index p,
                        Index r,
                        const Range& c);
  // -|||||-
  MatrixView operator()(const Range& l,
                        Index v,
                        Index s,
                        Index b,
                        Index p,
                        Index r,
                        const Range& c);
  // ||||--|
  MatrixView operator()(Index l,
                        Index v,
                        Index s,
                        Index b,
                        const Range& p,
                        const Range& r,
                        Index c);
  // |||-|-|
  MatrixView operator()(Index l,
                        Index v,
                        Index s,
                        const Range& b,
                        Index p,
                        const Range& r,
                        Index c);
  // ||-||-|
  MatrixView operator()(Index l,
                        Index v,
                        const Range& s,
                        Index b,
                        Index p,
                        const Range& r,
                        Index c);
  // |-|||-|
  MatrixView operator()(Index l,
                        const Range& v,
                        Index s,
                        Index b,
                        Index p,
                        const Range& r,
                        Index c);
  // -||||-|
  MatrixView operator()(const Range& l,
                        Index v,
                        Index s,
                        Index b,
                        Index p,
                        const Range& r,
                        Index c);
  // |||--||
  MatrixView operator()(Index l,
                        Index v,
                        Index s,
                        const Range& b,
                        const Range& p,
                        Index r,
                        Index c);
  // ||-|-||
  MatrixView operator()(Index l,
                        Index v,
                        const Range& s,
                        Index b,
                        const Range& p,
                        Index r,
                        Index c);
  // |-||-||
  MatrixView operator()(Index l,
                        const Range& v,
                        Index s,
                        Index b,
                        const Range& p,
                        Index r,
                        Index c);
  // -|||-||
  MatrixView operator()(const Range& l,
                        Index v,
                        Index s,
                        Index b,
                        const Range& p,
                        Index r,
                        Index c);
  // ||--|||
  MatrixView operator()(Index l,
                        Index v,
                        const Range& s,
                        const Range& b,
                        Index p,
                        Index r,
                        Index c);
  // |-|-|||
  MatrixView operator()(Index l,
                        const Range& v,
                        Index s,
                        const Range& b,
                        Index p,
                        Index r,
                        Index c);
  // -||-|||
  MatrixView operator()(const Range& l,
                        Index v,
                        Index s,
                        const Range& b,
                        Index p,
                        Index r,
                        Index c);
  // |--||||
  MatrixView operator()(Index l,
                        const Range& v,
                        const Range& s,
                        Index b,
                        Index p,
                        Index r,
                        Index c);
  // -|-||||
  MatrixView operator()(const Range& l,
                        Index v,
                        const Range& s,
                        Index b,
                        Index p,
                        Index r,
                        Index c);
  // --|||||
  MatrixView operator()(const Range& l,
                        const Range& v,
                        Index s,
                        Index b,
                        Index p,
                        Index r,
                        Index c);

  // Result 1D (7 combinations)
  // ||||||-
  VectorView operator()(
      Index l, Index v, Index s, Index b, Index p, Index r, const Range& c);
  // |||||-|
  VectorView operator()(
      Index l, Index v, Index s, Index b, Index p, const Range& r, Index c);
  // ||||-||
  VectorView operator()(
      Index l, Index v, Index s, Index b, const Range& p, Index r, Index c);
  // |||-|||
  VectorView operator()(
      Index l, Index v, Index s, const Range& b, Index p, Index r, Index c);
  // ||-||||
  VectorView operator()(
      Index l, Index v, const Range& s, Index b, Index p, Index r, Index c);
  // |-|||||
  VectorView operator()(
      Index l, const Range& v, Index s, Index b, Index p, Index r, Index c);
  // -||||||
  VectorView operator()(
      const Range& l, Index v, Index s, Index b, Index p, Index r, Index c);

#define GETFUN(l, v, s, b, p, r, c)                                     \
  *(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + \
    OFFSET(r) + OFFSET(c))
  // Result scalar (1 combination)
  // |||||||
  Numeric& operator()(
      Index l, Index v, Index s, Index b, Index p, Index r, Index c) {
    CHECK(l);
    CHECK(v);
    CHECK(s);
    CHECK(b);
    CHECK(p);
    CHECK(r);
    CHECK(c);
    return GETFUN(l, v, s, b, p, r, c);
  }

  /** Get element implementation without assertions. */
  Numeric& get(Index l, Index v, Index s, Index b, Index p, Index r, Index c) {
    return GETFUN(l, v, s, b, p, r, c);
  }
#undef GETFUN

  // Conversion to a plain C-array
  [[nodiscard]] const Numeric* get_c_array() const ARTS_NOEXCEPT;
  Numeric* get_c_array() ARTS_NOEXCEPT;

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

  //! Destructor.
  ~Tensor7View() override = default;

  // Friends:

  // Special constructor to make a Tensor7 view of a Tensor6.
  Tensor7View(const Tensor6View& a);

 protected:
  // Constructors:
  Tensor7View() = default;
  Tensor7View(Numeric* data,
              const Range& l,
              const Range& v,
              const Range& s,
              const Range& b,
              const Range& p,
              const Range& r,
              const Range& c);
  Tensor7View(Numeric* data,
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
              const Range& nc);
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
  Tensor7() = default;
  Tensor7(Index l, Index v, Index s, Index b, Index p, Index r, Index c);
  Tensor7(Index l,
          Index v,
          Index s,
          Index b,
          Index p,
          Index r,
          Index c,
          Numeric fill);
  Tensor7(const ConstTensor7View& v);
  Tensor7(const Tensor7& v);
  Tensor7(Tensor7&& v) noexcept : Tensor7View(std::forward<Tensor7View>(v)) {
    v.mdata = nullptr;
  }

  /** Initialization from a tensor type. */
  explicit Tensor7(const matpack::tensor7_like_not_tensor7 auto &init)
      : Tensor7(matpack::library_size(init), matpack::vitrine_size(init),
                matpack::shelf_size(init), matpack::book_size(init),
                matpack::page_size(init), matpack::row_size(init),
                matpack::column_size(init)) {
    *this = init;
  }

  /** Set from a tensor type. */
  Tensor7 &operator=(const matpack::tensor7_like_not_tensor7 auto &init) {
    if (const auto s = matpack::shape<Index, 7>(init); shape().data not_eq s)
      resize(s[0], s[1], s[2], s[3], s[4], s[5], s[6]);

    auto [I, J, K, L, M, N, O] = shape().data;
    for (Index i = 0; i < I; i++)
      for (Index j = 0; j < J; j++)
        for (Index k = 0; k < K; k++)
          for (Index x = 0; x < L; x++)
            for (Index m = 0; m < M; m++)
              for (Index n = 0; n < N; n++)
                for (Index o = 0; o < O; o++)
                  operator()(i, j, k, x, m, n, o) = init(i, j, k, x, m, n, o);

    return *this;
  }

  // Assignment operators:
  Tensor7& operator=(const Tensor7& x);
  Tensor7& operator=(Tensor7&& x) ARTS_NOEXCEPT;
  Tensor7& operator=(Numeric x);

  // Resize function:
  void resize(Index l, Index v, Index s, Index b, Index p, Index r, Index c);

  // Swap function:
  friend void swap(Tensor7& t1, Tensor7& t2) noexcept;

  // Destructor:
  ~Tensor7() noexcept override;

  /*! Reduce a Tensor7 to a Vector and leave this in an empty state */
  template <std::size_t dim0>
      Vector reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim0 < 7, "Bad Dimension, Out-of-Bounds");

    Range r0(0,
             dim0 == 0   ? nlibraries()
             : dim0 == 1 ? nvitrines()
             : dim0 == 2 ? nshelves()
             : dim0 == 3 ? nbooks()
             : dim0 == 4 ? npages()
             : dim0 == 5 ? nrows()
                         : ncols());

    Vector out(mdata, r0);
    ARTS_ASSERT(size() not_eq out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor7 to a Matrix and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1>
      Matrix reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim1 < 7, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");

    const Range r1(0,
                   dim1 == 1   ? nvitrines()
                   : dim1 == 2 ? nshelves()
                   : dim1 == 3 ? nbooks()
                   : dim1 == 4 ? npages()
                   : dim1 == 5 ? nrows()
                               : ncols());
    const Range r0(0,
                   dim0 == 0   ? nlibraries()
                   : dim0 == 1 ? nvitrines()
                   : dim0 == 2 ? nshelves()
                   : dim0 == 3 ? nbooks()
                   : dim0 == 4 ? npages()
                               : nrows(),
                   r1.get_extent());

    Matrix out(mdata, r0, r1);
    ARTS_ASSERT(size() not_eq out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor7 to a Tensor3 and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1, std::size_t dim2>
      Tensor3 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim2 < 7, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");

    const Range r2(0,
                   dim2 == 2   ? nshelves()
                   : dim2 == 3 ? nbooks()
                   : dim2 == 4 ? npages()
                   : dim2 == 5 ? nrows()
                               : ncols());
    const Range r1(0,
                   dim1 == 1   ? nvitrines()
                   : dim1 == 2 ? nshelves()
                   : dim1 == 3 ? nbooks()
                   : dim1 == 4 ? npages()
                               : nrows(),
                   r2.get_extent());
    const Range r0(0,
                   dim0 == 0   ? nlibraries()
                   : dim0 == 1 ? nvitrines()
                   : dim0 == 2 ? nshelves()
                   : dim0 == 3 ? nbooks()
                               : npages(),
                   r1.get_extent() * r2.get_extent());

    Tensor3 out(mdata, r0, r1, r2);
    ARTS_ASSERT(size() not_eq out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor7 to a Tensor4 and leave this in an empty state */
  template <std::size_t dim0,
            std::size_t dim1,
            std::size_t dim2,
            std::size_t dim3>
      Tensor4 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim3 < 7, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");
    static_assert(dim2 < dim3, "Bad Dimensions, dim3 must be larger than dim2");

    const Range r3(0,
                   dim3 == 3   ? nbooks()
                   : dim3 == 4 ? npages()
                   : dim3 == 5 ? nrows()
                               : ncols());
    const Range r2(0,
                   dim2 == 2   ? nshelves()
                   : dim2 == 3 ? nbooks()
                   : dim2 == 4 ? npages()
                               : nrows(),
                   r3.get_extent());
    const Range r1(0,
                   dim1 == 1   ? nvitrines()
                   : dim1 == 2 ? nshelves()
                   : dim1 == 3 ? nbooks()
                               : npages(),
                   r2.get_extent() * r3.get_extent());
    const Range r0(0,
                   dim0 == 0   ? nlibraries()
                   : dim0 == 1 ? nvitrines()
                   : dim0 == 2 ? nshelves()
                               : nbooks(),
                   r1.get_extent() * r2.get_extent() * r3.get_extent());

    Tensor4 out(mdata, r0, r1, r2, r3);
    ARTS_ASSERT(size() not_eq out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor7 to a Tensor5 and leave this in an empty state */
  template <std::size_t dim0,
            std::size_t dim1,
            std::size_t dim2,
            std::size_t dim3,
            std::size_t dim4>
      Tensor5 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim4 < 7, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");
    static_assert(dim2 < dim3, "Bad Dimensions, dim3 must be larger than dim2");
    static_assert(dim3 < dim4, "Bad Dimensions, dim4 must be larger than dim3");

    const Range r4(0, dim4 == 4 ? npages() : dim4 == 5 ? nrows() : ncols());
    const Range r3(0,
                   dim3 == 3   ? nbooks()
                   : dim3 == 4 ? npages()
                               : nrows(),
                   r4.get_extent());
    const Range r2(0,
                   dim2 == 2   ? nshelves()
                   : dim2 == 3 ? nbooks()
                               : npages(),
                   r3.get_extent() * r4.get_extent());
    const Range r1(0,
                   dim1 == 1   ? nvitrines()
                   : dim1 == 2 ? nshelves()
                               : nbooks(),
                   r2.get_extent() * r3.get_extent() * r4.get_extent());
    const Range r0(
        0,
        dim0 == 0   ? nlibraries()
        : dim0 == 1 ? nvitrines()
                    : nshelves(),
        r1.get_extent() * r2.get_extent() * r3.get_extent() * r4.get_extent());

    Tensor5 out(mdata, r0, r1, r2, r3, r4);
    ARTS_ASSERT(size() not_eq out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor7 to a Tensor6 and leave this in an empty state */
  template <std::size_t dim0,
            std::size_t dim1,
            std::size_t dim2,
            std::size_t dim3,
            std::size_t dim4,
            std::size_t dim5>
      Tensor6 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim5 < 7, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");
    static_assert(dim2 < dim3, "Bad Dimensions, dim3 must be larger than dim2");
    static_assert(dim3 < dim4, "Bad Dimensions, dim4 must be larger than dim3");
    static_assert(dim4 < dim5, "Bad Dimensions, dim5 must be larger than dim4");

    const Range r5(0, dim5 == 5 ? nrows() : ncols());
    const Range r4(0, dim4 == 4 ? npages() : nrows(), r5.get_extent());
    const Range r3(
        0, dim3 == 3 ? nbooks() : npages(), r4.get_extent() * r5.get_extent());
    const Range r2(0,
                   dim2 == 2 ? nshelves() : nbooks(),
                   r3.get_extent() * r4.get_extent() * r5.get_extent());
    const Range r1(
        0,
        dim1 == 1 ? nvitrines() : nshelves(),
        r2.get_extent() * r3.get_extent() * r4.get_extent() * r5.get_extent());
    const Range r0(0,
                   dim0 == 0 ? nlibraries() : nvitrines(),
                   r1.get_extent() * r2.get_extent() * r3.get_extent() *
                       r4.get_extent() * r5.get_extent());

    Tensor6 out(mdata, r0, r1, r2, r3, r4, r5);
    ARTS_ASSERT(size() not_eq out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  template <class F>
  void transform_elementwise(F&& func) {
    std::transform(mdata, mdata + size(), mdata, func);
  }
};

// Function declarations:
// ----------------------

void copy(ConstIterator7D origin,
          const ConstIterator7D& end,
          Iterator7D target);

void copy(Numeric x, Iterator7D target, const Iterator7D& end);

void transform(Tensor7View y, double (&my_func)(double), ConstTensor7View x);

Numeric max(const ConstTensor7View& x);

Numeric min(const ConstTensor7View& x);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor7view_get_elem(Tensor7View& tv,
                                   Index l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   Index c);

#endif
////////////////////////////////

/** An array of Tensor7. */
using ArrayOfTensor7 = Array<Tensor7>;

using ArrayOfArrayOfTensor7 = Array<ArrayOfTensor7>;

#endif  // matpackVII_h
