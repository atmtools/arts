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

#include <utility>

#include "matpackV.h"

#define CHECK(x)       \
  ARTS_ASSERT(0 <= x); \
  ARTS_ASSERT(x < m##x##r.mextent)
#define OFFSET(x) m##x##r.mstart + x* m##x##r.mstride

/** The outermost iterator class for rank 6 tensors. This takes into
    account the defined strided. */
class Iterator6D {
 public:
  // Constructors:
  /** Default constructor. */
  Iterator6D() = default;

  /** Explicit constructor. */
  Iterator6D(const Tensor5View& x, Index stride)
      : msv(x), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  Iterator6D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy.
      FIXME: Is it really necessary to have such a complicated check
      here? It could be sufficient to just test
      msv.mdata!=other.msv.mdata. */
  bool operator!=(const Iterator6D& other) const {
    if (msv.mdata + msv.msr.mstart + msv.mbr.mstart + msv.mpr.mstart +
            msv.mrr.mstart + msv.mcr.mstart !=
        other.msv.mdata + other.msv.msr.mstart + other.msv.mbr.mstart +
            other.msv.mpr.mstart + other.msv.mrr.mstart + other.msv.mcr.mstart)
      return true;
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
  Index mstride{0};
};

/** Const version of Iterator6D. */
class ConstIterator6D {
 public:
  // Constructors:
  /** Default constructor. */
  ConstIterator6D() = default;

  /** Explicit constructor. */
  ConstIterator6D(ConstTensor5View x, Index stride)
      : msv(std::move(x)), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator6D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy.
      FIXME: Is it really necessaary to have such a complicated check
      here? It could be sufficient to just test
      msv.mdata!=other.msv.mdata. */
  bool operator!=(const ConstIterator6D& other) const {
    if (msv.mdata + msv.msr.mstart + msv.mbr.mstart + msv.mpr.mstart +
            msv.mrr.mstart + msv.mcr.mstart !=
        other.msv.mdata + other.msv.msr.mstart + other.msv.mbr.mstart +
            other.msv.mpr.mstart + other.msv.mrr.mstart + other.msv.mcr.mstart)
      return true;
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
  Index mstride{0};
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
  constexpr ConstTensor6View(const ConstTensor6View&) = default;
  constexpr ConstTensor6View(ConstTensor6View&&) = default;
  ConstTensor6View& operator=(const ConstTensor6View&) = default;
  ConstTensor6View& operator=(ConstTensor6View&&) = default;

  // Member functions:
  [[nodiscard]] Index nvitrines() const noexcept { return mvr.mextent; }
  [[nodiscard]] Index nshelves() const noexcept { return msr.mextent; }
  [[nodiscard]] Index nbooks() const noexcept { return mbr.mextent; }
  [[nodiscard]] Index npages() const noexcept { return mpr.mextent; }
  [[nodiscard]] Index nrows() const noexcept { return mrr.mextent; }
  [[nodiscard]] Index ncols() const noexcept { return mcr.mextent; }

  // Total size
  [[nodiscard]] Index size() const noexcept {
    return nvitrines() * nshelves() * nbooks() * npages() * nrows() * ncols();
  }
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /*! Returns the shape as an array (to allow templates to just look for shape on different matpack objects) */
  [[nodiscard]] Shape<6> shape() const {
    return {nvitrines(), nshelves(), nbooks(), npages(), nrows(), ncols()};
  }

  // Const index operators:

  // Result 6D (1 combination)
  // ------
  ConstTensor6View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 5D (6 combinations)
  // -----|
  ConstTensor5View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ----|-
  ConstTensor5View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ---|--
  ConstTensor5View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // --|---
  ConstTensor5View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // -|----
  ConstTensor5View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |-----
  ConstTensor5View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 4D (5+4+3+2+1 = 15 combinations)
  // ----||
  ConstTensor4View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // ---|-|
  ConstTensor4View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // --|--|
  ConstTensor4View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // -|---|
  ConstTensor4View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // |----|
  ConstTensor4View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ---||-
  ConstTensor4View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // --|-|-
  ConstTensor4View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // -|--|-
  ConstTensor4View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // |---|-
  ConstTensor4View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // --||--
  ConstTensor4View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -|-|--
  ConstTensor4View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |--|--
  ConstTensor4View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // -||---
  ConstTensor4View operator()(const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // |-|---
  ConstTensor4View operator()(Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  // ||----
  ConstTensor4View operator()(Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
  // ---|||
  ConstTensor3View operator()(const Range& v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              Index c) const;
  // --|-||
  ConstTensor3View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // -|--||
  ConstTensor3View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // |---||
  ConstTensor3View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              Index c) const;
  // --||-|
  ConstTensor3View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // -|-|-|
  ConstTensor3View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // |--|-|
  ConstTensor3View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              Index c) const;
  // -||--|
  ConstTensor3View operator()(const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // |-|--|
  ConstTensor3View operator()(Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // ||---|
  ConstTensor3View operator()(Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  // --|||-
  ConstTensor3View operator()(const Range& v,
                              const Range& s,
                              Index b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // -|-||-
  ConstTensor3View operator()(const Range& v,
                              Index s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // |--||-
  ConstTensor3View operator()(Index v,
                              const Range& s,
                              const Range& b,
                              Index p,
                              Index r,
                              const Range& c) const;
  // -||-|-
  ConstTensor3View operator()(const Range& v,
                              Index s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // |-|-|-
  ConstTensor3View operator()(Index v,
                              const Range& s,
                              Index b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // ||--|-
  ConstTensor3View operator()(Index v,
                              Index s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  // -|||--
  ConstTensor3View operator()(const Range& v,
                              Index s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |-||--
  ConstTensor3View operator()(Index v,
                              const Range& s,
                              Index b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // ||-|--
  ConstTensor3View operator()(Index v,
                              Index s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  // |||---
  ConstTensor3View operator()(Index v,
                              Index s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  // Result 2D (15 combinations)
  // IIII--
  ConstMatrixView operator()(
      Index v, Index s, Index b, Index p, const Range& r, const Range& c) const;
  // III-I-
  ConstMatrixView operator()(
      Index v, Index s, Index b, const Range& p, Index r, const Range& c) const;
  // II-II-
  ConstMatrixView operator()(
      Index v, Index s, const Range& b, Index p, Index r, const Range& c) const;
  // I-III-
  ConstMatrixView operator()(
      Index v, const Range& s, Index b, Index p, Index r, const Range& c) const;
  // -IIII-
  ConstMatrixView operator()(
      const Range& v, Index s, Index b, Index p, Index r, const Range& c) const;
  // III--I
  ConstMatrixView operator()(
      Index v, Index s, Index b, const Range& p, const Range& r, Index c) const;
  // II-I-I
  ConstMatrixView operator()(
      Index v, Index s, const Range& b, Index p, const Range& r, Index c) const;
  // I-II-I
  ConstMatrixView operator()(
      Index v, const Range& s, Index b, Index p, const Range& r, Index c) const;
  // -III-I
  ConstMatrixView operator()(
      const Range& v, Index s, Index b, Index p, const Range& r, Index c) const;
  // II--II
  ConstMatrixView operator()(
      Index v, Index s, const Range& b, const Range& p, Index r, Index c) const;
  // I-I-II
  ConstMatrixView operator()(
      Index v, const Range& s, Index b, const Range& p, Index r, Index c) const;
  // -II-II
  ConstMatrixView operator()(
      const Range& v, Index s, Index b, const Range& p, Index r, Index c) const;
  // I--III
  ConstMatrixView operator()(
      Index v, const Range& s, const Range& b, Index p, Index r, Index c) const;
  // -I-III
  ConstMatrixView operator()(
      const Range& v, Index s, const Range& b, Index p, Index r, Index c) const;
  // --IIII
  ConstMatrixView operator()(
      const Range& v, const Range& s, Index b, Index p, Index r, Index c) const;

  // Result 1D (6 combinations)
  // IIIII-
  ConstVectorView operator()(
      Index v, Index s, Index b, Index p, Index r, const Range& c) const;
  // IIII-I
  ConstVectorView operator()(
      Index v, Index s, Index b, Index p, const Range& r, Index c) const;
  // III-II
  ConstVectorView operator()(
      Index v, Index s, Index b, const Range& p, Index r, Index c) const;
  // II-III
  ConstVectorView operator()(
      Index v, Index s, const Range& b, Index p, Index r, Index c) const;
  // I-IIII
  ConstVectorView operator()(
      Index v, const Range& s, Index b, Index p, Index r, Index c) const;
  // -IIIII
  ConstVectorView operator()(
      const Range& v, Index s, Index b, Index p, Index r, Index c) const;

  // Result scalar (1 combination)
  // IIIIII
  Numeric operator()(
      Index v, Index s, Index b, Index p, Index r, Index c) const {
    CHECK(v);
    CHECK(s);
    CHECK(b);
    CHECK(p);
    CHECK(r);
    CHECK(c);
    return get(v, s, b, p, r, c);
  }

  /** Get element implementation without assertions. */
  [[nodiscard]] Numeric get(
      Index v, Index s, Index b, Index p, Index r, Index c) const {
    return *(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) +
             OFFSET(c));
  }

  // Functions returning iterators:
  [[nodiscard]] ConstIterator6D begin() const;
  [[nodiscard]] ConstIterator6D end() const;

  // Destructor:
  virtual ~ConstTensor6View() = default;

  // Friends:
  friend class ConstIterator7D;
  friend class Tensor6View;
  friend class ConstTensor7View;

  friend std::ostream& operator<<(std::ostream& os, const ConstTensor6View& v);

  // Special constructor to make a Tensor6 view of a Tensor5.
  ConstTensor6View(const ConstTensor5View& a);

 protected:
  // Constructors:
  ConstTensor6View() = default;
  ConstTensor6View(Numeric* data,
                   const Range& v,
                   const Range& s,
                   const Range& b,
                   const Range& p,
                   const Range& r,
                   const Range& c);
  ConstTensor6View(Numeric* data,
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
                   const Range& nc);

  // Data members:
  // -------------
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

/** The Tensor6View class

    This contains the main implementation of a Tensor6. It defines
    the concepts of Tensor6View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor6View
    from a subrange of a Tensor6View. 

    The class Tensor6 is just a special case of a Tensor6View
    which also allocates storage. */
class Tensor6View : public ConstTensor6View {
 public:
  // Make const methods visible from base class
  using ConstTensor6View::begin;
  using ConstTensor6View::end;
  using ConstTensor6View::operator();
  using ConstTensor6View::get;

  constexpr Tensor6View(const Tensor6View&) = default;

  // Non-const index operators:

  // Result 6D (1 combination)
  // ------
  Tensor6View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 5D (6 combinations)
  // -----|
  Tensor5View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ----|-
  Tensor5View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ---|--
  Tensor5View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // --|---
  Tensor5View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // -|----
  Tensor5View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |-----
  Tensor5View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 4D (5+4+3+2+1 = 15 combinations)
  // ----||
  Tensor4View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // ---|-|
  Tensor4View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // --|--|
  Tensor4View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // -|---|
  Tensor4View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // |----|
  Tensor4View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ---||-
  Tensor4View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // --|-|-
  Tensor4View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // -|--|-
  Tensor4View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // |---|-
  Tensor4View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // --||--
  Tensor4View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -|-|--
  Tensor4View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |--|--
  Tensor4View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // -||---
  Tensor4View operator()(const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // |-|---
  Tensor4View operator()(Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);
  // ||----
  Tensor4View operator()(Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
  // ---|||
  Tensor3View operator()(const Range& v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         Index c);
  // --|-||
  Tensor3View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         Index c);
  // -|--||
  Tensor3View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // |---||
  Tensor3View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         Index c);
  // --||-|
  Tensor3View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         Index c);
  // -|-|-|
  Tensor3View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // |--|-|
  Tensor3View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         Index c);
  // -||--|
  Tensor3View operator()(const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // |-|--|
  Tensor3View operator()(Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // ||---|
  Tensor3View operator()(Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         Index c);
  // --|||-
  Tensor3View operator()(const Range& v,
                         const Range& s,
                         Index b,
                         Index p,
                         Index r,
                         const Range& c);
  // -|-||-
  Tensor3View operator()(const Range& v,
                         Index s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // |--||-
  Tensor3View operator()(Index v,
                         const Range& s,
                         const Range& b,
                         Index p,
                         Index r,
                         const Range& c);
  // -||-|-
  Tensor3View operator()(const Range& v,
                         Index s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // |-|-|-
  Tensor3View operator()(Index v,
                         const Range& s,
                         Index b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // ||--|-
  Tensor3View operator()(Index v,
                         Index s,
                         const Range& b,
                         const Range& p,
                         Index r,
                         const Range& c);
  // -|||--
  Tensor3View operator()(const Range& v,
                         Index s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |-||--
  Tensor3View operator()(Index v,
                         const Range& s,
                         Index b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // ||-|--
  Tensor3View operator()(Index v,
                         Index s,
                         const Range& b,
                         Index p,
                         const Range& r,
                         const Range& c);
  // |||---
  Tensor3View operator()(Index v,
                         Index s,
                         Index b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  // Result 2D (15 combinations)
  // IIII--
  MatrixView operator()(
      Index v, Index s, Index b, Index p, const Range& r, const Range& c);
  // III-I-
  MatrixView operator()(
      Index v, Index s, Index b, const Range& p, Index r, const Range& c);
  // II-II-
  MatrixView operator()(
      Index v, Index s, const Range& b, Index p, Index r, const Range& c);
  // I-III-
  MatrixView operator()(
      Index v, const Range& s, Index b, Index p, Index r, const Range& c);
  // -IIII-
  MatrixView operator()(
      const Range& v, Index s, Index b, Index p, Index r, const Range& c);
  // III--I
  MatrixView operator()(
      Index v, Index s, Index b, const Range& p, const Range& r, Index c);
  // II-I-I
  MatrixView operator()(
      Index v, Index s, const Range& b, Index p, const Range& r, Index c);
  // I-II-I
  MatrixView operator()(
      Index v, const Range& s, Index b, Index p, const Range& r, Index c);
  // -III-I
  MatrixView operator()(
      const Range& v, Index s, Index b, Index p, const Range& r, Index c);
  // II--II
  MatrixView operator()(
      Index v, Index s, const Range& b, const Range& p, Index r, Index c);
  // I-I-II
  MatrixView operator()(
      Index v, const Range& s, Index b, const Range& p, Index r, Index c);
  // -II-II
  MatrixView operator()(
      const Range& v, Index s, Index b, const Range& p, Index r, Index c);
  // I--III
  MatrixView operator()(
      Index v, const Range& s, const Range& b, Index p, Index r, Index c);
  // -I-III
  MatrixView operator()(
      const Range& v, Index s, const Range& b, Index p, Index r, Index c);
  // --IIII
  MatrixView operator()(
      const Range& v, const Range& s, Index b, Index p, Index r, Index c);

  // Result 1D (6 combinations)
  // IIIII-
  VectorView operator()(
      Index v, Index s, Index b, Index p, Index r, const Range& c);
  // IIII-I
  VectorView operator()(
      Index v, Index s, Index b, Index p, const Range& r, Index c);
  // III-II
  VectorView operator()(
      Index v, Index s, Index b, const Range& p, Index r, Index c);
  // II-III
  VectorView operator()(
      Index v, Index s, const Range& b, Index p, Index r, Index c);
  // I-IIII
  VectorView operator()(
      Index v, const Range& s, Index b, Index p, Index r, Index c);
  // -IIIII
  VectorView operator()(
      const Range& v, Index s, Index b, Index p, Index r, Index c);

#define GETFUN(v, s, b, p, r, c)                                        \
  *(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + \
    OFFSET(c))
  // Result scalar (1 combination)
  // IIIIII
  Numeric& operator()(Index v, Index s, Index b, Index p, Index r, Index c) {
    CHECK(v);
    CHECK(s);
    CHECK(b);
    CHECK(p);
    CHECK(r);
    CHECK(c);
    return GETFUN(v, s, b, p, r, c);
  }

  /** Get element implementation without assertions. */
  Numeric& get(Index v, Index s, Index b, Index p, Index r, Index c) {
    return GETFUN(v, s, b, p, r, c);
  }
#undef GETFUN

  // Conversion to a plain C-array
  [[nodiscard]] const Numeric* get_c_array() const ARTS_NOEXCEPT;
  Numeric* get_c_array() ARTS_NOEXCEPT;

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
  ~Tensor6View() override = default;

  // Friends:
  friend class Iterator7D;
  friend class Tensor7View;

  // Special constructor to make a Tensor6 view of a Tensor5.
  Tensor6View(const Tensor5View& a);

 protected:
  // Constructors:
  Tensor6View() = default;
  Tensor6View(Numeric* data,
              const Range& v,
              const Range& s,
              const Range& b,
              const Range& p,
              const Range& r,
              const Range& c);
  Tensor6View(Numeric* data,
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
              const Range& nc);
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
  Tensor6() = default;
  Tensor6(Index v, Index s, Index b, Index p, Index r, Index c);
  Tensor6(Index v, Index s, Index b, Index p, Index r, Index c, Numeric fill);
  Tensor6(const ConstTensor6View& v);
  Tensor6(const Tensor6& v);
  Tensor6(Tensor6&& v) noexcept : Tensor6View(std::forward<Tensor6View>(v)) {
    v.mdata = nullptr;
  }

  /*! Construct from known data
   * 
   * Note that this will call delete on the pointer if it is still valid
   * at the end of the lifetime of this variable
   * 
   * @param[in] d - A pointer to some raw data
   * @param[in] r0 - The Range along the first dimension
   * @param[in] r1 - The Range along the second dimension
   * @param[in] r2 - The Range along the third dimension
   * @param[in] r3 - The Range along the fourth dimension
   * @param[in] r4 - The Range along the fifth dimension
   * @param[in] r5 - The Range along the sixth dimension
   */
  Tensor6(Numeric* d,
          const Range& r0,
          const Range& r1,
          const Range& r2,
          const Range& r3,
          const Range& r4,
          const Range& r5) ARTS_NOEXCEPT
      : Tensor6View(d, r0, r1, r2, r3, r4, r5) {
    ARTS_ASSERT(not(r0.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r1.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r2.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r3.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r4.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r5.get_extent() < 0), "Must have size");
  }

  // Assignment operators:
  Tensor6& operator=(const Tensor6& x);
  Tensor6& operator=(Tensor6&& x) ARTS_NOEXCEPT;
  Tensor6& operator=(Numeric x);

  // Resize function:
  void resize(Index v, Index s, Index b, Index p, Index r, Index c);

  // Swap function:
  friend void swap(Tensor6& t1, Tensor6& t2) noexcept;

  // Destructor:
  ~Tensor6() noexcept override;

  /*! Reduce a Tensor6 to a Vector and leave this in an empty state */
  template <std::size_t dim0>
      Vector reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim0 < 6, "Bad Dimension, Out-of-Bounds");

    Range r0(0,
             dim0 == 0   ? nvitrines()
             : dim0 == 1 ? nshelves()
             : dim0 == 2 ? nbooks()
             : dim0 == 3 ? npages()
             : dim0 == 4 ? nrows()
                         : ncols());

    Vector out(mdata, r0);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor6 to a Matrix and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1>
      Matrix reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim1 < 6, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");

    const Range r1(0,
                   dim1 == 1   ? nshelves()
                   : dim1 == 2 ? nbooks()
                   : dim1 == 3 ? npages()
                   : dim1 == 4 ? nrows()
                               : ncols());
    const Range r0(0,
                   dim0 == 0   ? nvitrines()
                   : dim0 == 1 ? nshelves()
                   : dim0 == 2 ? nbooks()
                   : dim0 == 3 ? npages()
                               : nrows(),
                   r1.get_extent());

    Matrix out(mdata, r0, r1);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor6 to a Tensor3 and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1, std::size_t dim2>
      Tensor3 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim2 < 6, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");

    const Range r2(0,
                   dim2 == 2   ? nbooks()
                   : dim2 == 3 ? npages()
                   : dim2 == 4 ? nrows()
                               : ncols());
    const Range r1(0,
                   dim1 == 1   ? nshelves()
                   : dim1 == 2 ? nbooks()
                   : dim1 == 3 ? npages()
                               : nrows(),
                   r2.get_extent());
    const Range r0(0,
                   dim0 == 0   ? nvitrines()
                   : dim0 == 1 ? nshelves()
                   : dim0 == 2 ? nbooks()
                               : npages(),
                   r1.get_extent() * r2.get_extent());

    Tensor3 out(mdata, r0, r1, r2);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor6 to a Tensor4 and leave this in an empty state */
  template <std::size_t dim0,
            std::size_t dim1,
            std::size_t dim2,
            std::size_t dim3>
      Tensor4 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim3 < 6, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");
    static_assert(dim2 < dim3, "Bad Dimensions, dim3 must be larger than dim2");

    const Range r3(0, dim3 == 3 ? npages() : dim3 == 4 ? nrows() : ncols());
    const Range r2(0,
                   dim2 == 2   ? nbooks()
                   : dim2 == 3 ? npages()
                               : nrows(),
                   r3.get_extent());
    const Range r1(0,
                   dim1 == 1   ? nshelves()
                   : dim1 == 2 ? nbooks()
                               : npages(),
                   r2.get_extent() * r3.get_extent());
    const Range r0(0,
                   dim0 == 0   ? nvitrines()
                   : dim0 == 1 ? nshelves()
                               : nbooks(),
                   r1.get_extent() * r2.get_extent() * r3.get_extent());

    Tensor4 out(mdata, r0, r1, r2, r3);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor6 to a Tensor5 and leave this in an empty state */
  template <std::size_t dim0,
            std::size_t dim1,
            std::size_t dim2,
            std::size_t dim3,
            std::size_t dim4>
      Tensor5 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim4 < 6, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");
    static_assert(dim2 < dim3, "Bad Dimensions, dim3 must be larger than dim2");
    static_assert(dim3 < dim4, "Bad Dimensions, dim4 must be larger than dim3");

    const Range r4(0, dim4 == 4 ? nrows() : ncols());
    const Range r3(0, dim3 == 3 ? npages() : nrows(), r4.get_extent());
    const Range r2(
        0, dim2 == 2 ? nbooks() : npages(), r3.get_extent() * r4.get_extent());
    const Range r1(0,
                   dim1 == 1 ? nshelves() : nbooks(),
                   r2.get_extent() * r3.get_extent() * r4.get_extent());
    const Range r0(
        0,
        dim0 == 0 ? nvitrines() : nshelves(),
        r1.get_extent() * r2.get_extent() * r3.get_extent() * r4.get_extent());

    Tensor5 out(mdata, r0, r1, r2, r3, r4);
    ARTS_ASSERT(size() == out.size(),
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

void copy(ConstIterator6D origin,
          const ConstIterator6D& end,
          Iterator6D target);

void copy(Numeric x, Iterator6D target, const Iterator6D& end);

void transform(Tensor6View y, double (&my_func)(double), ConstTensor6View x);

Numeric max(const ConstTensor6View& x);

Numeric min(const ConstTensor6View& x);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor6view_get_elem(
    Tensor6View& tv, Index v, Index s, Index b, Index p, Index r, Index c);

#endif
////////////////////////////////

/** An array of Tensor6. */
using ArrayOfTensor6 = Array<Tensor6>;

using ArrayOfArrayOfTensor6 = Array<ArrayOfTensor6>;

#endif  // matpackVI_h
