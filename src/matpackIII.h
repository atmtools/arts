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
  Implementation of Tensors of Rank 3.

  The three dimensions are called: page, row, column.
  
  \author Stefan Buehler
  \date   2001-11-22
 */

#ifndef matpackIII_h
#define matpackIII_h

#include <utility>

#include "matpackI.h"

/** The outermost iterator class for rank 3 tensors. This takes into
    account the defined strided. */
class Iterator3D {
 public:
  // Constructors:
  /** Default constructor. */
  Iterator3D() = default;

  /** Explicit constructor. */
  Iterator3D(const MatrixView& x, Index stride)
      : msv(x), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  Iterator3D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const Iterator3D& other) const {
    if (msv.mdata + msv.mrr.mstart + msv.mcr.mstart !=
        other.msv.mdata + other.msv.mrr.mstart + other.msv.mcr.mstart)
      return true;
    return false;
  }

  /** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
  MatrixView* operator->() { return &msv; }

  /** Dereferencing. */
  MatrixView& operator*() { return msv; }

 private:
  /** Current position. */
  MatrixView msv;
  /** Stride. */
  Index mstride{0};
};

/** Const version of Iterator3D. */
class ConstIterator3D {
 public:
  // Constructors:
  /** Default constructor. */
  ConstIterator3D() = default;

  /** Explicit constructor. */
  ConstIterator3D(ConstMatrixView x, Index stride)
      : msv(std::move(x)), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator3D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const ConstIterator3D& other) const {
    if (msv.mdata + msv.mrr.mstart + msv.mcr.mstart !=
        other.msv.mdata + other.msv.mrr.mstart + other.msv.mcr.mstart)
      return true;
    return false;
  }

  /** The -> operator is needed, so that we can write i->begin() to get
    the 1D iterators. */
  const ConstMatrixView* operator->() const { return &msv; }

  /** Dereferencing. */
  const ConstMatrixView& operator*() const { return msv; }

 private:
  /** Current position. */
  ConstMatrixView msv;
  /** Stride. */
  Index mstride{0};
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
  constexpr ConstTensor3View(const ConstTensor3View&) = default;
  constexpr ConstTensor3View(ConstTensor3View&&) = default;
  ConstTensor3View& operator=(const ConstTensor3View&) = default;
  ConstTensor3View& operator=(ConstTensor3View&&) = default;

  // Member functions:

  /** Returns the number of pages. */
  [[nodiscard]] Index npages() const { return mpr.mextent; }

  /** Returns the number of rows. */
  [[nodiscard]] Index nrows() const { return mrr.mextent; }

  /** Returns the number of columns. */
  [[nodiscard]] Index ncols() const { return mcr.mextent; }

  // Total size
  [[nodiscard]] Index size() const noexcept {
    return npages() * nrows() * ncols();
  }
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /*! Returns the shape as an array (to allow templates to just look for shape on different matpack objects) */
  [[nodiscard]] Shape<3> shape() const { return {npages(), nrows(), ncols()}; }

  // Const index operators:
  ConstTensor3View operator()(const Range& p,
                              const Range& r,
                              const Range& c) const;

  ConstMatrixView operator()(const Range& p, const Range& r, Index c) const;
  ConstMatrixView operator()(const Range& p, Index r, const Range& c) const;
  ConstMatrixView operator()(Index p, const Range& r, const Range& c) const;

  ConstVectorView operator()(Index p, Index r, const Range& c) const;
  ConstVectorView operator()(Index p, const Range& r, Index c) const;
  ConstVectorView operator()(const Range& p, Index r, Index c) const;

  /** Plain const index operator. */
  Numeric operator()(Index p,
                     Index r,
                     Index c) const {  // Check if indices are valid:
    ARTS_ASSERT(0 <= p);
    ARTS_ASSERT(0 <= r);
    ARTS_ASSERT(0 <= c);
    ARTS_ASSERT(p < mpr.mextent);
    ARTS_ASSERT(r < mrr.mextent);
    ARTS_ASSERT(c < mcr.mextent);

    return get(p, r, c);
  }

  /** Get element implementation without assertions. */
  [[nodiscard]] Numeric get(Index p, Index r, Index c) const {
    return *(mdata + mpr.mstart + p * mpr.mstride + mrr.mstart +
             r * mrr.mstride + mcr.mstart + c * mcr.mstride);
  }

  // Functions returning iterators:
  [[nodiscard]] ConstIterator3D begin() const;
  [[nodiscard]] ConstIterator3D end() const;

  //! Destructor
  virtual ~ConstTensor3View() = default;

  // Friends:
  friend class Tensor3View;
  friend class ConstIterator4D;
  friend class ConstTensor4View;
  friend class ConstTensor5View;
  friend class ConstTensor6View;
  friend class ConstTensor7View;

  friend std::ostream& operator<<(std::ostream& os, const ConstTensor3View& v);

  // Special constructor to make a Tensor3 view of a matrix.
  ConstTensor3View(const ConstMatrixView& a);

 protected:
  // Constructors:
  ConstTensor3View() = default;
  ConstTensor3View(Numeric* data,
                   const Range& p,
                   const Range& r,
                   const Range& c);
  ConstTensor3View(Numeric* data,
                   const Range& pp,
                   const Range& pr,
                   const Range& pc,
                   const Range& np,
                   const Range& nr,
                   const Range& nc);

  // Data members:
  // -------------
  /** The page range of mdata that is actually used. */
  Range mpr{0, 0, 1};
  /** The row range of mdata that is actually used. */
  Range mrr{0, 0, 1};
  /** The column range of mdata that is actually used. */
  Range mcr{0, 0, 1};
  /** Pointer to the plain C array that holds the data */
  Numeric* mdata{nullptr};
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
  // Make const methods visible from base class
  using ConstTensor3View::begin;
  using ConstTensor3View::end;
  using ConstTensor3View::operator();
  using ConstTensor3View::get;

  constexpr Tensor3View(const Tensor3View&) = default;

  // Non-const index operators:

  Tensor3View operator()(const Range& p, const Range& r, const Range& c);

  MatrixView operator()(const Range& p, const Range& r, Index c);
  MatrixView operator()(const Range& p, Index r, const Range& c);
  MatrixView operator()(Index p, const Range& r, const Range& c);

  VectorView operator()(Index p, Index r, const Range& c);
  VectorView operator()(Index p, const Range& r, Index c);
  VectorView operator()(const Range& p, Index r, Index c);

  /** Plain non-const index operator. */
  Numeric& operator()(Index p, Index r, Index c) {
    // Check if indices are valid:
    ARTS_ASSERT(0 <= p);
    ARTS_ASSERT(0 <= r);
    ARTS_ASSERT(0 <= c);
    ARTS_ASSERT(p < mpr.mextent);
    ARTS_ASSERT(r < mrr.mextent);
    ARTS_ASSERT(c < mcr.mextent);

    return get(p, r, c);
  }

  /** Get element implementation without assertions. */
  Numeric& get(Index p, Index r, Index c) {
    return *(mdata + mpr.mstart + p * mpr.mstride + mrr.mstart +
             r * mrr.mstride + mcr.mstart + c * mcr.mstride);
  }

  // Conversion to a plain C-array
  [[nodiscard]] const Numeric* get_c_array() const ARTS_NOEXCEPT;
  Numeric* get_c_array() ARTS_NOEXCEPT;

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

  //! Destructor
  ~Tensor3View() override = default;

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
  Tensor3View() = default;
  Tensor3View(Numeric* data, const Range& p, const Range& r, const Range& c);
  Tensor3View(Numeric* data,
              const Range& pp,
              const Range& pr,
              const Range& pc,
              const Range& np,
              const Range& nr,
              const Range& nc);
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
  Tensor3() = default;
  Tensor3(Index p, Index r, Index c);
  Tensor3(Index p, Index r, Index c, Numeric fill);
  Tensor3(const ConstTensor3View& v);
  Tensor3(const Tensor3& v);
  Tensor3(Tensor3&& v) noexcept : Tensor3View(std::forward<Tensor3View>(v)) {
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
   */
  Tensor3(Numeric* d, const Range& r0, const Range& r1, const Range& r2)
      ARTS_NOEXCEPT : Tensor3View(d, r0, r1, r2) {
    ARTS_ASSERT(not(r0.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r1.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r2.get_extent() < 0), "Must have size");
  }

  // Assignment operators:
  Tensor3& operator=(const Tensor3& x);
  Tensor3& operator=(Tensor3&& x) noexcept;
  Tensor3& operator=(Numeric x);

  // Resize function:
  void resize(Index p, Index r, Index c);

  // Swap function:
  friend void swap(Tensor3& t1, Tensor3& t2) noexcept;

  // Destructor:
  ~Tensor3() noexcept override;

  /*! Reduce a Tensor3 to a Vector and leave this in an empty state */
  template <std::size_t dim0>
      Vector reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim0 < 3, "Bad Dimension, Out-of-Bounds");

    Range r0(0, dim0 == 0 ? npages() : dim0 == 1 ? nrows() : ncols());

    Vector out(mdata, r0);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor3 to a Matrix and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1>
      Matrix reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim1 < 3, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");

    const Range r1(0, dim1 == 1 ? nrows() : ncols());
    const Range r0(0, dim0 == 0 ? npages() : nrows(), r1.get_extent());

    Matrix out(mdata, r0, r1);
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

void copy(ConstIterator3D origin,
          const ConstIterator3D& end,
          Iterator3D target);

void copy(Numeric x, Iterator3D target, const Iterator3D& end);

void transform(Tensor3View y, double (&my_func)(double), ConstTensor3View x);

Numeric max(const ConstTensor3View& x);

Numeric min(const ConstTensor3View& x);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor3view_get_elem(Tensor3View& tv, Index p, Index r, Index c);

#endif
////////////////////////////////

void mult(Tensor3View A, const ConstVectorView B, const ConstMatrixView C);

/** An array of Tensor3. */
using ArrayOfTensor3 = Array<Tensor3>;

using ArrayOfArrayOfTensor3 = Array<ArrayOfTensor3>;

#endif  // matpackIII_h
