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

#include <utility>

#include "matpackIV.h"
#include "matpack_concepts.h"

/** The outermost iterator class for rank 5 tensors. This takes into
    account the defined strided. */
class Iterator5D {
 public:
  // Constructors:
  // Functions for Iterator5D
  // ------------------------

  /** Default constructor. */
  Iterator5D() = default;

  /** Explicit constructor. */
  Iterator5D(const Tensor4View& x, Index stride)
      : msv(x), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  Iterator5D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const Iterator5D& other) const {
    if (msv.mdata + msv.mbr.mstart + msv.mpr.mstart + msv.mrr.mstart +
            msv.mcr.mstart !=
        other.msv.mdata + other.msv.mbr.mstart + other.msv.mpr.mstart +
            other.msv.mrr.mstart + other.msv.mcr.mstart)
      return true;
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
  Index mstride{0};
};

/** Const version of Iterator5D. */
class ConstIterator5D {
 public:
  // Constructors:
  /** Default constructor. */
  ConstIterator5D() = default;

  /** Explicit constructor. */
  ConstIterator5D(ConstTensor4View x, Index stride)
      : msv(std::move(x)), mstride(stride) { /* Nothing to do here. */
  }

  // Operators:
  /** Prefix increment operator. */
  ConstIterator5D& operator++() {
    msv.mdata += mstride;
    return *this;
  }

  /** Not equal operator, needed for algorithms like copy. */
  bool operator!=(const ConstIterator5D& other) const {
    if (msv.mdata + msv.mbr.mstart + msv.mpr.mstart + msv.mrr.mstart +
            msv.mcr.mstart !=
        other.msv.mdata + other.msv.mbr.mstart + other.msv.mpr.mstart +
            other.msv.mrr.mstart + other.msv.mcr.mstart)
      return true;
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
  Index mstride{0};
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
  static constexpr bool matpack_type{true};
  
  constexpr ConstTensor5View(const ConstTensor5View&) = default;
  constexpr ConstTensor5View(ConstTensor5View&&) = default;
  ConstTensor5View& operator=(const ConstTensor5View&) = default;
  ConstTensor5View& operator=(ConstTensor5View&&) = default;

  // Member functions:
  [[nodiscard]] Index nshelves() const noexcept { return msr.mextent; }
  [[nodiscard]] Index nbooks() const noexcept { return mbr.mextent; }
  [[nodiscard]] Index npages() const noexcept { return mpr.mextent; }
  [[nodiscard]] Index nrows() const noexcept { return mrr.mextent; }
  [[nodiscard]] Index ncols() const noexcept { return mcr.mextent; }

  // Total size
  [[nodiscard]] Index size() const noexcept {
    return nshelves() * nbooks() * npages() * nrows() * ncols();
  }
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /*! Returns the shape as an array (to allow templates to just look for shape on different matpack objects) */
  [[nodiscard]] Shape<5> shape() const {
    return {nshelves(), nbooks(), npages(), nrows(), ncols()};
  }

  // Const index operators:
  ConstTensor5View operator()(const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  ConstTensor4View operator()(const Range& s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              Index c) const;
  ConstTensor4View operator()(const Range& s,
                              const Range& b,
                              const Range& p,
                              Index r,
                              const Range& c) const;
  ConstTensor4View operator()(const Range& s,
                              const Range& b,
                              Index p,
                              const Range& r,
                              const Range& c) const;
  ConstTensor4View operator()(const Range& s,
                              Index b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;
  ConstTensor4View operator()(Index s,
                              const Range& b,
                              const Range& p,
                              const Range& r,
                              const Range& c) const;

  ConstTensor3View operator()(
      const Range& s, const Range& b, const Range& p, Index r, Index c) const;
  ConstTensor3View operator()(
      const Range& s, const Range& b, Index p, const Range& r, Index c) const;
  ConstTensor3View operator()(
      const Range& s, const Range& b, Index p, Index r, const Range& c) const;
  ConstTensor3View operator()(
      const Range& s, Index b, const Range& p, Index r, const Range& c) const;
  ConstTensor3View operator()(
      const Range& s, Index b, const Range& p, const Range& r, Index c) const;
  ConstTensor3View operator()(
      const Range& s, Index b, Index p, const Range& r, const Range& c) const;
  ConstTensor3View operator()(
      Index s, const Range& b, Index p, const Range& r, const Range& c) const;
  ConstTensor3View operator()(
      Index s, const Range& b, const Range& p, Index r, const Range& c) const;
  ConstTensor3View operator()(
      Index s, const Range& b, const Range& p, const Range& r, Index c) const;
  ConstTensor3View operator()(
      Index s, Index b, const Range& p, const Range& r, const Range& c) const;

  ConstMatrixView operator()(
      const Range& s, const Range& b, Index p, Index r, Index c) const;
  ConstMatrixView operator()(
      const Range& s, Index b, const Range& p, Index r, Index c) const;
  ConstMatrixView operator()(
      const Range& s, Index b, Index p, const Range& r, Index c) const;
  ConstMatrixView operator()(
      const Range& s, Index b, Index p, Index r, const Range& c) const;
  ConstMatrixView operator()(
      Index s, const Range& b, Index p, Index r, const Range& c) const;
  ConstMatrixView operator()(
      Index s, const Range& b, Index p, const Range& r, Index c) const;
  ConstMatrixView operator()(
      Index s, const Range& b, const Range& p, Index r, Index c) const;
  ConstMatrixView operator()(
      Index s, Index b, const Range& p, const Range& r, Index c) const;
  ConstMatrixView operator()(
      Index s, Index b, const Range& p, Index r, const Range& c) const;
  ConstMatrixView operator()(
      Index s, Index b, Index p, const Range& r, const Range& c) const;

  ConstVectorView operator()(
      const Range& s, Index b, Index p, Index r, Index c) const;
  ConstVectorView operator()(
      Index s, const Range& b, Index p, Index r, Index c) const;
  ConstVectorView operator()(
      Index s, Index b, const Range& p, Index r, Index c) const;
  ConstVectorView operator()(
      Index s, Index b, Index p, const Range& r, Index c) const;
  ConstVectorView operator()(
      Index s, Index b, Index p, Index r, const Range& c) const;

  /** Plain const index operator. */
  Numeric operator()(Index s,
                     Index b,
                     Index p,
                     Index r,
                     Index c) const {  // Check if indices are valid:
    ARTS_ASSERT(0 <= s);
    ARTS_ASSERT(0 <= b);
    ARTS_ASSERT(0 <= p);
    ARTS_ASSERT(0 <= r);
    ARTS_ASSERT(0 <= c);
    ARTS_ASSERT(s < msr.mextent);
    ARTS_ASSERT(b < mbr.mextent);
    ARTS_ASSERT(p < mpr.mextent);
    ARTS_ASSERT(r < mrr.mextent);
    ARTS_ASSERT(c < mcr.mextent);

    return get(s, b, p, r, c);
  }

  /** Get element implementation without assertions. */
  [[nodiscard]] Numeric get(Index s, Index b, Index p, Index r, Index c) const {
    return *(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
             b * mbr.mstride + mpr.mstart + p * mpr.mstride + mrr.mstart +
             r * mrr.mstride + mcr.mstart + c * mcr.mstride);
  }

  // Functions returning iterators:
  [[nodiscard]] ConstIterator5D begin() const;
  [[nodiscard]] ConstIterator5D end() const;

  //! Destructor
  virtual ~ConstTensor5View() = default;

  // Friends:
  friend class Tensor5View;
  friend class ConstIterator6D;
  friend class ConstTensor6View;
  friend class ConstTensor7View;

  friend std::ostream& operator<<(std::ostream& os, const ConstTensor5View& v);

  // Special constructor to make a Tensor5 view of a Tensor4.
  ConstTensor5View(const ConstTensor4View& a);

 protected:
  // Constructors:
  ConstTensor5View() = default;
  ConstTensor5View(Numeric* data,
                   const Range& s,
                   const Range& b,
                   const Range& p,
                   const Range& r,
                   const Range& c);
  ConstTensor5View(Numeric* data,
                   const Range& ps,
                   const Range& pb,
                   const Range& pp,
                   const Range& pr,
                   const Range& pc,
                   const Range& ns,
                   const Range& nb,
                   const Range& np,
                   const Range& nr,
                   const Range& nc);

  // Data members:
  // -------------
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

/** The Tensor5View class

    This contains the main implementation of a Tensor5. It defines
    the concepts of Tensor5View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor5View
    from a subrange of a Tensor5View.

    The class Tensor5 is just a special case of a Tensor5View
    which also allocates storage. */
class Tensor5View : public ConstTensor5View {
 public:
  // Make const methods visible from base class
  using ConstTensor5View::begin;
  using ConstTensor5View::end;
  using ConstTensor5View::operator();
  using ConstTensor5View::get;

  constexpr Tensor5View(const Tensor5View&) = default;

  // Non-const index operators:

  Tensor5View operator()(const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c);

  Tensor4View operator()(
      const Range& s, const Range& b, const Range& p, const Range& r, Index c);
  Tensor4View operator()(
      const Range& s, const Range& b, const Range& p, Index r, const Range& c);
  Tensor4View operator()(
      const Range& s, const Range& b, Index p, const Range& r, const Range& c);
  Tensor4View operator()(
      const Range& s, Index b, const Range& p, const Range& r, const Range& c);
  Tensor4View operator()(
      Index s, const Range& b, const Range& p, const Range& r, const Range& c);

  Tensor3View operator()(
      const Range& s, const Range& b, const Range& p, Index r, Index c);
  Tensor3View operator()(
      const Range& s, const Range& b, Index p, const Range& r, Index c);
  Tensor3View operator()(
      const Range& s, const Range& b, Index p, Index r, const Range& c);
  Tensor3View operator()(
      const Range& s, Index b, const Range& p, Index r, const Range& c);
  Tensor3View operator()(
      const Range& s, Index b, const Range& p, const Range& r, Index c);
  Tensor3View operator()(
      const Range& s, Index b, Index p, const Range& r, const Range& c);
  Tensor3View operator()(
      Index s, const Range& b, Index p, const Range& r, const Range& c);
  Tensor3View operator()(
      Index s, const Range& b, const Range& p, Index r, const Range& c);
  Tensor3View operator()(
      Index s, const Range& b, const Range& p, const Range& r, Index c);
  Tensor3View operator()(
      Index s, Index b, const Range& p, const Range& r, const Range& c);

  MatrixView operator()(
      const Range& s, const Range& b, Index p, Index r, Index c);
  MatrixView operator()(
      const Range& s, Index b, const Range& p, Index r, Index c);
  MatrixView operator()(
      const Range& s, Index b, Index p, const Range& r, Index c);
  MatrixView operator()(
      const Range& s, Index b, Index p, Index r, const Range& c);
  MatrixView operator()(
      Index s, const Range& b, Index p, Index r, const Range& c);
  MatrixView operator()(
      Index s, const Range& b, Index p, const Range& r, Index c);
  MatrixView operator()(
      Index s, const Range& b, const Range& p, Index r, Index c);
  MatrixView operator()(
      Index s, Index b, const Range& p, const Range& r, Index c);
  MatrixView operator()(
      Index s, Index b, const Range& p, Index r, const Range& c);
  MatrixView operator()(
      Index s, Index b, Index p, const Range& r, const Range& c);

  VectorView operator()(const Range& s, Index b, Index p, Index r, Index c);
  VectorView operator()(Index s, const Range& b, Index p, Index r, Index c);
  VectorView operator()(Index s, Index b, const Range& p, Index r, Index c);
  VectorView operator()(Index s, Index b, Index p, const Range& r, Index c);
  VectorView operator()(Index s, Index b, Index p, Index r, const Range& c);

#define GETFUN(s, b, p, r, c)                                                  \
  *(mdata + msr.mstart + s * msr.mstride + mbr.mstart + b * mbr.mstride +      \
    mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride + mcr.mstart + \
    c * mcr.mstride)
  /** Plain const index operator. */
  Numeric& operator()(Index s,
                      Index b,
                      Index p,
                      Index r,
                      Index c) {  // Check if indices are valid:
    ARTS_ASSERT(0 <= s);
    ARTS_ASSERT(0 <= b);
    ARTS_ASSERT(0 <= p);
    ARTS_ASSERT(0 <= r);
    ARTS_ASSERT(0 <= c);
    ARTS_ASSERT(s < msr.mextent);
    ARTS_ASSERT(b < mbr.mextent);
    ARTS_ASSERT(p < mpr.mextent);
    ARTS_ASSERT(r < mrr.mextent);
    ARTS_ASSERT(c < mcr.mextent);

    return GETFUN(s, b, p, r, c);
  }

  /** Get element implementation without assertions. */
  Numeric& get(Index s, Index b, Index p, Index r, Index c) {
    return GETFUN(s, b, p, r, c);
  }
#undef GETFUN

  // Conversion to a plain C-array
  [[nodiscard]] const Numeric* get_c_array() const ARTS_NOEXCEPT;
  Numeric* get_c_array() ARTS_NOEXCEPT;

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
  ~Tensor5View() override = default;

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
  Tensor5View() = default;
  Tensor5View(Numeric* data,
              const Range& s,
              const Range& b,
              const Range& p,
              const Range& r,
              const Range& c);
  Tensor5View(Numeric* data,
              const Range& ps,
              const Range& pb,
              const Range& pp,
              const Range& pr,
              const Range& pc,
              const Range& ns,
              const Range& nb,
              const Range& np,
              const Range& nr,
              const Range& nc);
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
  Tensor5() = default;
  Tensor5(Index s, Index b, Index p, Index r, Index c);
  Tensor5(Index s, Index b, Index p, Index r, Index c, Numeric fill);
  Tensor5(const ConstTensor5View& v);
  Tensor5(const Tensor5& v);
  Tensor5(Tensor5&& v) noexcept : Tensor5View(std::forward<Tensor5View>(v)) {
    v.mdata = nullptr;
  }

  /** Initialization from a tensor type. */
  explicit Tensor5(const matpack::tensor5_like_not_tensor5 auto &init)
      : Tensor5(matpack::shelf_size(init), matpack::book_size(init),
                matpack::page_size(init), matpack::row_size(init),
                matpack::column_size(init)) {
    auto [I, J, K, L, M] = shape().data;
    for (Index i = 0; i < I; i++)
      for (Index j = 0; j < J; j++)
        for (Index k = 0; k < K; k++)
          for (Index x = 0; x < L; x++)
            for (Index m = 0; m < M; m++)
              operator()(i, j, k, x, m) = init(i, j, k, x, m);
  }

  /** Initialization from a matrix type. */
   Tensor5& operator=(const matpack::tensor5_like_not_tensor5 auto& init) {
    return *this = Tensor5(init);
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
   */
  Tensor5(Numeric* d,
          const Range& r0,
          const Range& r1,
          const Range& r2,
          const Range& r3,
          const Range& r4) ARTS_NOEXCEPT : Tensor5View(d, r0, r1, r2, r3, r4) {
    ARTS_ASSERT(not(r0.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r1.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r2.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r3.get_extent() < 0), "Must have size");
    ARTS_ASSERT(not(r4.get_extent() < 0), "Must have size");
  }

  // Assignment operators:
  Tensor5& operator=(const Tensor5& x);
  Tensor5& operator=(Tensor5&& x) ARTS_NOEXCEPT;
  Tensor5& operator=(Numeric x);

  // Resize function:
  void resize(Index s, Index b, Index p, Index r, Index c);

  // Swap function:
  friend void swap(Tensor5& t1, Tensor5& t2) noexcept;

  // Destructor:
  ~Tensor5() noexcept override;

  /*! Reduce a Tensor5 to a Vector and leave this in an empty state */
  template <std::size_t dim0>
      Vector reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim0 < 5, "Bad Dimension, Out-of-Bounds");

    Range r0(0,
             dim0 == 0   ? nshelves()
             : dim0 == 1 ? nbooks()
             : dim0 == 2 ? npages()
             : dim0 == 3 ? nrows()
                         : ncols());

    Vector out(mdata, r0);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor5 to a Matrix and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1>
      Matrix reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim1 < 5, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");

    const Range r1(0,
                   dim1 == 1   ? nbooks()
                   : dim1 == 2 ? npages()
                   : dim1 == 3 ? nrows()
                               : ncols());
    const Range r0(0,
                   dim0 == 0   ? nshelves()
                   : dim0 == 1 ? nbooks()
                   : dim0 == 2 ? npages()
                               : nrows(),
                   r1.get_extent());

    Matrix out(mdata, r0, r1);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor5 to a Tensor3 and leave this in an empty state */
  template <std::size_t dim0, std::size_t dim1, std::size_t dim2>
      Tensor3 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim2 < 5, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");

    const Range r2(0, dim2 == 2 ? npages() : dim2 == 3 ? nrows() : ncols());
    const Range r1(0,
                   dim1 == 1   ? nbooks()
                   : dim1 == 2 ? npages()
                               : nrows(),
                   r2.get_extent());
    const Range r0(0,
                   dim0 == 0   ? nshelves()
                   : dim0 == 1 ? nbooks()
                               : npages(),
                   r1.get_extent() * r2.get_extent());

    Tensor3 out(mdata, r0, r1, r2);
    ARTS_ASSERT(size() == out.size(),
                "Can only reduce size on same size input");
    mdata = nullptr;
    return out;
  }

  /*! Reduce a Tensor5 to a Tensor4 and leave this in an empty state */
  template <std::size_t dim0,
            std::size_t dim1,
            std::size_t dim2,
            std::size_t dim3>
      Tensor4 reduce_rank() && ARTS_NOEXCEPT {
    static_assert(dim3 < 5, "Bad Dimension, Out-of-Bounds");
    static_assert(dim0 < dim1, "Bad Dimensions, dim1 must be larger than dim0");
    static_assert(dim1 < dim2, "Bad Dimensions, dim2 must be larger than dim1");
    static_assert(dim2 < dim3, "Bad Dimensions, dim3 must be larger than dim2");

    const Range r3(0, dim3 == 3 ? nrows() : ncols());
    const Range r2(0, dim2 == 2 ? npages() : nrows(), r3.get_extent());
    const Range r1(
        0, dim1 == 1 ? nbooks() : npages(), r2.get_extent() * r3.get_extent());
    const Range r0(0,
                   dim0 == 0 ? nshelves() : nbooks(),
                   r1.get_extent() * r2.get_extent() * r3.get_extent());

    Tensor4 out(mdata, r0, r1, r2, r3);
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

void copy(ConstIterator5D origin,
          const ConstIterator5D& end,
          Iterator5D target);

void copy(Numeric x, Iterator5D target, const Iterator5D& end);

void transform(Tensor5View y, double (&my_func)(double), ConstTensor5View x);

Numeric max(const ConstTensor5View& x);

Numeric min(const ConstTensor5View& x);

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

Numeric debug_tensor5view_get_elem(
    Tensor5View& tv, Index s, Index b, Index p, Index r, Index c);

#endif
////////////////////////////////

/** An array of Tensor5. */
using ArrayOfTensor5 = Array<Tensor5>;

using ArrayOfArrayOfTensor5 = Array<ArrayOfTensor5>;

#endif  // matpackV_h
