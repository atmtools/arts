/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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
   \file   matpackIV.cc

   \author Oliver Lemke
   \date   2002-11-21
*/

#include "matpackIV.h"

#include "exceptions.h"

/** The -> operator is needed, so that we can write i->begin() to get
    the 3D iterators. */
Tensor3View* Iterator4D::operator->() { return &msv; }

/** Dereferencing. */
Tensor3View& Iterator4D::operator*() { return msv; }

/** The -> operator is needed, so that we can write i->begin() to get
    the 3D iterators. */
const ConstTensor3View* ConstIterator4D::operator->() const { return &msv; }

/** Dereferencing. */
const ConstTensor3View& ConstIterator4D::operator*() const { return msv; }

// Functions for ConstTensor4View:
// ------------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor4. This allows
    correct recursive behavior.  */
ConstTensor4View ConstTensor4View::operator()(const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  return ConstTensor4View(mdata, mbr, mpr, mrr, mcr, b, p, r, c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
ConstTensor3View ConstTensor4View::operator()(const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstTensor3View(
      mdata + mcr.mstart + c * mcr.mstride, mbr, mpr, mrr, b, p, r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
ConstTensor3View ConstTensor4View::operator()(const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstTensor3View(
      mdata + mrr.mstart + r * mrr.mstride, mbr, mpr, mcr, b, p, c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
ConstTensor3View ConstTensor4View::operator()(const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  // Check that p is valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstTensor3View(
      mdata + mpr.mstart + p * mpr.mstride, mbr, mrr, mcr, b, r, c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
ConstTensor3View ConstTensor4View::operator()(Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  // Check that b is valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(b < mbr.mextent);

  return ConstTensor3View(
      mdata + mbr.mstart + b * mbr.mstride, mpr, mrr, mcr, p, r, c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
ConstMatrixView ConstTensor4View::operator()(const Range& b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  // Check that r and c are valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(
      mdata + mrr.mstart + r * mrr.mstride + mcr.mstart + c * mcr.mstride,
      mbr,
      mpr,
      b,
      p);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
ConstMatrixView ConstTensor4View::operator()(const Range& b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  // Check that p and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(
      mdata + mpr.mstart + p * mpr.mstride + mcr.mstart + c * mcr.mstride,
      mbr,
      mrr,
      b,
      r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
ConstMatrixView ConstTensor4View::operator()(const Range& b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  // Check that p and r are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstMatrixView(
      mdata + mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride,
      mbr,
      mcr,
      b,
      c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
ConstMatrixView ConstTensor4View::operator()(Index b,
                                             const Range& p,
                                             Index r,
                                             const Range& c) const {
  // Check that b and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstMatrixView(
      mdata + mbr.mstart + b * mbr.mstride + mrr.mstart + r * mrr.mstride,
      mpr,
      mcr,
      p,
      c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
ConstMatrixView ConstTensor4View::operator()(Index b,
                                             const Range& p,
                                             const Range& r,
                                             Index c) const {
  // Check that b and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(
      mdata + mbr.mstart + b * mbr.mstride + mcr.mstart + c * mcr.mstride,
      mpr,
      mrr,
      p,
      r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
ConstMatrixView ConstTensor4View::operator()(Index b,
                                             Index p,
                                             const Range& r,
                                             const Range& c) const {
  // Check that b and p are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstMatrixView(
      mdata + mbr.mstart + b * mbr.mstride + mpr.mstart + p * mpr.mstride,
      mrr,
      mcr,
      r,
      c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) */
ConstVectorView ConstTensor4View::operator()(const Range& b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  // Check that p, r and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + mpr.mstart + p * mpr.mstride + mrr.mstart +
                             r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                         mbr,
                         b);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) */
ConstVectorView ConstTensor4View::operator()(Index b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  // Check that b, r and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + mbr.mstart + b * mbr.mstride + mrr.mstart +
                             r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                         mpr,
                         p);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) */
ConstVectorView ConstTensor4View::operator()(Index b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  // Check that b, p and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                             p * mpr.mstride + mcr.mstart + c * mcr.mstride,
                         mrr,
                         r);
}

/** Const index operator returning an object of type
    ConstVectorView. Reducing the dimension by three.) */
ConstVectorView ConstTensor4View::operator()(Index b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  // Check that b, p and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstVectorView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                             p * mpr.mstride + mrr.mstart + r * mrr.mstride,
                         mcr,
                         c);
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  Tensor4View is not pointing to the beginning of a Tensor4 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric* Tensor4View::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT(mbr.mstart == 0 and
                  (mbr.mstride == mcr.mextent * mrr.mextent * mpr.mextent or
                   size() == 0),
              "Book ",
              mbr)
  ARTS_ASSERT(mpr.mstart == 0 and
                  (mpr.mstride == mcr.mextent * mrr.mextent or size() == 0),
              "Page ",
              mpr)
  ARTS_ASSERT(mrr.mstart == 0 and (mrr.mstride == mcr.mextent or size() == 0),
              "Row ",
              mrr)
  ARTS_ASSERT(
      mcr.mstart == 0 and (mcr.mstride == 1 or size() == 0), "Column ", mcr)
  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  Tensor4View is not pointing to the beginning of a Tensor4 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric* Tensor4View::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT(mbr.mstart == 0 and
                  (mbr.mstride == mcr.mextent * mrr.mextent * mpr.mextent or
                   size() == 0),
              "Book ",
              mbr)
  ARTS_ASSERT(mpr.mstart == 0 and
                  (mpr.mstride == mcr.mextent * mrr.mextent or size() == 0),
              "Page ",
              mpr)
  ARTS_ASSERT(mrr.mstart == 0 and (mrr.mstride == mcr.mextent or size() == 0),
              "Row ",
              mrr)
  ARTS_ASSERT(
      mcr.mstart == 0 and (mcr.mstride == 1 or size() == 0), "Column ", mcr)
  return mdata;
}

/** Return const iterator to first book. */
ConstIterator4D ConstTensor4View::begin() const {
  return ConstIterator4D(ConstTensor3View(mdata + mbr.mstart, mpr, mrr, mcr),
                         mbr.mstride);
}

/** Return const iterator behind last book. */
ConstIterator4D ConstTensor4View::end() const {
  return ConstIterator4D(
      ConstTensor3View(
          mdata + mbr.mstart + (mbr.mextent) * mbr.mstride, mpr, mrr, mcr),
      mbr.mstride);
}

/** Special constructor to make a Tensor4 view of a Tensor3. */
ConstTensor4View::ConstTensor4View(const ConstTensor3View& a)
    : mbr(0, 1, a.mpr.mextent * a.mrr.mextent * a.mcr.mextent),
      mpr(a.mpr),
      mrr(a.mrr),
      mcr(a.mcr),
      mdata(a.mdata) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor4 to initialize
    its own Tensor4View part. The page range pr must have a stride to
    account for the length of one page. The book range br must have a
    stride to account for the length of one book. */
ConstTensor4View::ConstTensor4View(Numeric* data,
                                   const Range& br,
                                   const Range& pr,
                                   const Range& rr,
                                   const Range& cr)
    : mbr(br), mpr(pr), mrr(rr), mcr(cr), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
ConstTensor4View::ConstTensor4View(Numeric* data,
                                   const Range& pb,
                                   const Range& pp,
                                   const Range& pr,
                                   const Range& pc,
                                   const Range& nb,
                                   const Range& np,
                                   const Range& nr,
                                   const Range& nc)
    : mbr(pb, nb), mpr(pp, np), mrr(pr, nr), mcr(pc, nc), mdata(data) {
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Tensor to print each book in turn. */
std::ostream& operator<<(std::ostream& os, const ConstTensor4View& v) {
  // Page iterators:
  ConstIterator4D ib = v.begin();
  const ConstIterator4D end_book = v.end();

  if (ib != end_book) {
    os << *ib;
    ++ib;
  }

  for (; ib != end_book; ++ib) {
    os << "\n\n";
    os << *ib;
  }

  return os;
}

// Functions for Tensor4View:
// -------------------------

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor4. This allows
    correct recursive behavior.  */
Tensor4View Tensor4View::operator()(const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  return Tensor4View(mdata, mbr, mpr, mrr, mcr, b, p, r, c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
Tensor3View Tensor4View::operator()(const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return Tensor3View(
      mdata + mcr.mstart + c * mcr.mstride, mbr, mpr, mrr, b, p, r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
Tensor3View Tensor4View::operator()(const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return Tensor3View(
      mdata + mrr.mstart + r * mrr.mstride, mbr, mpr, mcr, b, p, c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
Tensor3View Tensor4View::operator()(const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  // Check that p is valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(p < mpr.mextent);

  return Tensor3View(
      mdata + mpr.mstart + p * mpr.mstride, mbr, mrr, mcr, b, r, c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
Tensor3View Tensor4View::operator()(Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  // Check that b is valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(b < mbr.mextent);

  return Tensor3View(
      mdata + mbr.mstart + b * mbr.mstride, mpr, mrr, mcr, p, r, c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
MatrixView Tensor4View::operator()(const Range& b,
                                   const Range& p,
                                   Index r,
                                   Index c) {
  // Check that r and c are valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(
      mdata + mrr.mstart + r * mrr.mstride + mcr.mstart + c * mcr.mstride,
      mbr,
      mpr,
      b,
      p);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
MatrixView Tensor4View::operator()(const Range& b,
                                   Index p,
                                   const Range& r,
                                   Index c) {
  // Check that p and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(
      mdata + mpr.mstart + p * mpr.mstride + mcr.mstart + c * mcr.mstride,
      mbr,
      mrr,
      b,
      r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
MatrixView Tensor4View::operator()(const Range& b,
                                   Index p,
                                   Index r,
                                   const Range& c) {
  // Check that p and r are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return MatrixView(
      mdata + mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride,
      mbr,
      mcr,
      b,
      c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
MatrixView Tensor4View::operator()(Index b,
                                   const Range& p,
                                   Index r,
                                   const Range& c) {
  // Check that b and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return MatrixView(
      mdata + mbr.mstart + b * mbr.mstride + mrr.mstart + r * mrr.mstride,
      mpr,
      mcr,
      p,
      c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
MatrixView Tensor4View::operator()(Index b,
                                   const Range& p,
                                   const Range& r,
                                   Index c) {
  // Check that b and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(
      mdata + mbr.mstart + b * mbr.mstride + mcr.mstart + c * mcr.mstride,
      mpr,
      mrr,
      p,
      r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
MatrixView Tensor4View::operator()(Index b,
                                   Index p,
                                   const Range& r,
                                   const Range& c) {
  // Check that b and p are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return MatrixView(
      mdata + mbr.mstart + b * mbr.mstride + mpr.mstart + p * mpr.mstride,
      mrr,
      mcr,
      r,
      c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by three.) */
VectorView Tensor4View::operator()(const Range& b, Index p, Index r, Index c) {
  // Check that p, r and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + mpr.mstart + p * mpr.mstride + mrr.mstart +
                        r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                    mbr,
                    b);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by three.) */
VectorView Tensor4View::operator()(Index b, const Range& p, Index r, Index c) {
  // Check that b, r and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + mbr.mstart + b * mbr.mstride + mrr.mstart +
                        r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                    mpr,
                    p);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by three.) */
VectorView Tensor4View::operator()(Index b, Index p, const Range& r, Index c) {
  // Check that b, p and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                        p * mpr.mstride + mcr.mstart + c * mcr.mstride,
                    mrr,
                    r);
}

/** Index operator returning an object of type
    VectorView. Reducing the dimension by three.) */
VectorView Tensor4View::operator()(Index b, Index p, Index r, const Range& c) {
  // Check that b, p and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return VectorView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                        p * mpr.mstride + mrr.mstart + r * mrr.mstride,
                    mcr,
                    c);
}

/** Return iterator to first book. */
Iterator4D Tensor4View::begin() {
  return Iterator4D(Tensor3View(mdata + mbr.mstart, mpr, mrr, mcr),
                    mbr.mstride);
}

/** Return iterator behind last book. */
Iterator4D Tensor4View::end() {
  return Iterator4D(
      Tensor3View(
          mdata + mbr.mstart + (mbr.mextent) * mbr.mstride, mpr, mrr, mcr),
      mbr.mstride);
}

/** Assignment operator. This copies the data from another Tensor4View
    to this Tensor4View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor4View by
    setting its range. */
Tensor4View& Tensor4View::operator=(const ConstTensor4View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from Tensor4View to Tensor4View. This is a tricky
    one. The problem is that since Tensor4View is derived from
    ConstTensor4View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
Tensor4View& Tensor4View::operator=(const Tensor4View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a Tensor4. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
Tensor4View& Tensor4View::operator=(const Tensor4& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assigning a scalar to a Tensor4View will set all elements to this
    value. */
Tensor4View& Tensor4View::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
Tensor4View& Tensor4View::operator*=(Numeric x) {
  const Iterator4D eb = end();
  for (Iterator4D b = begin(); b != eb; ++b) {
    *b *= x;
  }
  return *this;
}

/** Division by scalar. */
Tensor4View& Tensor4View::operator/=(Numeric x) {
  const Iterator4D eb = end();
  for (Iterator4D b = begin(); b != eb; ++b) {
    *b /= x;
  }
  return *this;
}

/** Addition of scalar. */
Tensor4View& Tensor4View::operator+=(Numeric x) {
  const Iterator4D eb = end();
  for (Iterator4D b = begin(); b != eb; ++b) {
    *b += x;
  }
  return *this;
}

/** Subtraction of scalar. */
Tensor4View& Tensor4View::operator-=(Numeric x) {
  const Iterator4D eb = end();
  for (Iterator4D b = begin(); b != eb; ++b) {
    *b -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor4. */
Tensor4View& Tensor4View::operator*=(const ConstTensor4View& x) {
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator4D xb = x.begin();
  Iterator4D b = begin();
  const Iterator4D eb = end();
  for (; b != eb; ++b, ++xb) {
    *b *= *xb;
  }
  return *this;
}

/** Element-vise division by another Tensor4. */
Tensor4View& Tensor4View::operator/=(const ConstTensor4View& x) {
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator4D xb = x.begin();
  Iterator4D b = begin();
  const Iterator4D eb = end();
  for (; b != eb; ++b, ++xb) {
    *b /= *xb;
  }
  return *this;
}

/** Element-vise addition of another Tensor4. */
Tensor4View& Tensor4View::operator+=(const ConstTensor4View& x) {
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator4D xb = x.begin();
  Iterator4D b = begin();
  const Iterator4D eb = end();
  for (; b != eb; ++b, ++xb) {
    *b += *xb;
  }
  return *this;
}

/** Element-vise subtraction of another Tensor4. */
Tensor4View& Tensor4View::operator-=(const ConstTensor4View& x) {
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator4D xb = x.begin();
  Iterator4D b = begin();
  const Iterator4D eb = end();
  for (; b != eb; ++b, ++xb) {
    *b -= *xb;
  }
  return *this;
}

/** Special constructor to make a Tensor4 view of a Tensor3. */
Tensor4View::Tensor4View(const Tensor3View& a)
    : ConstTensor4View(
          a.mdata,
          Range(0, 1, a.mpr.mextent * a.mrr.mextent * a.mcr.mextent),
          a.mpr,
          a.mrr,
          a.mcr) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor4 to initialize its
    own Tensor4View part. The row range rr must have a
    stride to account for the length of one row. */
Tensor4View::Tensor4View(Numeric* data,
                         const Range& br,
                         const Range& pr,
                         const Range& rr,
                         const Range& cr)
    : ConstTensor4View(data, br, pr, rr, cr) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pb Previous range.
    \param pp Previous range.
    \param pr Previous range.
    \param pc Previous range.
    \param nb New Range.
    \param np New Range.
    \param nr New Range.
    \param nc New Range.
  */
Tensor4View::Tensor4View(Numeric* data,
                         const Range& pb,
                         const Range& pp,
                         const Range& pr,
                         const Range& pc,
                         const Range& nb,
                         const Range& np,
                         const Range& nr,
                         const Range& nc)
    : ConstTensor4View(data, pb, pp, pr, pc, nb, np, nr, nc) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
void copy(ConstIterator4D origin,
          const ConstIterator4D& end,
          Iterator4D target) {
  for (; origin != end; ++origin, ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(origin->begin(), origin->end(), target->begin());
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, Iterator4D target, const Iterator4D& end) {
  for (; target != end; ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(x, target->begin(), target->end());
  }
}

// Functions for Tensor4:
// ---------------------

/** Constructor setting size. This constructor has to set the strides
    in the book, page and row ranges correctly! */
Tensor4::Tensor4(Index b, Index p, Index r, Index c)
    : Tensor4View(new Numeric[b * p * r * c],
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Tensor4::Tensor4(Index b, Index p, Index r, Index c, Numeric fill)
    : Tensor4View(new Numeric[b * p * r * c],
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, b * p * r * c, fill);
}

/** Copy constructor from Tensor4View. This automatically sets the size
    and copies the data. */
Tensor4::Tensor4(const ConstTensor4View& m)
    : Tensor4View(new Numeric[m.nbooks() * m.npages() * m.nrows() * m.ncols()],
                  Range(0, m.nbooks(), m.npages() * m.nrows() * m.ncols()),
                  Range(0, m.npages(), m.nrows() * m.ncols()),
                  Range(0, m.nrows(), m.ncols()),
                  Range(0, m.ncols())) {
  copy(m.begin(), m.end(), begin());
}

/** Copy constructor from Tensor4. This automatically sets the size
    and copies the data. */
Tensor4::Tensor4(const Tensor4& m)
    : Tensor4View(new Numeric[m.nbooks() * m.npages() * m.nrows() * m.ncols()],
                  Range(0, m.nbooks(), m.npages() * m.nrows() * m.ncols()),
                  Range(0, m.npages(), m.nrows() * m.ncols()),
                  Range(0, m.nrows(), m.ncols()),
                  Range(0, m.ncols())) {
  // There is a catch here: If m is an empty tensor, then it will have
  // dimensions of size 0. But these are used to initialize the stride
  // for higher dimensions! Thus, this method has to be consistent
  // with the behaviour of Range::Range. For now, Range::Range allows
  // also stride 0.
  std::memcpy(mdata,
              m.mdata,
              nbooks() * npages() * nrows() * ncols() * sizeof(Numeric));
}

//! Assignment operator from another tensor.
/*! 
  While dimensions of views can not be adjusted, dimensions of
  tensors *can* be adjusted. Hence, the behavior of the assignment
  operator is different.

  In this case the size of the target is automatically adjusted. This
  is important, so that structures containing tensors are copied
  correctly. 
  
  This is a deviation from the old ARTS paradigm that sizes must match
  exactly before copying!

  Note: It is sufficient to have only this one version of the
  assignment (Tensor = Tensor). It implicitly covers the cases
  Tensor=TensorView, etc, because there is a default constructor for
  Tensor from TensorView. (See C++ Primer Plus, page 571ff.)

  \param m The other tensor to assign to this one.
  \return This tensor, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
Tensor4& Tensor4::operator=(const Tensor4& x) {
  if (this != &x) {
    resize(x.nbooks(), x.npages(), x.nrows(), x.ncols());
    std::memcpy(mdata,
                x.mdata,
                nbooks() * npages() * nrows() * ncols() * sizeof(Numeric));
  }
  return *this;
}

//! Move assignment operator from another tensor.
Tensor4& Tensor4::operator=(Tensor4&& x) noexcept {
  if (this != &x) {
    delete[] mdata;
    mdata = x.mdata;
    mbr = x.mbr;
    mpr = x.mpr;
    mrr = x.mrr;
    mcr = x.mcr;
    x.mbr = Range(0, 0);
    x.mpr = Range(0, 0);
    x.mrr = Range(0, 0);
    x.mcr = Range(0, 0);
    x.mdata = nullptr;
  }
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
Tensor4& Tensor4::operator=(Numeric x) {
  std::fill_n(mdata, nbooks() * npages() * nrows() * ncols(), x);
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
void Tensor4::resize(Index b, Index p, Index r, Index c) {
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  if (mbr.mextent != b || mpr.mextent != p || mrr.mextent != r ||
      mcr.mextent != c) {
    delete[] mdata;
    mdata = new Numeric[b * p * r * c];

    mbr.mstart = 0;
    mbr.mextent = b;
    mbr.mstride = p * r * c;

    mpr.mstart = 0;
    mpr.mextent = p;
    mpr.mstride = r * c;

    mrr.mstart = 0;
    mrr.mextent = r;
    mrr.mstride = c;

    mcr.mstart = 0;
    mcr.mextent = c;
    mcr.mstride = 1;
  }
}

/** Swaps two objects. */
void swap(Tensor4& t1, Tensor4& t2) noexcept {
  using std::swap;
  swap(t1.mbr, t2.mbr);
  swap(t1.mpr, t2.mpr);
  swap(t1.mrr, t2.mrr);
  swap(t1.mcr, t2.mcr);
  swap(t1.mdata, t2.mdata);
}

/** Destructor for Tensor4. This is important, since Tensor4 uses new to
    allocate storage. */
Tensor4::~Tensor4() {
  //   cout << "Destroying a Tensor4:\n"
  //        << *this << "\n........................................\n";
  delete[] mdata;
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
void transform(Tensor4View y, double (&my_func)(double), ConstTensor4View x) {
  // Check dimensions:
  ARTS_ASSERT(y.nbooks() == x.nbooks());
  ARTS_ASSERT(y.npages() == x.npages());
  ARTS_ASSERT(y.nrows() == x.nrows());
  ARTS_ASSERT(y.ncols() == x.ncols());

  const ConstIterator4D xe = x.end();
  ConstIterator4D xi = x.begin();
  Iterator4D yi = y.begin();
  for (; xi != xe; ++xi, ++yi) {
    // Use the transform function of lower dimensional tensors
    // recursively:
    transform(*yi, my_func, *xi);
  }
}

/** Max function, tensor version. */
Numeric max(const ConstTensor4View& x) {
  const ConstIterator4D xe = x.end();
  ConstIterator4D xi = x.begin();

  // Initial value for max:
  Numeric themax = max(*xi);
  ++xi;

  for (; xi != xe; ++xi) {
    // Use the max function of lower dimensional tensors
    // recursively:
    Numeric maxi = max(*xi);
    if (maxi > themax) themax = maxi;
  }

  return themax;
}

/** Min function, tensor version. */
Numeric min(const ConstTensor4View& x) {
  const ConstIterator4D xe = x.end();
  ConstIterator4D xi = x.begin();

  // Initial value for min:
  Numeric themin = min(*xi);
  ++xi;

  for (; xi != xe; ++xi) {
    // Use the min function of lower dimensional tensors
    // recursively:
    Numeric mini = min(*xi);
    if (mini < themin) themin = mini;
  }

  return themin;
}

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

/** Helper function to access tensor elements.

    Because of function inlining the operator() is not
    accessible from the debuggger. This function helps to access
    Tensor elements from within the debugger.

    \param tv TensorView
    \param b  Book index
    \param p  Page index
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Numeric debug_tensor4view_get_elem(
    Tensor4View& tv, Index b, Index p, Index r, Index c) {
  return tv(b, p, r, c);
}

#endif
////////////////////////////////

Vector Tensor4::flatten() && ARTS_NOEXCEPT {
  Vector out(mdata, Range(0, size()));
  mdata = nullptr;
  return out;
}
