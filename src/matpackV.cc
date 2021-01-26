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
   \file   matpackV.cc

   \author Oliver Lemke
   \date   2002-11-21
*/

#include "matpackV.h"
#include "exceptions.h"

using std::runtime_error;

// Functions for ConstTensor5View:
// ------------------------------

//! Check if variable is empty.
/*!
 \param[in]  x The variable to check.
 \return True if the size of any dimension of x is 0.
 */
bool ConstTensor5View::empty() const {
  return (nshelves() == 0 || nbooks() == 0 || npages() == 0 || nrows() == 0 ||
          ncols() == 0);
}

/** Returns the number of shelves. */
Index ConstTensor5View::nshelves() const { return msr.mextent; }

/** Returns the number of books. */
Index ConstTensor5View::nbooks() const { return mbr.mextent; }

/** Returns the number of pages. */
Index ConstTensor5View::npages() const { return mpr.mextent; }

/** Returns the number of rows. */
Index ConstTensor5View::nrows() const { return mrr.mextent; }

/** Returns the number of columns. */
Index ConstTensor5View::ncols() const { return mcr.mextent; }

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor5. This allows
    correct recursive behavior.  */
ConstTensor5View ConstTensor5View::operator()(const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  return ConstTensor5View(mdata, msr, mbr, mpr, mrr, mcr, s, b, p, r, c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstTensor4View(
      mdata + mcr.mstart + c * mcr.mstride, msr, mbr, mpr, mrr, s, b, p, r);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstTensor4View(
      mdata + mrr.mstart + r * mrr.mstride, msr, mbr, mpr, mcr, s, b, p, c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  // Check that p is valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstTensor4View(
      mdata + mpr.mstart + p * mpr.mstride, msr, mbr, mrr, mcr, s, b, r, c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
ConstTensor4View ConstTensor5View::operator()(const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  // Check that b is valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(b < mbr.mextent);

  return ConstTensor4View(
      mdata + mbr.mstart + b * mbr.mstride, msr, mpr, mrr, mcr, s, p, r, c);
}

/** Const index operator returning an object of type
    ConstTensor4View. (Reducing the dimension by one.) */
ConstTensor4View ConstTensor5View::operator()(Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  // Check that s is valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(s < msr.mextent);

  return ConstTensor4View(
      mdata + msr.mstart + s * msr.mstride, mbr, mpr, mrr, mcr, b, p, r, c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    const Range& s, const Range& b, const Range& p, Index r, Index c) const {
  // Check that r and c is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstTensor3View(
      mdata + mrr.mstart + r * mrr.mstride + mcr.mstart + c * mcr.mstride,
      msr,
      mbr,
      mpr,
      s,
      b,
      p);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    const Range& s, const Range& b, Index p, const Range& r, Index c) const {
  // Check that p and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstTensor3View(
      mdata + mpr.mstart + p * mpr.mstride + mcr.mstart + c * mcr.mstride,
      msr,
      mbr,
      mrr,
      s,
      b,
      r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    const Range& s, const Range& b, Index p, Index r, const Range& c) const {
  // Check that p and r are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstTensor3View(
      mdata + mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride,
      msr,
      mbr,
      mcr,
      s,
      b,
      c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    const Range& s, Index b, const Range& p, Index r, const Range& c) const {
  // Check that b and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstTensor3View(
      mdata + mbr.mstart + b * mbr.mstride + mrr.mstart + r * mrr.mstride,
      msr,
      mpr,
      mcr,
      s,
      p,
      c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    const Range& s, Index b, const Range& p, const Range& r, Index c) const {
  // Check that b and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstTensor3View(
      mdata + mbr.mstart + b * mbr.mstride + mcr.mstart + c * mcr.mstride,
      msr,
      mpr,
      mrr,
      s,
      p,
      r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    const Range& s, Index b, Index p, const Range& r, const Range& c) const {
  // Check that b and p are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstTensor3View(
      mdata + mbr.mstart + b * mbr.mstride + mpr.mstart + p * mpr.mstride,
      msr,
      mrr,
      mcr,
      s,
      r,
      c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    Index s, const Range& b, Index p, const Range& r, const Range& c) const {
  // Check that s and p are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstTensor3View(
      mdata + msr.mstart + s * msr.mstride + mpr.mstart + p * mpr.mstride,
      mbr,
      mrr,
      mcr,
      b,
      r,
      c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    Index s, const Range& b, const Range& p, Index r, const Range& c) const {
  // Check that s and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstTensor3View(
      mdata + msr.mstart + s * msr.mstride + mrr.mstart + r * mrr.mstride,
      mbr,
      mpr,
      mcr,
      b,
      p,
      c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    Index s, const Range& b, const Range& p, const Range& r, Index c) const {
  // Check that s and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstTensor3View(
      mdata + msr.mstart + s * msr.mstride + mcr.mstart + c * mcr.mstride,
      mbr,
      mpr,
      mrr,
      b,
      p,
      r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by two.) */
ConstTensor3View ConstTensor5View::operator()(
    Index s, Index b, const Range& p, const Range& r, const Range& c) const {
  // Check that s and b are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);

  return ConstTensor3View(
      mdata + msr.mstart + s * msr.mstride + mbr.mstart + b * mbr.mstride,
      mpr,
      mrr,
      mcr,
      p,
      r,
      c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    const Range& s, const Range& b, Index p, Index r, Index c) const {
  // Check that p, r and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + mpr.mstart + p * mpr.mstride + mrr.mstart +
                             r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                         msr,
                         mbr,
                         s,
                         b);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    const Range& s, Index b, const Range& p, Index r, Index c) const {
  // Check that b, r and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + mbr.mstart + b * mbr.mstride + mrr.mstart +
                             r * mpr.mstride + mcr.mstart + c * mcr.mstride,
                         msr,
                         mpr,
                         s,
                         p);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    const Range& s, Index b, Index p, const Range& r, Index c) const {
  // Check that b, p and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                             p * mpr.mstride + mcr.mstart + c * mcr.mstride,
                         msr,
                         mrr,
                         s,
                         r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    const Range& s, Index b, Index p, Index r, const Range& c) const {
  // Check that b, p and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstMatrixView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                             p * mpr.mstride + mrr.mstart + r * mrr.mstride,
                         msr,
                         mcr,
                         s,
                         c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    Index s, const Range& b, Index p, Index r, const Range& c) const {
  // Check that s, p and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstMatrixView(mdata + msr.mstart + s * msr.mstride + mpr.mstart +
                             p * mpr.mstride + mrr.mstart + r * mrr.mstride,
                         mbr,
                         mcr,
                         b,
                         c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    Index s, const Range& b, Index p, const Range& r, Index c) const {
  // Check that s, p and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + msr.mstart + s * msr.mstride + mpr.mstart +
                             p * mpr.mstride + mcr.mstart + c * mcr.mstride,
                         mbr,
                         mrr,
                         b,
                         r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    Index s, const Range& b, const Range& p, Index r, Index c) const {
  // Check that s, r and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + msr.mstart + s * msr.mstride + mrr.mstart +
                             r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                         mbr,
                         mpr,
                         b,
                         p);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    Index s, Index b, const Range& p, const Range& r, Index c) const {
  // Check that s, b and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                             b * mbr.mstride + mcr.mstart + c * mcr.mstride,
                         mpr,
                         mrr,
                         p,
                         r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    Index s, Index b, const Range& p, Index r, const Range& c) const {
  // Check that s, b and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstMatrixView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                             b * mbr.mstride + mrr.mstart + r * mrr.mstride,
                         mpr,
                         mcr,
                         p,
                         c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by three.) */
ConstMatrixView ConstTensor5View::operator()(
    Index s, Index b, Index p, const Range& r, const Range& c) const {
  // Check that s, b and p are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstMatrixView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                             b * mbr.mstride + mpr.mstart + p * mpr.mstride,
                         mrr,
                         mcr,
                         r,
                         c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
ConstVectorView ConstTensor5View::operator()(
    const Range& s, Index b, Index p, Index r, Index c) const {
  // Check that b, p, r and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                             p * mpr.mstride + mrr.mstart + r * mrr.mstride +
                             mcr.mstart + c * mcr.mstride,
                         msr,
                         s);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
ConstVectorView ConstTensor5View::operator()(
    Index s, const Range& b, Index p, Index r, Index c) const {
  // Check that s, p, r and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + msr.mstart + s * msr.mstride + mpr.mstart +
                             p * mpr.mstride + mrr.mstart + r * mrr.mstride +
                             mcr.mstart + c * mcr.mstride,
                         mbr,
                         b);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
ConstVectorView ConstTensor5View::operator()(
    Index s, Index b, const Range& p, Index r, Index c) const {
  // Check that s, b, r and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                             b * mbr.mstride + mrr.mstart + r * mrr.mstride +
                             mcr.mstart + c * mcr.mstride,
                         mpr,
                         p);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
ConstVectorView ConstTensor5View::operator()(
    Index s, Index b, Index p, const Range& r, Index c) const {
  // Check that s, b, p and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                             b * mbr.mstride + mpr.mstart + p * mpr.mstride +
                             mcr.mstart + c * mcr.mstride,
                         mrr,
                         r);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
ConstVectorView ConstTensor5View::operator()(
    Index s, Index b, Index p, Index r, const Range& c) const {
  // Check that s, b, p and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstVectorView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                             b * mbr.mstride + mpr.mstart + p * mpr.mstride +
                             mrr.mstart + r * mrr.mstride,
                         mcr,
                         c);
}

/** Conversion to plain C-array.
    
This function returns a pointer to the raw data. It fails if the
Tensor5View is not pointing to the beginning of a Tensor5 or the stride
is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric* Tensor5View::get_c_array() {
  if (msr.mstart != 0 ||
      msr.mstride != mbr.mextent * mpr.mextent * mrr.mextent * mcr.mextent ||
      mbr.mstart != 0 ||
      mbr.mstride != mpr.mextent * mrr.mextent * mcr.mextent ||
      mpr.mstart != 0 || mpr.mstride != mrr.mextent * mcr.mextent ||
      mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1)
    throw std::runtime_error(
        "A Tensor5View can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  Tensor5View is not pointing to the beginning of a Tensor5 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric* Tensor5View::get_c_array() const {
  if (msr.mstart != 0 ||
      msr.mstride != mbr.mextent * mpr.mextent * mrr.mextent * mcr.mextent ||
      mbr.mstart != 0 ||
      mbr.mstride != mpr.mextent * mrr.mextent * mcr.mextent ||
      mpr.mstart != 0 || mpr.mstride != mrr.mextent * mcr.mextent ||
      mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1)
    throw std::runtime_error(
        "A Tensor5View can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** Return const iterator to first shelf. */
ConstIterator5D ConstTensor5View::begin() const {
  return ConstIterator5D(
      ConstTensor4View(mdata + msr.mstart, mbr, mpr, mrr, mcr), msr.mstride);
}

/** Return const iterator behind last shelf. */
ConstIterator5D ConstTensor5View::end() const {
  return ConstIterator5D(
      ConstTensor4View(
          mdata + msr.mstart + (msr.mextent) * msr.mstride, mbr, mpr, mrr, mcr),
      msr.mstride);
}

/** Special constructor to make a Tensor5 view of a Tensor4. */
ConstTensor5View::ConstTensor5View(const ConstTensor4View& a)
    : msr(0, 1, a.mbr.mextent * a.mpr.mextent * a.mrr.mextent * a.mcr.mextent),
      mbr(a.mbr),
      mpr(a.mpr),
      mrr(a.mrr),
      mcr(a.mcr),
      mdata(a.mdata) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor5 to initialize
    its own Tensor5View part. The book range br must have a stride to
    account for the length of one book. The shelf range sr must have a
    stride to account for the length of one shelf. */
ConstTensor5View::ConstTensor5View(Numeric* data,
                                   const Range& sr,
                                   const Range& br,
                                   const Range& pr,
                                   const Range& rr,
                                   const Range& cr)
    : msr(sr), mbr(br), mpr(pr), mrr(rr), mcr(cr), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
ConstTensor5View::ConstTensor5View(Numeric* data,
                                   const Range& ps,
                                   const Range& pb,
                                   const Range& pp,
                                   const Range& pr,
                                   const Range& pc,
                                   const Range& ns,
                                   const Range& nb,
                                   const Range& np,
                                   const Range& nr,
                                   const Range& nc)
    : msr(ps, ns),
      mbr(pb, nb),
      mpr(pp, np),
      mrr(pr, nr),
      mcr(pc, nc),
      mdata(data) {
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Tensor to print each book in turn. */
std::ostream& operator<<(std::ostream& os, const ConstTensor5View& v) {
  // Page iterators:
  ConstIterator5D is = v.begin();
  const ConstIterator5D end_shelf = v.end();

  if (is != end_shelf) {
    os << *is;
    ++is;
  }

  for (; is != end_shelf; ++is) {
    os << "\n\n";
    os << *is;
  }

  return os;
}

// Functions for Tensor5View:
// -------------------------

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor5. This allows
    correct recursive behavior.  */
Tensor5View Tensor5View::operator()(const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  return Tensor5View(mdata, msr, mbr, mpr, mrr, mcr, s, b, p, r, c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
Tensor4View Tensor5View::operator()(
    const Range& s, const Range& b, const Range& p, const Range& r, Index c) {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return Tensor4View(
      mdata + mcr.mstart + c * mcr.mstride, msr, mbr, mpr, mrr, s, b, p, r);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
Tensor4View Tensor5View::operator()(
    const Range& s, const Range& b, const Range& p, Index r, const Range& c) {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return Tensor4View(
      mdata + mrr.mstart + r * mrr.mstride, msr, mbr, mpr, mcr, s, b, p, c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
Tensor4View Tensor5View::operator()(
    const Range& s, const Range& b, Index p, const Range& r, const Range& c) {
  // Check that p is valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(p < mpr.mextent);

  return Tensor4View(
      mdata + mpr.mstart + p * mpr.mstride, msr, mbr, mrr, mcr, s, b, r, c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
Tensor4View Tensor5View::operator()(
    const Range& s, Index b, const Range& p, const Range& r, const Range& c) {
  // Check that b is valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(b < mbr.mextent);

  return Tensor4View(
      mdata + mbr.mstart + b * mbr.mstride, msr, mpr, mrr, mcr, s, p, r, c);
}

/** Index operator returning an object of type
    Tensor4View. (Reducing the dimension by one.) */
Tensor4View Tensor5View::operator()(
    Index s, const Range& b, const Range& p, const Range& r, const Range& c) {
  // Check that s is valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(s < msr.mextent);

  return Tensor4View(
      mdata + msr.mstart + s * msr.mstride, mbr, mpr, mrr, mcr, b, p, r, c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    const Range& s, const Range& b, const Range& p, Index r, Index c) {
  // Check that r and c is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return Tensor3View(
      mdata + mrr.mstart + r * mrr.mstride + mcr.mstart + c * mcr.mstride,
      msr,
      mbr,
      mpr,
      s,
      b,
      p);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    const Range& s, const Range& b, Index p, const Range& r, Index c) {
  // Check that p and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return Tensor3View(
      mdata + mpr.mstart + p * mpr.mstride + mcr.mstart + c * mcr.mstride,
      msr,
      mbr,
      mrr,
      s,
      b,
      r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    const Range& s, const Range& b, Index p, Index r, const Range& c) {
  // Check that p and r are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return Tensor3View(
      mdata + mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride,
      msr,
      mbr,
      mcr,
      s,
      b,
      c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    const Range& s, Index b, const Range& p, Index r, const Range& c) {
  // Check that b and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return Tensor3View(
      mdata + mbr.mstart + b * mbr.mstride + mrr.mstart + r * mrr.mstride,
      msr,
      mpr,
      mcr,
      s,
      p,
      c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    const Range& s, Index b, const Range& p, const Range& r, Index c) {
  // Check that b and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return Tensor3View(
      mdata + mbr.mstart + b * mbr.mstride + mcr.mstart + c * mcr.mstride,
      msr,
      mpr,
      mrr,
      s,
      p,
      r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    const Range& s, Index b, Index p, const Range& r, const Range& c) {
  // Check that b and p are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return Tensor3View(
      mdata + mbr.mstart + b * mbr.mstride + mpr.mstart + p * mpr.mstride,
      msr,
      mrr,
      mcr,
      s,
      r,
      c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    Index s, const Range& b, Index p, const Range& r, const Range& c) {
  // Check that s and p are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return Tensor3View(
      mdata + msr.mstart + s * msr.mstride + mpr.mstart + p * mpr.mstride,
      mbr,
      mrr,
      mcr,
      b,
      r,
      c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    Index s, const Range& b, const Range& p, Index r, const Range& c) {
  // Check that s and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return Tensor3View(
      mdata + msr.mstart + s * msr.mstride + mrr.mstart + r * mrr.mstride,
      mbr,
      mpr,
      mcr,
      b,
      p,
      c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    Index s, const Range& b, const Range& p, const Range& r, Index c) {
  // Check that s and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return Tensor3View(
      mdata + msr.mstart + s * msr.mstride + mcr.mstart + c * mcr.mstride,
      mbr,
      mpr,
      mrr,
      b,
      p,
      r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by two.) */
Tensor3View Tensor5View::operator()(
    Index s, Index b, const Range& p, const Range& r, const Range& c) {
  // Check that s and b are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);

  return Tensor3View(
      mdata + msr.mstart + s * msr.mstride + mbr.mstart + b * mbr.mstride,
      mpr,
      mrr,
      mcr,
      p,
      r,
      c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    const Range& s, const Range& b, Index p, Index r, Index c) {
  // Check that p, r and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + mpr.mstart + p * mpr.mstride + mrr.mstart +
                        r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                    msr,
                    mbr,
                    s,
                    b);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    const Range& s, Index b, const Range& p, Index r, Index c) {
  // Check that b, r and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + mbr.mstart + b * mbr.mstride + mrr.mstart +
                        r * mpr.mstride + mcr.mstart + c * mcr.mstride,
                    msr,
                    mpr,
                    s,
                    p);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    const Range& s, Index b, Index p, const Range& r, Index c) {
  // Check that b, p and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                        p * mpr.mstride + mcr.mstart + c * mcr.mstride,
                    msr,
                    mrr,
                    s,
                    r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    const Range& s, Index b, Index p, Index r, const Range& c) {
  // Check that b, p and r are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return MatrixView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                        p * mpr.mstride + mrr.mstart + r * mrr.mstride,
                    msr,
                    mcr,
                    s,
                    c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    Index s, const Range& b, Index p, Index r, const Range& c) {
  // Check that s, p and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return MatrixView(mdata + msr.mstart + s * msr.mstride + mpr.mstart +
                        p * mpr.mstride + mrr.mstart + r * mrr.mstride,
                    mbr,
                    mcr,
                    b,
                    c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    Index s, const Range& b, Index p, const Range& r, Index c) {
  // Check that s, p and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + msr.mstart + s * msr.mstride + mpr.mstart +
                        p * mpr.mstride + mcr.mstart + c * mcr.mstride,
                    mbr,
                    mrr,
                    b,
                    r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    Index s, const Range& b, const Range& p, Index r, Index c) {
  // Check that s, r and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + msr.mstart + s * msr.mstride + mrr.mstart +
                        r * mrr.mstride + mcr.mstart + c * mcr.mstride,
                    mbr,
                    mpr,
                    b,
                    p);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    Index s, Index b, const Range& p, const Range& r, Index c) {
  // Check that s, b and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                        b * mbr.mstride + mcr.mstart + c * mcr.mstride,
                    mpr,
                    mrr,
                    p,
                    r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    Index s, Index b, const Range& p, Index r, const Range& c) {
  // Check that s, b and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return MatrixView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                        b * mbr.mstride + mrr.mstart + r * mrr.mstride,
                    mpr,
                    mcr,
                    p,
                    c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by three.) */
MatrixView Tensor5View::operator()(
    Index s, Index b, Index p, const Range& r, const Range& c) {
  // Check that s, b and p are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);

  return MatrixView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                        b * mbr.mstride + mpr.mstart + p * mpr.mstride,
                    mrr,
                    mcr,
                    r,
                    c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
VectorView Tensor5View::operator()(
    const Range& s, Index b, Index p, Index r, Index c) {
  // Check that b, p, r and c are valid:
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + mbr.mstart + b * mbr.mstride + mpr.mstart +
                        p * mpr.mstride + mrr.mstart + r * mrr.mstride +
                        mcr.mstart + c * mcr.mstride,
                    msr,
                    s);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
VectorView Tensor5View::operator()(
    Index s, const Range& b, Index p, Index r, Index c) {
  // Check that s, p, r and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + msr.mstart + s * msr.mstride + mpr.mstart +
                        p * mpr.mstride + mrr.mstart + r * mrr.mstride +
                        mcr.mstart + c * mcr.mstride,
                    mbr,
                    b);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
VectorView Tensor5View::operator()(
    Index s, Index b, const Range& p, Index r, Index c) {
  // Check that s, b, r and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                        b * mbr.mstride + mrr.mstart + r * mrr.mstride +
                        mcr.mstart + c * mcr.mstride,
                    mpr,
                    p);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by four.) */
VectorView Tensor5View::operator()(
    Index s, Index b, Index p, const Range& r, Index c) {
  // Check that s, b, p and c are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                        b * mbr.mstride + mpr.mstart + p * mpr.mstride +
                        mcr.mstart + c * mcr.mstride,
                    mrr,
                    r);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by four.) */
VectorView Tensor5View::operator()(
    Index s, Index b, Index p, Index r, const Range& c) {
  // Check that s, b, p and r are valid:
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(s < msr.mextent);
  ARTS_ASSERT(b < mbr.mextent);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return VectorView(mdata + msr.mstart + s * msr.mstride + mbr.mstart +
                        b * mbr.mstride + mpr.mstart + p * mpr.mstride +
                        mrr.mstart + r * mrr.mstride,
                    mcr,
                    c);
}

/** Return iterator to first shelf. */
Iterator5D Tensor5View::begin() {
  return Iterator5D(Tensor4View(mdata + msr.mstart, mbr, mpr, mrr, mcr),
                    msr.mstride);
}

/** Return iterator behind last shelf. */
Iterator5D Tensor5View::end() {
  return Iterator5D(
      Tensor4View(
          mdata + msr.mstart + (msr.mextent) * msr.mstride, mbr, mpr, mrr, mcr),
      msr.mstride);
}

/** Assignment operator. This copies the data from another Tensor5View
    to this Tensor5View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor5View by
    setting its range. */
Tensor5View& Tensor5View::operator=(const ConstTensor5View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(msr.mextent == m.msr.mextent);
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from Tensor5View to Tensor5View. This is a tricky
    one. The problem is that since Tensor5View is derived from
    ConstTensor5View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
Tensor5View& Tensor5View::operator=(const Tensor5View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(msr.mextent == m.msr.mextent);
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a Tensor5. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
Tensor5View& Tensor5View::operator=(const Tensor5& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(msr.mextent == m.msr.mextent);
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assigning a scalar to a Tensor5View will set all elements to this
    value. */
Tensor5View& Tensor5View::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
Tensor5View& Tensor5View::operator*=(Numeric x) {
  const Iterator5D es = end();
  for (Iterator5D s = begin(); s != es; ++s) {
    *s *= x;
  }
  return *this;
}

/** Division by scalar. */
Tensor5View& Tensor5View::operator/=(Numeric x) {
  const Iterator5D es = end();
  for (Iterator5D s = begin(); s != es; ++s) {
    *s /= x;
  }
  return *this;
}

/** Addition of scalar. */
Tensor5View& Tensor5View::operator+=(Numeric x) {
  const Iterator5D es = end();
  for (Iterator5D s = begin(); s != es; ++s) {
    *s += x;
  }
  return *this;
}

/** Subtraction of scalar. */
Tensor5View& Tensor5View::operator-=(Numeric x) {
  const Iterator5D es = end();
  for (Iterator5D s = begin(); s != es; ++s) {
    *s -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor5. */
Tensor5View& Tensor5View::operator*=(const ConstTensor5View& x) {
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator5D xs = x.begin();
  Iterator5D s = begin();
  const Iterator5D es = end();
  for (; s != es; ++s, ++xs) {
    *s *= *xs;
  }
  return *this;
}

/** Element-vise division by another Tensor5. */
Tensor5View& Tensor5View::operator/=(const ConstTensor5View& x) {
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator5D xs = x.begin();
  Iterator5D s = begin();
  const Iterator5D es = end();
  for (; s != es; ++s, ++xs) {
    *s /= *xs;
  }
  return *this;
}

/** Element-vise addition of another Tensor5. */
Tensor5View& Tensor5View::operator+=(const ConstTensor5View& x) {
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator5D xs = x.begin();
  Iterator5D s = begin();
  const Iterator5D es = end();
  for (; s != es; ++s, ++xs) {
    *s += *xs;
  }
  return *this;
}

/** Element-vise subtraction of another Tensor5. */
Tensor5View& Tensor5View::operator-=(const ConstTensor5View& x) {
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator5D xs = x.begin();
  Iterator5D s = begin();
  const Iterator5D es = end();
  for (; s != es; ++s, ++xs) {
    *s -= *xs;
  }
  return *this;
}

/** Special constructor to make a Tensor5 view of a Tensor4. */
Tensor5View::Tensor5View(const Tensor4View& a)
    : ConstTensor5View(
          a.mdata,
          Range(0,
                1,
                a.mbr.mextent * a.mpr.mextent * a.mrr.mextent * a.mcr.mextent),
          a.mbr,
          a.mpr,
          a.mrr,
          a.mcr) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor5 to initialize its
    own Tensor5View part. The row range rr must have a
    stride to account for the length of one row. */
Tensor5View::Tensor5View(Numeric* data,
                         const Range& sr,
                         const Range& br,
                         const Range& pr,
                         const Range& rr,
                         const Range& cr)
    : ConstTensor5View(data, sr, br, pr, rr, cr) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param ps Previous range.
    \param pb Previous range.
    \param pp Previous range.
    \param pr Previous range.
    \param pc Previous range.
    \param ns New Range.
    \param nb New Range.
    \param np New Range.
    \param nr New Range.
    \param nc New Range.
  */
Tensor5View::Tensor5View(Numeric* data,
                         const Range& ps,
                         const Range& pb,
                         const Range& pp,
                         const Range& pr,
                         const Range& pc,
                         const Range& ns,
                         const Range& nb,
                         const Range& np,
                         const Range& nr,
                         const Range& nc)
    : ConstTensor5View(data, ps, pb, pp, pr, pc, ns, nb, np, nr, nc) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
void copy(ConstIterator5D origin,
          const ConstIterator5D& end,
          Iterator5D target) {
  for (; origin != end; ++origin, ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(origin->begin(), origin->end(), target->begin());
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, Iterator5D target, const Iterator5D& end) {
  for (; target != end; ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(x, target->begin(), target->end());
  }
}

// Functions for Tensor5:
// ---------------------

/** Constructor setting size. This constructor has to set the strides
    in the shelf, book, page and row ranges correctly! */
Tensor5::Tensor5(Index s, Index b, Index p, Index r, Index c)
    : Tensor5View(new Numeric[s * b * p * r * c],
                  Range(0, s, b * p * r * c),
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Tensor5::Tensor5(Index s, Index b, Index p, Index r, Index c, Numeric fill)
    : Tensor5View(new Numeric[s * b * p * r * c],
                  Range(0, s, b * p * r * c),
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, s * b * p * r * c, fill);
}

/** Copy constructor from Tensor5View. This automatically sets the size
    and copies the data. */
Tensor5::Tensor5(const ConstTensor5View& m)
    : Tensor5View(
          new Numeric[m.nshelves() * m.nbooks() * m.npages() * m.nrows() *
                      m.ncols()],
          Range(
              0, m.nshelves(), m.nbooks() * m.npages() * m.nrows() * m.ncols()),
          Range(0, m.nbooks(), m.npages() * m.nrows() * m.ncols()),
          Range(0, m.npages(), m.nrows() * m.ncols()),
          Range(0, m.nrows(), m.ncols()),
          Range(0, m.ncols())) {
  copy(m.begin(), m.end(), begin());
}

/** Copy constructor from Tensor5. This automatically sets the size
    and copies the data. */
Tensor5::Tensor5(const Tensor5& m)
    : Tensor5View(
          new Numeric[m.nshelves() * m.nbooks() * m.npages() * m.nrows() *
                      m.ncols()],
          Range(
              0, m.nshelves(), m.nbooks() * m.npages() * m.nrows() * m.ncols()),
          Range(0, m.nbooks(), m.npages() * m.nrows() * m.ncols()),
          Range(0, m.npages(), m.nrows() * m.ncols()),
          Range(0, m.nrows(), m.ncols()),
          Range(0, m.ncols())) {
  // There is a catch here: If m is an empty tensor, then it will have
  // dimensions of size 0. But these are used to initialize the stride
  // for higher dimensions! Thus, this method has to be consistent
  // with the behaviour of Range::Range. For now, Range::Range allows
  // also stride 0.
  std::memcpy(
      mdata,
      m.mdata,
      nshelves() * nbooks() * npages() * nrows() * ncols() * sizeof(Numeric));
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
Tensor5& Tensor5::operator=(const Tensor5& x) {
  if (this != &x) {
    resize(x.nshelves(), x.nbooks(), x.npages(), x.nrows(), x.ncols());
    std::memcpy(
        mdata,
        x.mdata,
        nshelves() * nbooks() * npages() * nrows() * ncols() * sizeof(Numeric));
  }
  return *this;
}

//! Move assignment operator from another tensor.
Tensor5& Tensor5::operator=(Tensor5&& x) noexcept {
  if (this != &x) {
    delete[] mdata;
    mdata = x.mdata;
    msr = x.msr;
    mbr = x.mbr;
    mpr = x.mpr;
    mrr = x.mrr;
    mcr = x.mcr;
    x.msr = Range(0, 0);
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
Tensor5& Tensor5::operator=(Numeric x) {
  std::fill_n(mdata, nshelves() * nbooks() * npages() * nrows() * ncols(), x);
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
void Tensor5::resize(Index s, Index b, Index p, Index r, Index c) {
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  if (msr.mextent != s || mbr.mextent != b || mpr.mextent != p ||
      mrr.mextent != r || mcr.mextent != c) {
    delete[] mdata;
    mdata = new Numeric[s * b * p * r * c];

    msr.mstart = 0;
    msr.mextent = s;
    msr.mstride = b * p * r * c;

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
void swap(Tensor5& t1, Tensor5& t2) {
  std::swap(t1.msr, t2.msr);
  std::swap(t1.mbr, t2.mbr);
  std::swap(t1.mpr, t2.mpr);
  std::swap(t1.mrr, t2.mrr);
  std::swap(t1.mcr, t2.mcr);
  std::swap(t1.mdata, t2.mdata);
}

/** Destructor for Tensor5. This is important, since Tensor5 uses new to
    allocate storage. */
Tensor5::~Tensor5() {
  //   cout << "Destroying a Tensor5:\n"
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
void transform(Tensor5View y, double (&my_func)(double), ConstTensor5View x) {
  // Check dimensions:
  ARTS_ASSERT(y.nshelves() == x.nshelves());
  ARTS_ASSERT(y.nbooks() == x.nbooks());
  ARTS_ASSERT(y.npages() == x.npages());
  ARTS_ASSERT(y.nrows() == x.nrows());
  ARTS_ASSERT(y.ncols() == x.ncols());

  const ConstIterator5D xe = x.end();
  ConstIterator5D xi = x.begin();
  Iterator5D yi = y.begin();
  for (; xi != xe; ++xi, ++yi) {
    // Use the transform function of lower dimensional tensors
    // recursively:
    transform(*yi, my_func, *xi);
  }
}

/** Max function, tensor version. */
Numeric max(const ConstTensor5View& x) {
  const ConstIterator5D xe = x.end();
  ConstIterator5D xi = x.begin();

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
Numeric min(const ConstTensor5View& x) {
  const ConstIterator5D xe = x.end();
  ConstIterator5D xi = x.begin();

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
    \param s  Shelf index
    \param b  Book index
    \param p  Page index
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Numeric debug_tensor5view_get_elem(
    Tensor5View& tv, Index s, Index b, Index p, Index r, Index c) {
  return tv(s, b, p, r, c);
}

#endif
////////////////////////////////
