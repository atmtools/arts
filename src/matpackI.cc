/* Copyright (C) 2001-2019 Stefan Buehler <sbuehler@ltu.se>

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
   \file   matpackI.cc
   \author Stefan Buehler
   \date   2001-09-15
*/

#include "matpackI.h"
#include <cmath>
#include <cstring>
#include "blas.h"
#include "exceptions.h"

using std::cout;
using std::endl;
using std::runtime_error;
using std::setw;

// Define the global joker object:
extern const Joker joker = Joker();

static const Numeric RAD2DEG = 57.295779513082323;

std::ostream& operator<<(std::ostream& os, const Range& r) {
  os << "Range(" << r.get_start() << ", " << r.get_extent() << ", "
     << r.get_stride() << ")";
  return os;
}

// Functions for ConstVectorView:
// ------------------------------

bool ConstVectorView::empty() const ARTS_NOEXCEPT { return (nelem() == 0); }

Index ConstVectorView::nelem() const ARTS_NOEXCEPT { return mrange.mextent; }

Index ConstVectorView::size() const ARTS_NOEXCEPT { return mrange.mextent; }

Numeric ConstVectorView::sum() const ARTS_NOEXCEPT {
  Numeric s = 0;
  ConstIterator1D i = begin();
  const ConstIterator1D e = end();

  for (; i != e; ++i) s += *i;

  return s;
}

ConstVectorView ConstVectorView::operator[](const Range& r) const ARTS_NOEXCEPT {
  return ConstVectorView(mdata, mrange, r);
}

ConstIterator1D ConstVectorView::begin() const ARTS_NOEXCEPT {
  return ConstIterator1D(mdata + mrange.mstart, mrange.mstride);
}

ConstIterator1D ConstVectorView::end() const ARTS_NOEXCEPT {
  return ConstIterator1D(
      mdata + mrange.mstart + (mrange.mextent) * mrange.mstride,
      mrange.mstride);
}

ConstVectorView::operator ConstMatrixView() const {
  // The old version (before 2013-01-18) of this was:
  //    return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
  // Bus this was a bug! The problem is that the matrix index operator adds
  // the mstart from both row and columm range object to mdata

  return ConstMatrixView(mdata, mrange, Range(0, 1));
}

/* This one is a bit tricky: We have to cast away the arguments const
   qualifier, because mdata is not const. This should be safe, since
   there are no non-const methods for ConstVectorView. */
ConstVectorView::ConstVectorView(const Numeric& a) ARTS_NOEXCEPT
    : mrange(0, 1), mdata(&const_cast<Numeric&>(a)) {
  // Nothing to do here.
}

ConstVectorView::ConstVectorView(Numeric* data, const Range& range) ARTS_NOEXCEPT
    : mrange(range), mdata(data) {
  // Nothing to do here.
}

ConstVectorView::ConstVectorView(Numeric* data, const Range& p, const Range& n) ARTS_NOEXCEPT
    : mrange(p, n), mdata(data) {
  // Nothing to do here.
}

/* Output operator. This demonstrates how iterators can be used to
   traverse the vector. The iterators know which part of the vector
   is `active', and also the stride. */
std::ostream& operator<<(std::ostream& os, const ConstVectorView& v) {
  ConstIterator1D i = v.begin();
  const ConstIterator1D end = v.end();

  if (i != end) {
    os << setw(3) << *i;
    ++i;
  }
  for (; i != end; ++i) {
    os << " " << setw(3) << *i;
  }

  return os;
}

// Functions for VectorView:
// ------------------------

VectorView::VectorView(const Vector&) {
  ARTS_ASSERT (false,
      "Creating a VectorView from a const Vector is not allowed.\n"
      "This is not really a runtime error, but I don't want to start\n"
      "producing direct output from inside matpack. And just exiting is\n"
      "not so nice.\n"
      "If you see this error, there is a bug in the code, not in the\n"
      "ARTS input.")
}

VectorView::VectorView(Vector& v) ARTS_NOEXCEPT {
  mdata = v.mdata;
  mrange = v.mrange;
}

VectorView VectorView::operator[](const Range& r) ARTS_NOEXCEPT {
  return VectorView(mdata, mrange, r);
}

Iterator1D VectorView::begin() ARTS_NOEXCEPT {
  return Iterator1D(mdata + mrange.mstart, mrange.mstride);
}

Iterator1D VectorView::end() ARTS_NOEXCEPT {
  return Iterator1D(mdata + mrange.mstart + (mrange.mextent) * mrange.mstride,
                    mrange.mstride);
}

VectorView& VectorView::operator=(const ConstVectorView& v) {
  //  cout << "Assigning VectorView from ConstVectorView.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.mrange.mextent);

  copy(v.begin(), v.end(), begin());

  return *this;
}

VectorView& VectorView::operator=(const VectorView& v) {
  //  cout << "Assigning VectorView from VectorView.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.mrange.mextent);

  copy(v.begin(), v.end(), begin());

  return *this;
}

VectorView& VectorView::operator=(const Vector& v) {
  //  cout << "Assigning VectorView from Vector.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.mrange.mextent);

  copy(v.begin(), v.end(), begin());

  return *this;
}

VectorView& VectorView::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

VectorView VectorView::operator*=(Numeric x) ARTS_NOEXCEPT {
  const Iterator1D e = end();
  for (Iterator1D i = begin(); i != e; ++i) *i *= x;
  return *this;
}

VectorView VectorView::operator/=(Numeric x) ARTS_NOEXCEPT {
  const Iterator1D e = end();
  for (Iterator1D i = begin(); i != e; ++i) *i /= x;
  return *this;
}

VectorView VectorView::operator+=(Numeric x) ARTS_NOEXCEPT {
  const Iterator1D e = end();
  for (Iterator1D i = begin(); i != e; ++i) *i += x;
  return *this;
}

VectorView VectorView::operator-=(Numeric x) ARTS_NOEXCEPT {
  const Iterator1D e = end();
  for (Iterator1D i = begin(); i != e; ++i) *i -= x;
  return *this;
}

VectorView VectorView::operator*=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  Iterator1D i = begin();
  const Iterator1D e = end();

  for (; i != e; ++i, ++s) *i *= *s;
  return *this;
}

VectorView VectorView::operator/=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  Iterator1D i = begin();
  const Iterator1D e = end();

  for (; i != e; ++i, ++s) *i /= *s;
  return *this;
}

VectorView VectorView::operator+=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  Iterator1D i = begin();
  const Iterator1D e = end();

  for (; i != e; ++i, ++s) *i += *s;
  return *this;
}

VectorView VectorView::operator-=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  Iterator1D i = begin();
  const Iterator1D e = end();

  for (; i != e; ++i, ++s) *i -= *s;
  return *this;
}

VectorView::operator MatrixView() ARTS_NOEXCEPT {
  // The old version (before 2013-01-18) of this was:
  //    return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
  // Bus this was a bug! The problem is that the matrix index operator adds
  // the mstart from both row and columm range object to mdata

  return MatrixView(mdata, mrange, Range(0, 1));
}

const Numeric* VectorView::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT(not (mrange.mstart != 0 || mrange.mstride != 1),
        "A VectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

Numeric* VectorView::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT(not (mrange.mstart != 0 || mrange.mstride != 1),
        "A VectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

VectorView::VectorView(Numeric& a) ARTS_NOEXCEPT : ConstVectorView(a) {
  // Nothing to do here.
}

VectorView::VectorView(Numeric* data, const Range& range) ARTS_NOEXCEPT
    : ConstVectorView(data, range) {
  // Nothing to do here.
}

VectorView::VectorView(Numeric* data, const Range& p, const Range& n) ARTS_NOEXCEPT
    : ConstVectorView(data, p, n) {
  // Nothing to do here.
}

void copy(ConstIterator1D origin,
          const ConstIterator1D& end,
          Iterator1D target) {
  if (origin.mstride == 1 && target.mstride == 1)
    memcpy((void*)target.mx,
           (void*)origin.mx,
           sizeof(Numeric) * (end.mx - origin.mx));
  else
    for (; origin != end; ++origin, ++target) *target = *origin;
}

void copy(Numeric x, Iterator1D target, const Iterator1D& end) ARTS_NOEXCEPT {
  for (; target != end; ++target) *target = x;
}

// Functions for Vector:
// ---------------------

Vector::Vector(std::initializer_list<Numeric> init)
    : VectorView(new Numeric[init.size()], Range(0, init.size())) {
  std::copy(init.begin(), init.end(), begin());
}

Vector::Vector(Index n) : VectorView(new Numeric[n], Range(0, n)) {
  // Nothing to do here.
}

Vector::Vector(Index n, Numeric fill)
    : VectorView(new Numeric[n], Range(0, n)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, n, fill);
}

Vector::Vector(Numeric start, Index extent, Numeric stride)
    : VectorView(new Numeric[extent], Range(0, extent)) {
  // Fill with values:
  Numeric x = start;
  Iterator1D i = begin();
  const Iterator1D e = end();
  for (; i != e; ++i) {
    *i = x;
    x += stride;
  }
}

Vector::Vector(const ConstVectorView& v)
    : VectorView(new Numeric[v.nelem()], Range(0, v.nelem())) {
  copy(v.begin(), v.end(), begin());
}

Vector::Vector(const Vector& v)
    : VectorView(new Numeric[v.nelem()], Range(0, v.nelem())) {
  std::memcpy(mdata, v.mdata, nelem() * sizeof(Numeric));
}

Vector::Vector(const std::vector<Numeric>& v)
    : VectorView(new Numeric[v.size()], Range(0, v.size())) {
  std::vector<Numeric>::const_iterator vec_it_end = v.end();
  Iterator1D this_it = this->begin();
  for (std::vector<Numeric>::const_iterator vec_it = v.begin();
       vec_it != vec_it_end;
       ++vec_it, ++this_it)
    *this_it = *vec_it;
}

Vector& Vector::operator=(std::initializer_list<Numeric> v) {
  resize(v.size());
  std::copy(v.begin(), v.end(), begin());
  return *this;
}

Vector& Vector::operator=(const Vector& v) {
  if (this != &v) {
    resize(v.nelem());
    std::memcpy(mdata, v.mdata, nelem() * sizeof(Numeric));
  }
  return *this;
}

Vector& Vector::operator=(Vector&& v) noexcept {
  if (this != &v) {
    delete[] mdata;
    mdata = v.mdata;
    mrange = v.mrange;
    v.mrange = Range(0, 0);
    v.mdata = nullptr;
  }
  return *this;
}

Vector& Vector::operator=(const Array<Numeric>& x) {
  resize(x.nelem());
  VectorView::operator=(x);
  return *this;
}

Vector& Vector::operator=(Numeric x) {
  std::fill_n(mdata, nelem(), x);
  return *this;
}

void Vector::resize(Index n) {
  ARTS_ASSERT(0 <= n);
  if (mrange.mextent != n) {
    delete[] mdata;
    mdata = new Numeric[n];
    mrange.mstart = 0;
    mrange.mextent = n;
    mrange.mstride = 1;
  }
}

void swap(Vector& v1, Vector& v2) {
  std::swap(v1.mrange, v2.mrange);
  std::swap(v1.mdata, v2.mdata);
}

Vector::~Vector() { delete[] mdata; }

// Functions for ConstMatrixView:
// ------------------------------

//! Returns true if variable size is zero.
bool ConstMatrixView::empty() const ARTS_NOEXCEPT { return (nrows() == 0 || ncols() == 0); }

/** Returns the number of rows. */
Index ConstMatrixView::nrows() const ARTS_NOEXCEPT { return mrr.mextent; }

/** Returns the number of columns. */
Index ConstMatrixView::ncols() const ARTS_NOEXCEPT { return mcr.mextent; }

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior.  */
ConstMatrixView ConstMatrixView::operator()(const Range& r,
                                            const Range& c) const ARTS_NOEXCEPT {
  return ConstMatrixView(mdata, mrr, mcr, r, c);
}

/** Const index operator returning a column as an object of type
    ConstVectorView.

    \param r A range of rows.
    \param c Index of selected column */
ConstVectorView ConstMatrixView::operator()(const Range& r, Index c) const ARTS_NOEXCEPT {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(mdata + mcr.mstart + c * mcr.mstride, mrr, r);
}

/** Const index operator returning a row as an object of type
    ConstVectorView.

    \param r Index of selected row.
    \param c Range of columns */
ConstVectorView ConstMatrixView::operator()(Index r, const Range& c) const ARTS_NOEXCEPT {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstVectorView(mdata + mrr.mstart + r * mrr.mstride, mcr, c);
}

/** Return const iterator to first row. */
ConstIterator2D ConstMatrixView::begin() const ARTS_NOEXCEPT {
  return ConstIterator2D(ConstVectorView(mdata + mrr.mstart, mcr), mrr.mstride);
}

/** Return const iterator behind last row. */
ConstIterator2D ConstMatrixView::end() const ARTS_NOEXCEPT {
  return ConstIterator2D(
      ConstVectorView(mdata + mrr.mstart + (mrr.mextent) * mrr.mstride, mcr),
      mrr.mstride);
}

//! Matrix diagonal as vector.
/*!
  Returns a ConstMatrixView on the diagonal entries of the matrix. For a given
  (n,m) matrix M the diagonal vector v is the vector of length min{n,m} with entries

       v[i] = M(i,i)

  \return The diagonal vector v.
*/
ConstVectorView ConstMatrixView::diagonal() const ARTS_NOEXCEPT {
  Index n = std::min(mrr.mextent, mcr.mextent);
  return ConstVectorView(mdata + mrr.mstart + mcr.mstart,
                         Range(0, n, mrr.mstride + mcr.mstride));
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
ConstMatrixView::ConstMatrixView(Numeric* data,
                                 const Range& rr,
                                 const Range& cr) ARTS_NOEXCEPT
    : mrr(rr), mcr(cr), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pr Previous range.
    \param pc Previous range.
    \param nr New Range.
    \param nc New Range.
  */
ConstMatrixView::ConstMatrixView(Numeric* data,
                                 const Range& pr,
                                 const Range& pc,
                                 const Range& nr,
                                 const Range& nc) ARTS_NOEXCEPT
    : mrr(pr, nr), mcr(pc, nc), mdata(data) {
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the matrix. The iterators know which part of the matrix
    is `active', and also the strides in both directions. This
    function is a bit more complicated than necessary to illustrate
    the concept, because the formating should look nice. This means
    that the first row, and the first element in each row, have to be
    treated individually. */
std::ostream& operator<<(std::ostream& os, const ConstMatrixView& v) {
  // Row iterators:
  ConstIterator2D ir = v.begin();
  const ConstIterator2D end_row = v.end();

  if (ir != end_row) {
    ConstIterator1D ic = ir->begin();
    const ConstIterator1D end_col = ir->end();

    if (ic != end_col) {
      os << setw(3) << *ic;
      ++ic;
    }
    for (; ic != end_col; ++ic) {
      os << " " << setw(3) << *ic;
    }
    ++ir;
  }
  for (; ir != end_row; ++ir) {
    ConstIterator1D ic = ir->begin();
    const ConstIterator1D end_col = ir->end();

    os << "\n";
    if (ic != end_col) {
      os << setw(3) << *ic;
      ++ic;
    }
    for (; ic != end_col; ++ic) {
      os << " " << setw(3) << *ic;
    }
  }

  return os;
}

// Functions for MatrixView:
// -------------------------

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Matrix. This allows correct
    recursive behavior.  */
MatrixView MatrixView::operator()(const Range& r, const Range& c) ARTS_NOEXCEPT {
  return MatrixView(mdata, mrr, mcr, r, c);
}

/** Index operator returning a column as an object of type VectorView.

    \param r A range of rows.
    \param c Index of selected column */
VectorView MatrixView::operator()(const Range& r, Index c) ARTS_NOEXCEPT {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(mdata + mcr.mstart + c * mcr.mstride, mrr, r);
}

/** Index operator returning a row as an object of type VectorView.

    \param r Index of selected row.
    \param c Range of columns */
VectorView MatrixView::operator()(Index r, const Range& c) ARTS_NOEXCEPT {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return VectorView(mdata + mrr.mstart + r * mrr.mstride, mcr, c);
}

///** Return const iterator to first row. Has to be redefined here, since it is
//    hiden by the non-const operator of the derived class.*/
//ConstIterator2D MatrixView::begin() const
//{
//  return ConstMatrixView::begin();
//}
//
///** Return const iterator behind last row. */
//ConstIterator2D MatrixView::end() const
//{
//  return ConstMatrixView::end();
//}

/** Return iterator to first row. */
Iterator2D MatrixView::begin() ARTS_NOEXCEPT {
  return Iterator2D(VectorView(mdata + mrr.mstart, mcr), mrr.mstride);
}

/** Return iterator behind last row. */
Iterator2D MatrixView::end() ARTS_NOEXCEPT {
  return Iterator2D(
      VectorView(mdata + mrr.mstart + (mrr.mextent) * mrr.mstride, mcr),
      mrr.mstride);
}

/** Assignment operator. This copies the data from another MatrixView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
MatrixView& MatrixView::operator=(const ConstMatrixView& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from MatrixView to MatrixView. This is a tricky
    one. The problem is that since MatrixView is derived from
    ConstMatrixView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
MatrixView& MatrixView::operator=(const MatrixView& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a Matrix. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
MatrixView& MatrixView::operator=(const Matrix& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a vector. This copies the data from a VectorView
    to this MatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this MatrixView by
    setting its range. */
MatrixView& MatrixView::operator=(const ConstVectorView& v) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == v.nelem());
  ARTS_ASSERT(mcr.mextent == 1);
  //  dummy = ConstMatrixView(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstMatrixView dummy(v);
  copy(dummy.begin(), dummy.end(), begin());
  return *this;
}

/** Assigning a scalar to a MatrixView will set all elements to this
    value. */
MatrixView& MatrixView::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

/** Multiplication by scalar. */
MatrixView& MatrixView::operator*=(Numeric x) ARTS_NOEXCEPT {
  const Iterator2D er = end();
  for (Iterator2D r = begin(); r != er; ++r) {
    const Iterator1D ec = r->end();
    for (Iterator1D c = r->begin(); c != ec; ++c) *c *= x;
  }
  return *this;
}

/** Division by scalar. */
MatrixView& MatrixView::operator/=(Numeric x) ARTS_NOEXCEPT {
  const Iterator2D er = end();
  for (Iterator2D r = begin(); r != er; ++r) {
    const Iterator1D ec = r->end();
    for (Iterator1D c = r->begin(); c != ec; ++c) *c /= x;
  }
  return *this;
}

/** Addition of scalar. */
MatrixView& MatrixView::operator+=(Numeric x) ARTS_NOEXCEPT {
  const Iterator2D er = end();
  for (Iterator2D r = begin(); r != er; ++r) {
    const Iterator1D ec = r->end();
    for (Iterator1D c = r->begin(); c != ec; ++c) *c += x;
  }
  return *this;
}

/** Subtraction of scalar. */
MatrixView& MatrixView::operator-=(Numeric x) ARTS_NOEXCEPT {
  const Iterator2D er = end();
  for (Iterator2D r = begin(); r != er; ++r) {
    const Iterator1D ec = r->end();
    for (Iterator1D c = r->begin(); c != ec; ++c) *c -= x;
  }
  return *this;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric* MatrixView::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT(not (mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 || mcr.mstride != 1),
    "A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric* MatrixView::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT(not (mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 || mcr.mstride != 1),
    "A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** Element-vise multiplication by another Matrix. */
MatrixView& MatrixView::operator*=(const ConstMatrixView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    Iterator1D c = r->begin();
    const Iterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c *= *sc;
  }
  return *this;
}

/** Element-vise division by another Matrix. */
MatrixView& MatrixView::operator/=(const ConstMatrixView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    Iterator1D c = r->begin();
    const Iterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c /= *sc;
  }
  return *this;
}

/** Element-vise addition of another Matrix. */
MatrixView& MatrixView::operator+=(const ConstMatrixView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    Iterator1D c = r->begin();
    const Iterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c += *sc;
  }
  return *this;
}

/** Element-vise subtraction of another Matrix. */
MatrixView& MatrixView::operator-=(const ConstMatrixView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    Iterator1D c = r->begin();
    const Iterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c -= *sc;
  }
  return *this;
}

/** Element-vise multiplication by a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator*=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nelem());
  ARTS_ASSERT(ncols() == 1);
  ConstIterator1D sc = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sc) {
    Iterator1D c = r->begin();
    *c *= *sc;
  }
  return *this;
}

/** Element-vise division by a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator/=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nelem());
  ARTS_ASSERT(ncols() == 1);
  ConstIterator1D sc = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sc) {
    Iterator1D c = r->begin();
    *c /= *sc;
  }
  return *this;
}

/** Element-vise addition of a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator+=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nelem());
  ARTS_ASSERT(ncols() == 1);
  ConstIterator1D sc = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sc) {
    Iterator1D c = r->begin();
    *c += *sc;
  }
  return *this;
}

/** Element-vise subtraction of a Vector (acting like a 1-column Matrix). */
MatrixView& MatrixView::operator-=(const ConstVectorView& x) ARTS_NOEXCEPT {
  ARTS_ASSERT(nrows() == x.nelem());
  ARTS_ASSERT(ncols() == 1);
  ConstIterator1D sc = x.begin();
  Iterator2D r = begin();
  const Iterator2D er = end();
  for (; r != er; ++r, ++sc) {
    Iterator1D c = r->begin();
    *c -= *sc;
  }
  return *this;
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
MatrixView::MatrixView(Numeric* data, const Range& rr, const Range& cr) ARTS_NOEXCEPT
    : ConstMatrixView(data, rr, cr) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pr Previous range.
    \param pc Previous range.
    \param nr New Range.
    \param nc New Range.
  */
MatrixView::MatrixView(Numeric* data,
                       const Range& pr,
                       const Range& pc,
                       const Range& nr,
                       const Range& nc) ARTS_NOEXCEPT
    : ConstMatrixView(data, pr, pc, nr, nc) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors.

    Origin, end, and target are 2D iterators, marking rows in a
    matrix. For each row the 1D iterator is obtained and used to copy
    the elements. 
*/
void copy(ConstIterator2D origin,
          const ConstIterator2D& end,
          Iterator2D target) {
  for (; origin != end; ++origin, ++target) {
    ConstIterator1D o = origin->begin();
    const ConstIterator1D e = origin->end();
    Iterator1D t = target->begin();
    for (; o != e; ++o, ++t) *t = *o;
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, Iterator2D target, const Iterator2D& end) ARTS_NOEXCEPT {
  for (; target != end; ++target) {
    Iterator1D t = target->begin();
    const Iterator1D e = target->end();
    for (; t != e; ++t) *t = x;
  }
}

// Functions for Matrix:
// ---------------------

/** Constructor setting size. This constructor has to set the stride
    in the row range correctly! */
Matrix::Matrix(Index r, Index c)
    : MatrixView(new Numeric[r * c], Range(0, r, c), Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Matrix::Matrix(Index r, Index c, Numeric fill)
    : MatrixView(new Numeric[r * c], Range(0, r, c), Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, r * c, fill);
}

/** Copy constructor from MatrixView. This automatically sets the size
    and copies the data. */
Matrix::Matrix(const ConstMatrixView& m)
    : MatrixView(new Numeric[m.nrows() * m.ncols()],
                 Range(0, m.nrows(), m.ncols()),
                 Range(0, m.ncols())) {
  copy(m.begin(), m.end(), begin());
}

/** Copy constructor from Matrix. This automatically sets the size
    and copies the data. */
Matrix::Matrix(const Matrix& m)
    : MatrixView(new Numeric[m.nrows() * m.ncols()],
                 Range(0, m.nrows(), m.ncols()),
                 Range(0, m.ncols())) {
  // There is a catch here: If m is an empty matrix, then it will have
  // 0 colunns. But this is used to initialize the stride of the row
  // Range! Thus, this method has to be consistent with the behaviour
  // of Range::Range. For now, Range::Range allows also stride 0.
  std::memcpy(mdata, m.mdata, nrows() * ncols() * sizeof(Numeric));
}

//! Assignment operator from another matrix.
/*! 
  While dimensions of MatrixViews can not be adjusted, dimensions of
  matrices *can* be adjusted. Hence, the behavior of the assignment
  operator is different.

  In this case the size of the target is automatically adjusted. This
  is important, so that structures containing matrices are copied
  correctly. 
  
  This is a deviation from the old ARTS paradigm that sizes must match
  exactly before copying!

  Note: It is sufficient to have only this one version of the
  assignment (Matrix = Matrix). It implicitly covers the cases
  Matrix=MatrixView, etc, because there is a default constructor for
  Matrix from MatrixView. (See C++ Primer Plus, page 571ff.)

  \param m The other matrix to assign to this one.

  \return This matrix, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
Matrix& Matrix::operator=(const Matrix& m) {
  if (this != &m) {
    resize(m.nrows(), m.ncols());
    std::memcpy(mdata, m.mdata, nrows() * ncols() * sizeof(Numeric));
  }
  return *this;
}

//! Move assignment operator from another matrix.
Matrix& Matrix::operator=(Matrix&& m) noexcept {
  if (this != &m) {
    delete[] mdata;
    mdata = m.mdata;
    mrr = m.mrr;
    mcr = m.mcr;
    m.mrr = Range(0, 0);
    m.mcr = Range(0, 0);
    m.mdata = nullptr;
  }
  return *this;
}

/** Assignment operator from scalar. Assignment operators also seem to
    be not inherited. */
Matrix& Matrix::operator=(Numeric x) {
  std::fill_n(mdata, nrows() * ncols(), x);
  return *this;
}

//! Assignment from a vector.
/*! 
  This copies the data from a VectorView to this MatrixView.

  The dimension is adjusted automatically.

  \param v The vector to assign to this matrix.

  \return This matrix, by tradition.

  \author Stefan Buehler
  \date   2002-12-19
*/
Matrix& Matrix::operator=(const ConstVectorView& v) {
  resize(v.nelem(), 1);
  ConstMatrixView dummy(v);
  copy(dummy.begin(), dummy.end(), begin());
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Matrix is not
    initialized, so it will contain random values.*/
void Matrix::resize(Index r, Index c) {
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  if (mrr.mextent != r || mcr.mextent != c) {
    delete[] mdata;
    mdata = new Numeric[r * c];

    mrr.mstart = 0;
    mrr.mextent = r;
    mrr.mstride = c;

    mcr.mstart = 0;
    mcr.mextent = c;
    mcr.mstride = 1;
  }
}

/** Swaps two objects. */
void swap(Matrix& m1, Matrix& m2) {
  std::swap(m1.mrr, m2.mrr);
  std::swap(m1.mcr, m2.mcr);
  std::swap(m1.mdata, m2.mdata);
}

/** Destructor for Matrix. This is important, since Matrix uses new to
    allocate storage. */
Matrix::~Matrix() {
  //   cout << "Destroying a Matrix:\n"
  //        << *this << "\n........................................\n";
  delete[] mdata;
}

// Some general Matrix Vector functions:

/** Scalar product. The two vectors may be identical. */
Numeric operator*(const ConstVectorView& a, const ConstVectorView& b) ARTS_NOEXCEPT {
  // Check dimensions:
  ARTS_ASSERT(a.nelem() == b.nelem());

  const ConstIterator1D ae = a.end();
  ConstIterator1D ai = a.begin();
  ConstIterator1D bi = b.begin();

  Numeric res = 0;
  for (; ai != ae; ++ai, ++bi) res += (*ai) * (*bi);

  return res;
}

//! Matrix-Vector Multiplication
/*!

  Computes the Matrix-Vector product y = M * x, for a m times n matrix M, a
  length-m vector y and a length-n vector x.

  The product is computed using the dgemv_ routine from the BLAS library if
  the matrix is contiguous in memory. If this is not the case, the mult_general
  method is used to compute the product.

  No memory is allocated for the computation and the matrix and vector views
  may not overlap.

  \param[out] y The length-m VectorView where the result is stored.
  \param[in] M Reference to the m-times-n ConstMatrixView holding the matrix M.
  \param[in] x Reference to the length-n ConstVectorView holding the vector x.
*/
void mult(VectorView y, const ConstMatrixView& M, const ConstVectorView& x) {
  ARTS_ASSERT(y.mrange.get_extent() == M.mrr.get_extent());
  ARTS_ASSERT(M.mcr.get_extent() == x.mrange.get_extent());
  ARTS_ASSERT((M.mcr.get_extent() != 0) && (M.mrr.get_extent() != 0));

  if ((M.mcr.get_stride() == 1) || (M.mrr.get_stride() == 1)) {
    char trans;
    int m, n;
    double zero = 0.0;
    double one = 1.0;
    int LDA, incx, incy;

    if (M.mcr.get_stride() != 1) {
      trans = 'n';
      m = (int)M.mrr.get_extent();
      n = (int)M.mcr.get_extent();
      LDA = (int)M.mcr.get_stride();
    } else {
      trans = 't';
      m = (int)M.mcr.get_extent();
      n = (int)M.mrr.get_extent();
      LDA = (int)M.mrr.get_stride();
      if (M.mrr.get_stride() == 1) LDA = m;
    }

    incx = (int)x.mrange.get_stride();
    incy = (int)y.mrange.get_stride();

    double* mstart = M.mdata + M.mcr.get_start() + M.mrr.get_start();
    double* ystart = y.mdata + y.mrange.get_start();
    double* xstart = x.mdata + x.mrange.get_start();

    dgemv_(&trans,
           &m,
           &n,
           &one,
           mstart,
           &LDA,
           xstart,
           &incx,
           &zero,
           ystart,
           &incy);

  } else {
    mult_general(y, M, x);
  }
}

/** Matrix Vector multiplication. y = M*x. Note that the order is different
    from MTL, output comes first! Dimensions of y, M, and x must
    match. No memory reallocation takes place, only the data is
    copied. Using this function on overlapping Matrix and VectorViews belonging
    to the same Matrix will lead to unpredictable results.

    The implementation here is different from the other multiplication
    routines. It does not use iterators but a more drastic approach to gain
    maximum performance.  */
void mult_general(VectorView y,
                  const ConstMatrixView& M,
                  const ConstVectorView& x) ARTS_NOEXCEPT {
  // Check dimensions:
  ARTS_ASSERT(y.mrange.mextent == M.mrr.mextent);
  ARTS_ASSERT(M.mcr.mextent == x.mrange.mextent);
  ARTS_ASSERT(M.mcr.mextent != 0 && M.mrr.mextent != 0);

  // Let's first find the pointers to the starting positions
  Numeric* mdata = M.mdata + M.mcr.mstart + M.mrr.mstart;
  Numeric* xdata = x.mdata + x.mrange.mstart;
  Numeric* yelem = y.mdata + y.mrange.mstart;

  Index i = M.mrr.mextent;
  while (i--) {
    Numeric* melem = mdata;
    Numeric* xelem = xdata;  // Reset xelem to first element of source vector

    // Multiply first element of matrix row with first element of source
    // vector. This is done outside the loop to avoid initialization of the
    // target vector's element with zero (saves one assignment)
    *yelem = *melem * *xelem;

    Index j = M.mcr.mextent;  // --j (instead of j-- like in the outer loop)
    while (--j)               // is important here because we only want
    {                         // mextent-1 cycles
      melem += M.mcr.mstride;
      xelem += x.mrange.mstride;
      *yelem += *melem * *xelem;
    }

    mdata += M.mrr.mstride;     // Jump to next matrix row
    yelem += y.mrange.mstride;  // Jump to next element in target vector
  }
}

//! Matrix-Matrix Multiplication
/*!
  Performs the matrix multiplication A = B * C. The dimensions must match, i.e.
  A must be a m times n matrix, B a m times k matrix and C a k times c matrix.
  No memory reallocation takes place, only the data is copied. Using this function
  on overlapping MatrixViews belonging to the same Matrix will lead to unpredictable
  results. In particular, this means that A and B must not be the same matrix!

  If the memory layout allows it, the multiplication is performed using BLAS'
  _dgemm, which leads to a significant speed up of the operation. To be compatible
  with BLAS the matrix views A, B and C must satisfy:

  - A must have column or row stride 1
  - B must have column or row stride 1
  - C must have column stride 1

  That means that A and B can be ConstMatrixView objects corresponding to
  transposed/non-transposed submatrices of a matrix, that are continuous along their
  first/second dimension. C must correspond to a non-transposed submatrix of a
  matrix, that is continuous along its second dimension.

  \param[in,out] A The matrix A, that will hold the result of the multiplication.
  \param[in] B The matrix B
  \param[in] C The matrix C
*/
void mult(MatrixView A, const ConstMatrixView& B, const ConstMatrixView& C) {
  // Check dimensions:
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  // Catch trivial case if one of the matrices is empty.
  if ((B.nrows() == 0) || (B.ncols() == 0) || (C.ncols() == 0)) return;

  // Matrices B and C must be continuous in at least on dimension,  C
  // must be continuous along the second dimension.
  if (((B.mrr.get_stride() == 1) || (B.mcr.get_stride() == 1)) &&
      ((C.mrr.get_stride() == 1) || (C.mcr.get_stride() == 1)) &&
      (A.mcr.get_stride() == 1)) {
    // BLAS uses column-major order while arts uses row-major order.
    // Hence instead of C = A * B we compute C^T = A^T * B^T!

    int k, m, n;

    k = (int)B.ncols();
    m = (int)C.ncols();
    n = (int)B.nrows();

    // Note also the clash in nomenclature: BLAS uses C = A * B while
    // arts uses A = B * C. Taking into accout this and the difference in
    // memory layouts, we need to map the MatrixViews A, B and C to the BLAS
    // arguments as follows:
    // A (arts) -> C (BLAS)
    // B (arts) -> B (BLAS)
    // C (arts) -> A (BLAS)

    // Char indicating whether A (BLAS) or B (BLAS) should be transposed.
    char transa, transb;
    // Sizes of the matrices along the direction in which they are
    // traversed.
    int lda, ldb, ldc;

    // Check if C (arts) is transposed.
    if (C.mrr.get_stride() == 1) {
      transa = 'T';
      lda = (int)C.mcr.get_stride();
    } else {
      transa = 'N';
      lda = (int)C.mrr.get_stride();
    }

    // Check if B (arts) is transposed.
    if (B.mrr.get_stride() == 1) {
      transb = 'T';
      ldb = (int)B.mcr.get_stride();
    } else {
      transb = 'N';
      ldb = (int)B.mrr.get_stride();
    }

    // In case B (arts) has only one column, column and row stride are 1.
    // We therefore need to set ldb to k, since dgemm_ requires lda to be at
    // least k / m if A is non-transposed / transposed.
    if ((B.mcr.get_stride() == 1) && (B.mrr.get_stride() == 1)) {
      transb = 'N';
      ldb = k;
    }

    // The same holds for C (arts).
    if ((C.mcr.get_stride() == 1) && (C.mrr.get_stride() == 1)) {
      transa = 'N';
      lda = m;
    }

    ldc = (int)A.mrr.get_stride();
    // The same holds for A (arts).
    if ((A.mcr.get_stride() == 1) && (A.mrr.get_stride() == 1)) {
      ldc = m;
    }
    double alpha = 1.0, beta = 0.0;

    dgemm_(&transa,
           &transb,
           &m,
           &n,
           &k,
           &alpha,
           C.mdata + C.mrr.get_start() + C.mcr.get_start(),
           &lda,
           B.mdata + B.mrr.get_start() + B.mcr.get_start(),
           &ldb,
           &beta,
           A.mdata + A.mrr.get_start() + A.mcr.get_start(),
           &ldc);

  } else {
    mult_general(A, B, C);
  }
}

//! General matrix multiplication.
/*!
  This is the fallback matrix multiplication which works for all
  ConstMatrixView objects.

  \param[in,out] A The matrix A, that will hold the result of the multiplication.
  \param[in] B The matrix B
  \param[in] C The matrix C
*/
void mult_general(MatrixView A,
                  const ConstMatrixView& B,
                  const ConstMatrixView& C) ARTS_NOEXCEPT {
  // Check dimensions:
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  // Let's get the transpose of C, so that we can use 2D iterators to
  // access the columns (= rows of the transpose).
  ConstMatrixView CT = transpose(C);

  const Iterator2D ae = A.end();
  Iterator2D ai = A.begin();
  ConstIterator2D bi = B.begin();

  // This walks through the rows of A and B:
  for (; ai != ae; ++ai, ++bi) {
    const Iterator1D ace = ai->end();
    Iterator1D aci = ai->begin();
    ConstIterator2D cti = CT.begin();

    // This walks through the columns of A with a 1D iterator, and
    // at the same time through the rows of CT, which are the columns of
    // C, with a 2D iterator:
    for (; aci != ace; ++aci, ++cti) {
      // The operator * is used to compute the scalar product
      // between rows of B and rows of C.transpose().
      *aci = (*bi) * (*cti);
    }
  }
}

//! cross3
/*!
    Calculates the cross product between two vectors of length 3.

    c = a x b, for 3D vectors. The vector c must have length 3 and can not be
    the same variable as a or b.

    param    c   Out: The cross product vector
    \param   a   In: A vector of length 3.
    \param   b   In: A vector of length 3.

    \author Patrick Eriksson
    \date   2012-02-12
*/
void cross3(VectorView c, const ConstVectorView& a, const ConstVectorView& b) ARTS_NOEXCEPT {
  ARTS_ASSERT(a.nelem() == 3);
  ARTS_ASSERT(b.nelem() == 3);
  ARTS_ASSERT(c.nelem() == 3);

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/*!
    Returns numeric angle between two vectors in degrees.

    \param  a   In:  A vector of length N.
    \param  b   In:  A vector of length N.

    \author Richard Larsson
    \date   2012-07-10
*/
Numeric vector_angle(ConstVectorView a, ConstVectorView b) {
  ARTS_ASSERT(a.nelem() == b.nelem());
  Numeric arg = (a * b) / sqrt(a * a) / sqrt(b * b);

  // Due to numerical inaccuracy, arg might be slightly larger than 1.
  // We catch those cases to avoid spurious returns of NaNs
  return fabs(arg) > 1. ? 0. : acos(arg) * RAD2DEG;
}

/*!
    Calculates the projection of two vectors of equal length.

    c = proj_a(b). Projecting b on a. The vector c must have the same length
    but can not be the same variable as a or b.

    \param   c   Out: The projection of b on a.
    \param   a   In:  A vector of length N.
    \param   b   In:  A vector of length N.

    \author Richard Larsson
    \date   2012-07-10
*/
void proj(Vector& c, ConstVectorView a, ConstVectorView b) ARTS_NOEXCEPT {
  ARTS_ASSERT(a.nelem() == b.nelem());
  ARTS_ASSERT(a.nelem() == c.nelem());

  const Numeric C = (a * b) / (a * a);
  c = a;
  c *= C;
};

/** Const version of transpose. */
ConstMatrixView transpose(ConstMatrixView m) ARTS_NOEXCEPT {
  return ConstMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
    transpose. The original is not changed! */
MatrixView transpose(MatrixView m) ARTS_NOEXCEPT { return MatrixView(m.mdata, m.mcr, m.mrr); }

/** A generic transform function for vectors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    Although the matrix version of this can also be used for vectors,
    thanks to the automatic interpretation of a vector as a one column
    matrix, this one is slightly more efficient. However, the
    difference is very small (only a few percent). 

    The two views may be the same one, in which case the
    conversion happens in place. 

    \param   y Output:   The results of the function acting on each element of x.
    \param    my_func A function (e.g., sqrt).
    \param    x   A vector. */
void transform(VectorView y, double (&my_func)(double), ConstVectorView x) {
  // Check dimensions:
  ARTS_ASSERT(y.nelem() == x.nelem());

  const ConstIterator1D xe = x.end();
  ConstIterator1D xi = x.begin();
  Iterator1D yi = y.begin();
  for (; xi != xe; ++xi, ++yi) *yi = my_func(*xi);
}

/** A generic transform function for matrices, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    This function can also be used for Vectors, because there is a
    conversion to MatrixView.

    The two Matrix views may be the same one, in which case the
    conversion happens in place. 

   \param   y Output:   The results of the function acting on each element of x.
   \param    my_func A function (e.g., sqrt).
   \param    x   A matrix. */
void transform(MatrixView y, double (&my_func)(double), ConstMatrixView x) {
  // Check dimensions:
  ARTS_ASSERT(y.nrows() == x.nrows());
  ARTS_ASSERT(y.ncols() == x.ncols());

  const ConstIterator2D rxe = x.end();
  ConstIterator2D rx = x.begin();
  Iterator2D ry = y.begin();
  for (; rx != rxe; ++rx, ++ry) {
    const ConstIterator1D cxe = rx->end();
    ConstIterator1D cx = rx->begin();
    Iterator1D cy = ry->begin();
    for (; cx != cxe; ++cx, ++cy) *cy = my_func(*cx);
  }
}

/** Max function, vector version. */
Numeric max(const ConstVectorView& x) ARTS_NOEXCEPT {
  // Initial value for max:
  Numeric max = x[0];

  const ConstIterator1D xe = x.end();
  ConstIterator1D xi = x.begin();

  for (; xi != xe; ++xi) {
    if (*xi > max) max = *xi;
  }

  return max;
}

/** Max function, matrix version. */
Numeric max(const ConstMatrixView& x) ARTS_NOEXCEPT {
  // Initial value for max:
  Numeric max = x(0, 0);

  const ConstIterator2D rxe = x.end();
  ConstIterator2D rx = x.begin();

  for (; rx != rxe; ++rx) {
    const ConstIterator1D cxe = rx->end();
    ConstIterator1D cx = rx->begin();

    for (; cx != cxe; ++cx)
      if (*cx > max) max = *cx;
  }

  return max;
}

/** Min function, vector version. */
Numeric min(const ConstVectorView& x) ARTS_NOEXCEPT {
  // Initial value for min:
  Numeric min = x[0];

  const ConstIterator1D xe = x.end();
  ConstIterator1D xi = x.begin();

  for (; xi != xe; ++xi) {
    if (*xi < min) min = *xi;
  }

  return min;
}

/** Min function, matrix version. */
Numeric min(const ConstMatrixView& x) ARTS_NOEXCEPT {
  // Initial value for min:
  Numeric min = x(0, 0);

  const ConstIterator2D rxe = x.end();
  ConstIterator2D rx = x.begin();

  for (; rx != rxe; ++rx) {
    const ConstIterator1D cxe = rx->end();
    ConstIterator1D cx = rx->begin();

    for (; cx != cxe; ++cx)
      if (*cx < min) min = *cx;
  }

  return min;
}

/** Mean function, vector version. */
Numeric mean(const ConstVectorView& x) ARTS_NOEXCEPT {
  // Initial value for mean:
  Numeric mean = 0;

  const ConstIterator1D xe = x.end();
  ConstIterator1D xi = x.begin();

  for (; xi != xe; ++xi) mean += *xi;

  mean /= (Numeric)x.nelem();

  return mean;
}

/** Mean function, matrix version. */
Numeric mean(const ConstMatrixView& x) ARTS_NOEXCEPT {
  // Initial value for mean:
  Numeric mean = 0;

  const ConstIterator2D rxe = x.end();
  ConstIterator2D rx = x.begin();

  for (; rx != rxe; ++rx) {
    const ConstIterator1D cxe = rx->end();
    ConstIterator1D cx = rx->begin();

    for (; cx != cxe; ++cx) mean += *cx;
  }

  mean /= (Numeric)(x.nrows() * x.ncols());

  return mean;
}

/** Assignment operator from Array<Numeric>. This copies the data from
    an Array<Numeric> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Numeric> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy. */
VectorView& VectorView::operator=(const Array<Numeric>& v) {
  //  cout << "Assigning VectorView from Array<Numeric>.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.nelem());

  // Iterators for Array:
  Array<Numeric>::const_iterator i = v.begin();
  const Array<Numeric>::const_iterator e = v.end();
  // Iterator for Vector:
  Iterator1D target = begin();

  for (; i != e; ++i, ++target) *target = *i;

  return *this;
}

// Const

// Converts constant matrix to constant eigen map
ConstMatrixViewMap MapToEigen(const ConstMatrixView& A) {
  return ConstMatrixViewMap(A.mdata + A.mrr.get_start() + A.mcr.get_start(),
                            A.nrows(),
                            A.ncols(),
                            StrideType(A.mrr.get_stride(), A.mcr.get_stride()));
}

// Converts constant vector to constant eigen row-view
ConstMatrixViewMap MapToEigen(const ConstVectorView& A) {
  return ConstMatrixViewMap(A.mdata + A.mrange.get_start(),
                            A.nelem(),
                            1,
                            StrideType(A.mrange.get_stride(), 1));
}

// Converts constant vector to constant eigen row-view
ConstMatrixViewMap MapToEigenRow(const ConstVectorView& A) {
  return MapToEigen(A);
}

// Converts constant vector to constant eigen column-view
ConstMatrixViewMap MapToEigenCol(const ConstVectorView& A) {
  return ConstMatrixViewMap(A.mdata + A.mrange.get_start(),
                            1,
                            A.nelem(),
                            StrideType(1, A.mrange.get_stride()));
}

// Non- const

// Converts matrix to eigen map
MatrixViewMap MapToEigen(MatrixView& A) {
  return MatrixViewMap(A.mdata + A.mrr.get_start() + A.mcr.get_start(),
                       A.nrows(),
                       A.ncols(),
                       StrideType(A.mrr.get_stride(), A.mcr.get_stride()));
}

// Converts vector to eigen map row-view
MatrixViewMap MapToEigen(VectorView& A) {
  return MatrixViewMap(A.mdata + A.mrange.get_start(),
                       A.nelem(),
                       1,
                       StrideType(A.mrange.get_stride(), 1));
}

// Converts vector to eigen map row-view
MatrixViewMap MapToEigenRow(VectorView& A) { return MapToEigen(A); }

// Converts vector to eigen map column-view
MatrixViewMap MapToEigenCol(VectorView& A) {
  return MatrixViewMap(A.mdata + A.mrange.get_start(),
                       1,
                       A.nelem(),
                       StrideType(1, A.mrange.get_stride()));
}

// Special 4x4

// Converts matrix to eigen map
Matrix4x4ViewMap MapToEigen4x4(MatrixView& A) {
  return Matrix4x4ViewMap(A.mdata + A.mrr.get_start() + A.mcr.get_start(),
                          4,
                          4,
                          StrideType(A.mrr.get_stride(), A.mcr.get_stride()));
}

// Converts constant matrix to constant eigen map
ConstMatrix4x4ViewMap MapToEigen4x4(const ConstMatrixView& A) {
  return ConstMatrix4x4ViewMap(
      A.mdata + A.mrr.get_start() + A.mcr.get_start(),
      4,
      4,
      StrideType(A.mrr.get_stride(), A.mcr.get_stride()));
}

////////////////////////////////
// Helper function for debugging
#ifndef NDEBUG

/** Helper function to access matrix elements.

    Because of function inlining the operator() is not
    accessible from the debuggger. This function helps to access
    Matrix elements from within the debugger.

    \param mv MatrixView
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Numeric debug_matrixview_get_elem(MatrixView& mv, Index r, Index c) {
  return mv(r, c);
}

#endif
////////////////////////////////
