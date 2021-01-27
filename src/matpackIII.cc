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
   \file   matpackIII.cc

   \author Oliver Lemke
   \date   2002-11-21
*/

#include "matpackIII.h"
#include "exceptions.h"

// Functions for ConstTensor3View:
// ------------------------------

//! Check if variable is empty.
/*!
 \param[in]  x The variable to check.
 \return True if the size of any dimension of x is 0.
 */
bool ConstTensor3View::empty() const {
  return (npages() == 0 || nrows() == 0 || ncols() == 0);
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior.  */
ConstTensor3View ConstTensor3View::operator()(const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  return ConstTensor3View(mdata, mpr, mrr, mcr, p, r, c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
ConstMatrixView ConstTensor3View::operator()(const Range& p,
                                             const Range& r,
                                             Index c) const {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstMatrixView(mdata + mcr.mstart + c * mcr.mstride, mpr, mrr, p, r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
ConstMatrixView ConstTensor3View::operator()(const Range& p,
                                             Index r,
                                             const Range& c) const {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstMatrixView(mdata + mrr.mstart + r * mrr.mstride, mpr, mcr, p, c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by one.) */
ConstMatrixView ConstTensor3View::operator()(Index p,
                                             const Range& r,
                                             const Range& c) const {
  // Check that p is valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(p < mpr.mextent);

  return ConstMatrixView(mdata + mpr.mstart + p * mpr.mstride, mrr, mcr, r, c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) */
ConstVectorView ConstTensor3View::operator()(Index p,
                                             Index r,
                                             const Range& c) const {
  // Check that p and r are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstVectorView(
      mdata + mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride,
      mcr,
      c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) */
ConstVectorView ConstTensor3View::operator()(Index p,
                                             const Range& r,
                                             Index c) const {
  // Check that p and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(
      mdata + mpr.mstart + p * mpr.mstride + mcr.mstart + c * mcr.mstride,
      mrr,
      r);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by two.) */
ConstVectorView ConstTensor3View::operator()(const Range& p,
                                             Index r,
                                             Index c) const {
  // Check that r and c are valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstVectorView(
      mdata + mrr.mstart + r * mrr.mstride + mcr.mstart + c * mcr.mstride,
      mpr,
      p);
}

/** Return const iterator to first page. */
ConstIterator3D ConstTensor3View::begin() const {
  return ConstIterator3D(ConstMatrixView(mdata + mpr.mstart, mrr, mcr),
                         mpr.mstride);
}

/** Return const iterator behind last page. */
ConstIterator3D ConstTensor3View::end() const {
  return ConstIterator3D(
      ConstMatrixView(
          mdata + mpr.mstart + (mpr.mextent) * mpr.mstride, mrr, mcr),
      mpr.mstride);
}

/** Special constructor to make a Tensor3 view of a matrix. */
ConstTensor3View::ConstTensor3View(const ConstMatrixView& a)
    : mpr(0, 1, a.mrr.mextent * a.mcr.mextent),
      mrr(a.mrr),
      mcr(a.mcr),
      mdata(a.mdata) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor3 to initialize
    its own Tensor3View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
ConstTensor3View::ConstTensor3View(Numeric* data,
                                   const Range& pr,
                                   const Range& rr,
                                   const Range& cr)
    : mpr(pr), mrr(rr), mcr(cr), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
ConstTensor3View::ConstTensor3View(Numeric* data,
                                   const Range& pp,
                                   const Range& pr,
                                   const Range& pc,
                                   const Range& np,
                                   const Range& nr,
                                   const Range& nc)
    : mpr(pp, np), mrr(pr, nr), mcr(pc, nc), mdata(data) {
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Matrix to print each page in turn. */
std::ostream& operator<<(std::ostream& os, const ConstTensor3View& v) {
  // Page iterators:
  ConstIterator3D ip = v.begin();
  const ConstIterator3D end_page = v.end();

  if (ip != end_page) {
    os << *ip;
    ++ip;
  }

  for (; ip != end_page; ++ip) {
    os << "\n\n";
    os << *ip;
  }

  return os;
}

// Functions for Tensor3View:
// -------------------------

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor3. This allows
    correct recursive behavior.  */
Tensor3View Tensor3View::operator()(const Range& p,
                                    const Range& r,
                                    const Range& c) {
  return Tensor3View(mdata, mpr, mrr, mcr, p, r, c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
MatrixView Tensor3View::operator()(const Range& p, const Range& r, Index c) {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return MatrixView(mdata + mcr.mstart + c * mcr.mstride, mpr, mrr, p, r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
MatrixView Tensor3View::operator()(const Range& p, Index r, const Range& c) {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return MatrixView(mdata + mrr.mstart + r * mrr.mstride, mpr, mcr, p, c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by one.) */
MatrixView Tensor3View::operator()(Index p, const Range& r, const Range& c) {
  // Check that p is valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(p < mpr.mextent);

  return MatrixView(mdata + mpr.mstart + p * mpr.mstride, mrr, mcr, r, c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by two.) */
VectorView Tensor3View::operator()(Index p, Index r, const Range& c) {
  // Check that p and r are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(r < mrr.mextent);

  return VectorView(
      mdata + mpr.mstart + p * mpr.mstride + mrr.mstart + r * mrr.mstride,
      mcr,
      c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by two.) */
VectorView Tensor3View::operator()(Index p, const Range& r, Index c) {
  // Check that p and c are valid:
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(p < mpr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(
      mdata + mpr.mstart + p * mpr.mstride + mcr.mstart + c * mcr.mstride,
      mrr,
      r);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by two.) */
VectorView Tensor3View::operator()(const Range& p, Index r, Index c) {
  // Check that r and r are valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < mrr.mextent);
  ARTS_ASSERT(c < mcr.mextent);

  return VectorView(
      mdata + mrr.mstart + r * mrr.mstride + mcr.mstart + c * mcr.mstride,
      mpr,
      p);
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  Tensor3View is not pointing to the beginning of a Tensor3 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric* Tensor3View::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT (not (mpr.mstart != 0 || mpr.mstride != mrr.mextent * mcr.mextent ||
      mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1),
        "A Tensor3View can only be converted to a plain C-array if it's pointing to a continuous block of data");
  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  Tensor3View is not pointing to the beginning of a Tensor3 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric* Tensor3View::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT (not (mpr.mstart != 0 || mpr.mstride != mrr.mextent * mcr.mextent ||
  mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
  mcr.mstride != 1),
  "A Tensor3View can only be converted to a plain C-array if it's pointing to a continuous block of data");
  return mdata;
}

/** Return iterator to first page. */
Iterator3D Tensor3View::begin() {
  return Iterator3D(MatrixView(mdata + mpr.mstart, mrr, mcr), mpr.mstride);
}

/** Return iterator behind last page. */
Iterator3D Tensor3View::end() {
  return Iterator3D(
      MatrixView(mdata + mpr.mstart + (mpr.mextent) * mpr.mstride, mrr, mcr),
      mpr.mstride);
}

/** Assignment operator. This copies the data from another Tensor3View
    to this Tensor3View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor3View by
    setting its range. */
Tensor3View& Tensor3View::operator=(const ConstTensor3View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from Tensor3View to Tensor3View. This is a tricky
    one. The problem is that since Tensor3View is derived from
    ConstTensor3View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
Tensor3View& Tensor3View::operator=(const Tensor3View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a Tensor3. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
Tensor3View& Tensor3View::operator=(const Tensor3& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assigning a scalar to a Tensor3View will set all elements to this
    value. */
Tensor3View& Tensor3View::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

// Some little helper functions:
//------------------------------

Numeric add(Numeric x, Numeric y) { return x + y; }

/** Multiplication by scalar. */
Tensor3View& Tensor3View::operator*=(Numeric x) {
  const Iterator3D ep = end();
  for (Iterator3D p = begin(); p != ep; ++p) {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
Tensor3View& Tensor3View::operator/=(Numeric x) {
  const Iterator3D ep = end();
  for (Iterator3D p = begin(); p != ep; ++p) {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
Tensor3View& Tensor3View::operator+=(Numeric x) {
  const Iterator3D ep = end();
  for (Iterator3D p = begin(); p != ep; ++p) {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
Tensor3View& Tensor3View::operator-=(Numeric x) {
  const Iterator3D ep = end();
  for (Iterator3D p = begin(); p != ep; ++p) {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor3. */
Tensor3View& Tensor3View::operator*=(const ConstTensor3View& x) {
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator3D xp = x.begin();
  Iterator3D p = begin();
  const Iterator3D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p *= *xp;
  }
  return *this;
}

/** Element-vise division by another Tensor3. */
Tensor3View& Tensor3View::operator/=(const ConstTensor3View& x) {
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator3D xp = x.begin();
  Iterator3D p = begin();
  const Iterator3D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p /= *xp;
  }
  return *this;
}

/** Element-vise addition of another Tensor3. */
Tensor3View& Tensor3View::operator+=(const ConstTensor3View& x) {
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator3D xp = x.begin();
  Iterator3D p = begin();
  const Iterator3D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p += *xp;
  }
  return *this;
}

/** Element-vise subtraction of another Tensor3. */
Tensor3View& Tensor3View::operator-=(const ConstTensor3View& x) {
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator3D xp = x.begin();
  Iterator3D p = begin();
  const Iterator3D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p -= *xp;
  }
  return *this;
}

/** Special constructor to make a Tensor3 view of a matrix. */
Tensor3View::Tensor3View(const MatrixView& a)
    : ConstTensor3View(
          a.mdata, Range(0, 1, a.mrr.mextent * a.mcr.mextent), a.mrr, a.mcr) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor3 to initialize its
    own Tensor3View part. The row range rr must have a
    stride to account for the length of one row. */
Tensor3View::Tensor3View(Numeric* data,
                         const Range& pr,
                         const Range& rr,
                         const Range& cr)
    : ConstTensor3View(data, pr, rr, cr) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pp Previous range.
    \param pr Previous range.
    \param pc Previous range.
    \param np New Range.
    \param nr New Range.
    \param nc New Range.
  */
Tensor3View::Tensor3View(Numeric* data,
                         const Range& pp,
                         const Range& pr,
                         const Range& pc,
                         const Range& np,
                         const Range& nr,
                         const Range& nc)
    : ConstTensor3View(data, pp, pr, pc, np, nr, nc) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
void copy(ConstIterator3D origin,
          const ConstIterator3D& end,
          Iterator3D target) {
  for (; origin != end; ++origin, ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(origin->begin(), origin->end(), target->begin());
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, Iterator3D target, const Iterator3D& end) {
  for (; target != end; ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(x, target->begin(), target->end());
  }
}

// Functions for Tensor3:
// ---------------------

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
Tensor3::Tensor3(Index p, Index r, Index c)
    : Tensor3View(new Numeric[p * r * c],
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Tensor3::Tensor3(Index p, Index r, Index c, Numeric fill)
    : Tensor3View(new Numeric[p * r * c],
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, p * r * c, fill);
}

/** Copy constructor from Tensor3View. This automatically sets the size
    and copies the data. */
Tensor3::Tensor3(const ConstTensor3View& m)
    : Tensor3View(new Numeric[m.npages() * m.nrows() * m.ncols()],
                  Range(0, m.npages(), m.nrows() * m.ncols()),
                  Range(0, m.nrows(), m.ncols()),
                  Range(0, m.ncols())) {
  copy(m.begin(), m.end(), begin());
}

/** Copy constructor from Tensor3. This automatically sets the size
    and copies the data. */
Tensor3::Tensor3(const Tensor3& m)
    : Tensor3View(new Numeric[m.npages() * m.nrows() * m.ncols()],
                  Range(0, m.npages(), m.nrows() * m.ncols()),
                  Range(0, m.nrows(), m.ncols()),
                  Range(0, m.ncols())) {
  // There is a catch here: If m is an empty tensor, then it will have
  // dimensions of size 0. But these are used to initialize the stride
  // for higher dimensions! Thus, this method has to be consistent
  // with the behaviour of Range::Range. For now, Range::Range allows
  // also stride 0.
  std::memcpy(mdata, m.mdata, npages() * nrows() * ncols() * sizeof(Numeric));
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
Tensor3& Tensor3::operator=(const Tensor3& x) {
  if (this != &x) {
    resize(x.npages(), x.nrows(), x.ncols());
    std::memcpy(mdata, x.mdata, npages() * nrows() * ncols() * sizeof(Numeric));
  }
  return *this;
}

//! Move assignment operator from another tensor.
Tensor3& Tensor3::operator=(Tensor3&& x) noexcept {
  if (this != &x) {
    delete[] mdata;
    mdata = x.mdata;
    mpr = x.mpr;
    mrr = x.mrr;
    mcr = x.mcr;
    x.mpr = Range(0, 0);
    x.mrr = Range(0, 0);
    x.mcr = Range(0, 0);
    x.mdata = nullptr;
  }
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
Tensor3& Tensor3::operator=(Numeric x) {
  std::fill_n(mdata, npages() * nrows() * ncols(), x);
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
void Tensor3::resize(Index p, Index r, Index c) {
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  if (mpr.mextent != p || mrr.mextent != r || mcr.mextent != c) {
    delete[] mdata;
    mdata = new Numeric[p * r * c];

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
void swap(Tensor3& t1, Tensor3& t2) {
  std::swap(t1.mpr, t2.mpr);
  std::swap(t1.mrr, t2.mrr);
  std::swap(t1.mcr, t2.mcr);
  std::swap(t1.mdata, t2.mdata);
}

/** Destructor for Tensor3. This is important, since Tensor3 uses new to
    allocate storage. */
Tensor3::~Tensor3() {
  //   cout << "Destroying a Tensor3:\n"
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
void transform(Tensor3View y, double (&my_func)(double), ConstTensor3View x) {
  // Check dimensions:
  ARTS_ASSERT(y.npages() == x.npages());
  ARTS_ASSERT(y.nrows() == x.nrows());
  ARTS_ASSERT(y.ncols() == x.ncols());

  const ConstIterator3D xe = x.end();
  ConstIterator3D xi = x.begin();
  Iterator3D yi = y.begin();
  for (; xi != xe; ++xi, ++yi) {
    // Use the transform function of lower dimensional tensors
    // recursively:
    transform(*yi, my_func, *xi);
  }
}

/** Max function, tensor version. */
Numeric max(const ConstTensor3View& x) {
  const ConstIterator3D xe = x.end();
  ConstIterator3D xi = x.begin();

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
Numeric min(const ConstTensor3View& x) {
  const ConstIterator3D xe = x.end();
  ConstIterator3D xi = x.begin();

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
    \param p  Page index
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Numeric debug_tensor3view_get_elem(Tensor3View& tv, Index p, Index r, Index c) {
  return tv(p, r, c);
}

#endif
////////////////////////////////

//! mult Tensor3
/*!
    Pointwise multiplication of a vector element and matrix.

    mult(Tensor3& A, const ConstVectorView B, const ConstMatrixView C) for multiplying vector pointwise with matrix.
        Useful, e.g., for frequency gridded absorption
        vector multiplied by normalized Stokes extinction
        matrix to get extinction matrix as a function of frequency.

    \param   A   Out: Tensor3 with N pages, M rows and L columns
    \param   B   In: A Vector of length N.
    \param   C   In: A Matrix of size M x L.

    \author Richard Larsson
    \date   2012-07-17
*/
void mult(Tensor3View A, const ConstVectorView B, const ConstMatrixView C) {
  ARTS_ASSERT(A.npages() == B.nelem());
  ARTS_ASSERT(A.nrows() == C.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());

  for (Index ii = 0; ii < B.nelem(); ii++) {
    A(ii, joker, joker) = C;
    A(ii, joker, joker) *= B[ii];
  }
}
