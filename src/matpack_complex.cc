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

/*!
  \file   complex.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-12-19
  
  \brief  A class implementing complex numbers for ARTS.
*/

#include "matpack_complex.h"
#include "blas.h"
#include "exceptions.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Functions for ConstComplexVectorView:
// ------------------------------

//! Returns true if variable size is zero.
// bool ConstComplexVectorView::empty() const
// {
//     return (nelem() == 0);
// }

/** Returns the number of elements.  The names `size' and `length'
    are already used by STL functions returning size_t. To avoid
    confusion we choose the name `nelem'. This is also more
    consistent with `nrow' and `ncol' for matrices.
    
    The value range of long, which is used to store the index is on a
    PC from -2147483648 to 2147483647. This means that a 15GB large
    array of float can be addressed with this index. So the extra bit
    that size_t has compared to long is not needed. */
// Index ConstComplexVectorView::nelem() const
// {
//   return mrange.mextent;
// }

/** The sum of all elements of a Vector. */
Complex ConstComplexVectorView::sum() const {
  Complex s = 0;
  ConstComplexIterator1D i = begin();
  const ConstComplexIterator1D e = end();

  for (; i != e; ++i) s += *i;

  return s;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Vector. This allows
    correct recursive behavior.  */
ConstComplexVectorView ConstComplexVectorView::operator[](
    const Range& r) const {
  return ConstComplexVectorView(mdata, mrange, r);
}

/** Return const iterator to first element. */
ConstComplexIterator1D ConstComplexVectorView::begin() const {
  return ConstComplexIterator1D(mdata + mrange.mstart, mrange.mstride);
}

/** Return const iterator behind last element. */
ConstComplexIterator1D ConstComplexVectorView::end() const {
  return ConstComplexIterator1D(
      mdata + mrange.mstart + (mrange.mextent) * mrange.mstride,
      mrange.mstride);
}

/** Conversion to const 1 column matrix. */
ConstComplexVectorView::operator ConstComplexMatrixView() const {
  return ConstComplexMatrixView(mdata, mrange, Range(0, 1));
}

/** A special constructor, which allows to make a ConstComplexVectorView from
    a scalar.

    This one is a bit tricky: We have to cast away the arguments const
    qualifier, because mdata is not const. This should be safe, since
    there are no non-const methods for ConstComplexVectorView.
*/
ConstComplexVectorView::ConstComplexVectorView(const Complex& a)
    : mrange(0, 1), mdata(&const_cast<Complex&>(a)) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
ConstComplexVectorView::ConstComplexVectorView(Complex* data,
                                               const Range& range)
    : mrange(range), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub ranges from
    sub ranges. That means that the new range has to be interpreted
    relative to the original range. The new range may contain -1 for
    the extent which acts as a joker. However, the used Range
    constructor converts this to an explicit range, consistent with
    the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
ConstComplexVectorView::ConstComplexVectorView(Complex* data,
                                               const Range& p,
                                               const Range& n)
    : mrange(p, n), mdata(data) {
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the vector. The iterators know which part of the vector
    is `active', and also the stride. */
std::ostream& operator<<(std::ostream& os, const ConstComplexVectorView& v) {
  ConstComplexIterator1D i = v.begin();
  const ConstComplexIterator1D end = v.end();

  if (i != end) {
    os << std::setw(3) << *i;
    ++i;
  }
  for (; i != end; ++i) {
    os << " " << std::setw(3) << *i;
  }

  return os;
}

// Functions for ComplexVectorView:
// ------------------------

/** Bail out immediately if somebody tries to create a ComplexVectorView from
 a const Complex*Vector. */
ComplexVectorView::ComplexVectorView(const ComplexVector&) {
  ARTS_ASSERT (false,
      "Creating a ComplexVectorView from a const ComplexVector is not allowed.\n"
      "This is not really a runtime error, but I don't want to start\n"
      "producing direct output from inside matpack. And just exiting is\n"
      "not so nice.\n"
      "If you see this error, there is a bug in the code, not in the\n"
      "ARTS input.")
}

/** Create ComplexVectorView from a ComplexVector. */
ComplexVectorView::ComplexVectorView(ComplexVector& v) {
  mdata = v.mdata;
  mrange = v.mrange;
}

/** Index operator for subrange. We have to also account for the case,
    that *this is already a subrange of a Vector. This allows correct
    recursive behavior.  */
ComplexVectorView ComplexVectorView::operator[](const Range& r) {
  return ComplexVectorView(mdata, mrange, r);
}

/** Return iterator to first element. */
ComplexIterator1D ComplexVectorView::begin() {
  return ComplexIterator1D(mdata + mrange.mstart, mrange.mstride);
}

/** Return iterator behind last element. */
ComplexIterator1D ComplexVectorView::end() {
  return ComplexIterator1D(
      mdata + mrange.mstart + (mrange.mextent) * mrange.mstride,
      mrange.mstride);
}

/** Assignment operator. This copies the data from another VectorView
 t o this Complex*VectorView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this VectorView by
    setting its range. */
ComplexVectorView& ComplexVectorView::operator=(
    const ConstComplexVectorView& v) {
  //  cout << "Assigning VectorView from ConstVectorView.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.mrange.mextent);

  copy(v.begin(), v.end(), begin());

  return *this;
}

/** Assignment from ComplexVectorView to ComplexVectorView. This is a tricky
 one. The problem is that since ComplexVectorView is derived from
 ConstComplexVectorView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
ComplexVectorView& ComplexVectorView::operator=(const ComplexVectorView& v) {
  //  cout << "Assigning ComplexVectorView from ComplexVectorView.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.mrange.mextent);

  copy(v.begin(), v.end(), begin());

  return *this;
}

/** Assignment from ComplexVector. This is important to avoid a bug when
 assigning a Vector to a Complex*VectorView. */
ComplexVectorView& ComplexVectorView::operator=(const ComplexVector& v) {
  //  cout << "Assigning ComplexVectorView from Vector.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.mrange.mextent);

  copy(v.begin(), v.end(), begin());

  return *this;
}

/** Assigning a scalar to a ComplexVectorView will set all elements to this
    value. */
ComplexVectorView& ComplexVectorView::operator=(Complex x) {
  copy(x, begin(), end());
  return *this;
}

/** Multiplication by scalar. */
ComplexVectorView ComplexVectorView::operator*=(Complex x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i *= x;
  return *this;
}

/** Multiplication by scalar. */
ComplexVectorView ComplexVectorView::operator*=(Numeric x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i *= x;
  return *this;
}

/** Division by scalar. */
ComplexVectorView ComplexVectorView::operator/=(Complex x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i /= x;
  return *this;
}

/** Division by scalar. */
ComplexVectorView ComplexVectorView::operator/=(Numeric x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i /= x;
  return *this;
}

/** Addition of scalar. */
ComplexVectorView ComplexVectorView::operator+=(Complex x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i += x;
  return *this;
}

/** Addition of scalar. */
ComplexVectorView ComplexVectorView::operator+=(Numeric x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i += x;
  return *this;
}

/** Subtraction of scalar. */
ComplexVectorView ComplexVectorView::operator-=(Complex x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i -= x;
  return *this;
}

/** Subtraction of scalar. */
ComplexVectorView ComplexVectorView::operator-=(Numeric x) {
  const ComplexIterator1D e = end();
  for (ComplexIterator1D i = begin(); i != e; ++i) *i -= x;
  return *this;
}

/** Element-vise multiplication by another vector. */
ComplexVectorView ComplexVectorView::operator*=(
    const ConstComplexVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstComplexIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i *= *s;
  return *this;
}

/** Element-vise multiplication by another vector. */
ComplexVectorView ComplexVectorView::operator*=(const ConstVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i *= *s;
  return *this;
}

/** Element-vise division by another vector. */
ComplexVectorView ComplexVectorView::operator/=(
    const ConstComplexVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstComplexIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i /= *s;
  return *this;
}

/** Element-vise division by another vector. */
ComplexVectorView ComplexVectorView::operator/=(const ConstVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i /= *s;
  return *this;
}

/** Element-vise addition of another vector. */
ComplexVectorView ComplexVectorView::operator+=(
    const ConstComplexVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstComplexIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i += *s;
  return *this;
}

/** Element-vise addition of another vector. */
ComplexVectorView ComplexVectorView::operator+=(const ConstVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i += *s;
  return *this;
}

/** Element-vise subtraction of another vector. */
ComplexVectorView ComplexVectorView::operator-=(
    const ConstComplexVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstComplexIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i -= *s;
  return *this;
}

/** Element-vise subtraction of another vector. */
ComplexVectorView ComplexVectorView::operator-=(const ConstVectorView& x) {
  ARTS_ASSERT(nelem() == x.nelem());

  ConstIterator1D s = x.begin();

  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();

  for (; i != e; ++i, ++s) *i -= *s;
  return *this;
}

/** Conversion to 1 column matrix. */
ComplexVectorView::operator ComplexMatrixView() {
  // The old version (before 2013-01-18) of this was:
  //    return ConstMatrixView(mdata,mrange,Range(mrange.mstart,1));
  // Bus this was a bug! The problem is that the matrix index operator adds
  // the mstart from both row and columm range object to mdata

  return ComplexMatrixView(mdata, mrange, Range(0, 1));
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  VectorView is not pointing to the beginning of a Vector or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Complex* ComplexVectorView::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT (not (mrange.mstart != 0 || mrange.mstride != 1),
        "A ComplexVectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");
  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  VectorView is not pointing to the beginning of a Vector or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Complex* ComplexVectorView::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT (not (mrange.mstart != 0 || mrange.mstride != 1),
        "A ComplexVectorView can only be converted to a plain C-array if it's pointing to a continuous block of data");
  return mdata;
}

/** A special constructor, which allows to make a VectorView from
    a scalar.
*/
ComplexVectorView::ComplexVectorView(Complex& a) : ConstComplexVectorView(a) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Vector to initialize its
    own VectorView part. */
ComplexVectorView::ComplexVectorView(Complex* data, const Range& range)
    : ConstComplexVectorView(data, range) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub ranges from
    sub ranges. That means that the new range has to be interpreted
    relative to the original range. The new range may contain -1 for
    the extent which acts as a joker. However, the used Range
    constructor converts this to an explicit range, consistent with
    the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
ComplexVectorView::ComplexVectorView(Complex* data,
                                     const Range& p,
                                     const Range& n)
    : ConstComplexVectorView(data, p, n) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can for example copy data between different
    kinds of subvectors. */
void copy(ConstComplexIterator1D origin,
          const ConstComplexIterator1D& end,
          ComplexIterator1D target) {
  if (origin.mstride == 1 && target.mstride == 1)
    memcpy((void*)target.mx,
           (void*)origin.mx,
           sizeof(Complex) * (end.mx - origin.mx));
  else
    for (; origin != end; ++origin, ++target) *target = *origin;
}

/** Copy a scalar to all elements. */
void copy(Complex x, ComplexIterator1D target, const ComplexIterator1D& end) {
  for (; target != end; ++target) *target = x;
}

// Functions for Vector:
// ---------------------

/** Constructor setting size. */
ComplexVector::ComplexVector(Index n)
    : ComplexVectorView(new Complex[n], Range(0, n)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
ComplexVector::ComplexVector(Index n, Complex fill)
    : ComplexVectorView(new Complex[n], Range(0, n)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  const Complex* stop = mdata + n;
  for (Complex* x = mdata; x < stop; ++x) *x = fill;
}

/** Constructor setting size and filling with constant value. */
ComplexVector::ComplexVector(Index n, Numeric fill)
    : ComplexVectorView(new Complex[n], Range(0, n)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  const Complex* stop = mdata + n;
  for (Complex* x = mdata; x < stop; ++x) *x = fill;
}

/** Constructor filling with values. 
 * 
 *   Examples:
 * 
 *   Vector v(1,5,1);  // 1, 2, 3, 4, 5
 *   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
 *   Vector v(5,5,-1); // 5, 4, 3, 2, 1
 */
ComplexVector::ComplexVector(Complex start, Index extent, Complex stride)
    : ComplexVectorView(new Complex[extent], Range(0, extent)) {
  // Fill with values:
  Complex x = start;
  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();
  for (; i != e; ++i) {
    *i = x;
    x += stride;
  }
}

/** Constructor filling with values. 
* 
*   Examples:
* 
*   Vector v(1,5,1);  // 1, 2, 3, 4, 5
*   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
*   Vector v(5,5,-1); // 5, 4, 3, 2, 1
*/
ComplexVector::ComplexVector(Numeric start, Index extent, Complex stride)
    : ComplexVectorView(new Complex[extent], Range(0, extent)) {
  // Fill with values:
  Complex x = start;
  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();
  for (; i != e; ++i) {
    *i = x;
    x += stride;
  }
}

/** Constructor filling with values. 
 * 
 *   Examples:
 * 
 *   Vector v(1,5,1);  // 1, 2, 3, 4, 5
 *   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
 *   Vector v(5,5,-1); // 5, 4, 3, 2, 1
 */
ComplexVector::ComplexVector(Complex start, Index extent, Numeric stride)
    : ComplexVectorView(new Complex[extent], Range(0, extent)) {
  // Fill with values:
  Complex x = start;
  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();
  for (; i != e; ++i) {
    *i = x;
    x += stride;
  }
}

/** Constructor filling with values. 
 * 
 *   Examples:
 * 
 *   Vector v(1,5,1);  // 1, 2, 3, 4, 5
 *   Vector v(1,5,.5); // 1, 1.5, 2, 2.5, 3
 *   Vector v(5,5,-1); // 5, 4, 3, 2, 1
 */
ComplexVector::ComplexVector(Numeric start, Index extent, Numeric stride)
    : ComplexVectorView(new Complex[extent], Range(0, extent)) {
  // Fill with values:
  Complex x = start;
  ComplexIterator1D i = begin();
  const ComplexIterator1D e = end();
  for (; i != e; ++i) {
    *i = x;
    x += stride;
  }
}

/** Copy constructor from ComplexVectorView. This automatically sets the size
 *   and copies the data. The vector created will have start zero and
 *   stride 1, independent on how these parameters are set for the
 *   original. So, what is copied is the data, not the shape
 *   of the selection. */
ComplexVector::ComplexVector(const ConstComplexVectorView& v)
    : ComplexVectorView(new Complex[v.nelem()], Range(0, v.nelem())) {
  copy(v.begin(), v.end(), begin());
}

/** Copy constructor from ComplexVector. This is important to override the
 *   automatically generated shallow constructor. We want deep copies!  */
ComplexVector::ComplexVector(const ComplexVector& v)
    : ComplexVectorView(new Complex[v.nelem()], Range(0, v.nelem())) {
  copy(v.begin(), v.end(), begin());
}

ComplexVector::ComplexVector(const Vector& v) : ComplexVectorView(new Complex[v.nelem()], Range(0, v.nelem())) {
  for (Index i=0; i<nelem(); i++) operator[](i) = Complex(v[i], 0);
}

/** Converting constructor from std::vector. */
ComplexVector::ComplexVector(const std::vector<Complex>& v)
    : ComplexVectorView(new Complex[v.size()], Range(0, v.size())) {
  auto vec_it_end = v.end();
  ComplexIterator1D this_it = this->begin();
  for (auto vec_it = v.begin();
       vec_it != vec_it_end;
       ++vec_it, ++this_it)
    *this_it = *vec_it;
}

/** Converting constructor from std::vector. */
ComplexVector::ComplexVector(const std::vector<Numeric>& v)
    : ComplexVectorView(new Complex[v.size()], Range(0, v.size())) {
  auto vec_it_end = v.end();
  ComplexIterator1D this_it = this->begin();
  for (auto vec_it = v.begin();
       vec_it != vec_it_end;
       ++vec_it, ++this_it)
    *this_it = *vec_it;
}

//! Assignment from another Vector.
/*! 
 * While dimensions of VectorViews can not be adjusted, dimensions of
 * Vectors *can* be adjusted. Hence, the behavior of the assignment
 * operator is different.
 * 
 * In this case the size of the target is automatically adjusted. This
 * is important, so that structures containing Vectors are copied
 * correctly. 
 * 
 * This is a deviation from the old ARTS paradigm that sizes must match
 * exactly before copying!
 * 
 * Note: It is sufficient to have only this one version of the
 * assignment (Vector = Vector). It implicitly covers the cases
 * Vector=VectorView, etc, because there is a default constructor for
 * Vector from VectorView. (See C++ Primer Plus, page 571ff.)
 * 
 * \param v The other vector to copy to this one.
 * 
 * \return This vector, by tradition.
 * 
 * \author Stefan Buehler
 * \date   2002-12-19
 */
ComplexVector& ComplexVector::operator=(const ComplexVector& v) {
  const auto n = v.nelem();
  resize(n);

  auto* in = v.mdata + v.mrange.mstart;
  auto* end = v.mdata + n;
  auto* out = mdata;

  while (in < end) {
    out = in;
    out++;
    in += v.mrange.mstride;
  }

  return *this;
}

//! Assignment operator from Array<Numeric>.
/*!
 * This copies the data from a Array<Numeric> to this VectorView. The
 * size is adjusted automatically.
 * 
 * Array<Numeric> can be useful to collect things in, because there
 * is a .push_back method for it. Then, after collecting we usually
 * have to transfer the content to a Vector. With this assignment
 * operator that's easy.
 * 
 * \param x The array to assign to this.
 * 
 * \return This vector, by tradition.
 * 
 * \author Stefan Buehler
 * \date   2002-12-19
 */
ComplexVector& ComplexVector::operator=(const Array<Complex>& x) {
  resize(x.nelem());
  ComplexVectorView::operator=(x);
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
 *   inherited. */
ComplexVector& ComplexVector::operator=(Complex x) {
  ComplexVectorView::operator=(x);
  return *this;
}

// /** Assignment operator from VectorView. Assignment operators are not
//     inherited. */
// inline Vector& Vector::operator=(const VectorView v)
// {
//   cout << "Assigning Vector from Vector View.\n";
//  // Check that sizes are compatible:
//   ARTS_ASSERT(mrange.mextent==v.mrange.mextent);
//   VectorView::operator=(v);
//   return *this;
// }

/** Resize function. If the size is already correct this function does
 *   nothing. All data is lost after resizing! The new Vector is not
 *   initialized, so it will contain random values.  */
void ComplexVector::resize(Index n) {
  ARTS_ASSERT(0 <= n);
  if (mrange.mextent != n) {
    delete[] mdata;
    mdata = new Complex[n];
    mrange.mstart = 0;
    mrange.mextent = n;
    mrange.mstride = 1;
  }
}

/** Swaps two objects. */
void swap(ComplexVector& v1, ComplexVector& v2) {
  std::swap(v1.mrange, v2.mrange);
  std::swap(v1.mdata, v2.mdata);
}

/** Destructor for ComplexVector. This is important, since Vector uses new to
 *   allocate storage. */
ComplexVector::~ComplexVector() { delete[] mdata; }

// Functions for ConstMatrixView:
// ------------------------------

//! Returns true if variable size is zero.
// bool ConstComplexMatrixView::empty() const
// {
//     return (nrows() == 0 || ncols() == 0);
// }

/** Returns the number of rows. */
// Index ConstComplexMatrixView::nrows() const
// {
//   return mrr.mextent;
// }

/** Returns the number of columns. */
// Index ConstComplexMatrixView::ncols() const
// {
//   return mcr.mextent;
// }

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Matrix. This allows
    correct recursive behavior.  */
ConstComplexMatrixView ConstComplexMatrixView::operator()(
    const Range& r, const Range& c) const {
  return ConstComplexMatrixView(mdata, mrr, mcr, r, c);
}

/** Const index operator returning a column as an object of type
 C onstComplex*VectorView.

    \param r A range of rows.
    \param c Index of selected column */
ConstComplexVectorView ConstComplexMatrixView::operator()(const Range& r,
                                                          Index c) const {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return ConstComplexVectorView(mdata + mcr.mstart + c * mcr.mstride, mrr, r);
}

/** Const index operator returning a row as an object of type
 C onstComplex*VectorView.

    \param r Index of selected row.
    \param c Range of columns */
ConstComplexVectorView ConstComplexMatrixView::operator()(
    Index r, const Range& c) const {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return ConstComplexVectorView(mdata + mrr.mstart + r * mrr.mstride, mcr, c);
}

/** Return const iterator to first row. */
ConstComplexIterator2D ConstComplexMatrixView::begin() const {
  return ConstComplexIterator2D(ConstComplexVectorView(mdata + mrr.mstart, mcr),
                                mrr.mstride);
}

/** Return const iterator behind last row. */
ConstComplexIterator2D ConstComplexMatrixView::end() const {
  return ConstComplexIterator2D(
      ConstComplexVectorView(mdata + mrr.mstart + (mrr.mextent) * mrr.mstride,
                             mcr),
      mrr.mstride);
}

//! ComplexMatrix diagonal as vector.
/*!
 R eturns a ConstComplex*MatrixView on the diagonal entries of the matrix. For a given
  (n,m) matrix M the diagonal vector v is the vector of length min{n,m} with entries

       v[i] = M(i,i)

  \return The diagonal vector v.
*/
ConstComplexVectorView ConstComplexMatrixView::diagonal() const {
  Index n = std::min(mrr.mextent, mcr.mextent);
  return ConstComplexVectorView(mdata + mrr.mstart + mcr.mstart,
                                Range(0, n, mrr.mstride + mcr.mstride));
}

/** Explicit constructor. This one is used by Matrix to initialize its
    own MatrixView part. The row range rr must have a
    stride to account for the length of one row. */
ConstComplexMatrixView::ConstComplexMatrixView(Complex* data,
                                               const Range& rr,
                                               const Range& cr)
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
ConstComplexMatrixView::ConstComplexMatrixView(Complex* data,
                                               const Range& pr,
                                               const Range& pc,
                                               const Range& nr,
                                               const Range& nc)
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
std::ostream& operator<<(std::ostream& os, const ConstComplexMatrixView& v) {
  // Row iterators:
  ConstComplexIterator2D ir = v.begin();
  const ConstComplexIterator2D end_row = v.end();

  if (ir != end_row) {
    ConstComplexIterator1D ic = ir->begin();
    const ConstComplexIterator1D end_col = ir->end();

    if (ic != end_col) {
      os << std::setw(3) << *ic;
      ++ic;
    }
    for (; ic != end_col; ++ic) {
      os << " " << std::setw(3) << *ic;
    }
    ++ir;
  }
  for (; ir != end_row; ++ir) {
    ConstComplexIterator1D ic = ir->begin();
    const ConstComplexIterator1D end_col = ir->end();

    os << "\n";
    if (ic != end_col) {
      os << std::setw(3) << *ic;
      ++ic;
    }
    for (; ic != end_col; ++ic) {
      os << " " << std::setw(3) << *ic;
    }
  }

  return os;
}

// Functions for ComplexMatrixView:
// -------------------------

/** Index operator for subrange. We have to also account for the case,
 t hat *this is already a subrange of a Complex*Matrix. This allows correct
    recursive behavior.  */
ComplexMatrixView ComplexMatrixView::operator()(const Range& r,
                                                const Range& c) {
  return ComplexMatrixView(mdata, mrr, mcr, r, c);
}

/** Index operator returning a column as an object of type ComplexVectorView.

    \param r A range of rows.
    \param c Index of selected column */
ComplexVectorView ComplexMatrixView::operator()(const Range& r, Index c) {
  // Check that c is valid:
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(c < mcr.mextent);

  return ComplexVectorView(mdata + mcr.mstart + c * mcr.mstride, mrr, r);
}

/** Index operator returning a row as an object of type ComplexVectorView.

    \param r Index of selected row.
    \param c Range of columns */
ComplexVectorView ComplexMatrixView::operator()(Index r, const Range& c) {
  // Check that r is valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < mrr.mextent);

  return ComplexVectorView(mdata + mrr.mstart + r * mrr.mstride, mcr, c);
}

/** Return iterator to first row. */
ComplexIterator2D ComplexMatrixView::begin() {
  return ComplexIterator2D(ComplexVectorView(mdata + mrr.mstart, mcr),
                           mrr.mstride);
}

/** Return iterator behind last row. */
ComplexIterator2D ComplexMatrixView::end() {
  return ComplexIterator2D(
      ComplexVectorView(mdata + mrr.mstart + (mrr.mextent) * mrr.mstride, mcr),
      mrr.mstride);
}

//! ComplexMatrix diagonal as vector.
/*!
  Returns a ComplexMatrixView on the diagonal entries of the matrix. For a given
  (n,m) matrix M the diagonal vector v is the vector of length min{n,m} with entries

       v[i] = M(i,i)

  \return The diagonal vector v.
*/
ComplexVectorView ComplexMatrixView::diagonal() {
  Index n = std::min(mrr.mextent, mcr.mextent);
  return ComplexVectorView(mdata + mrr.mstart + mcr.mstart,
                           Range(0, n, mrr.mstride + mcr.mstride));
}

/** Assignment operator. This copies the data from another ComplexMatrixView
    to this ComplexMatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this ComplexMatrixView by
    setting its range. */
ComplexMatrixView& ComplexMatrixView::operator=(
    const ConstComplexMatrixView& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from ComplexMatrixView to ComplexMatrixView. This is a tricky
    one. The problem is that since ComplexMatrixView is derived from
    ConstMatrixView, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
ComplexMatrixView& ComplexMatrixView::operator=(const ComplexMatrixView& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a ComplexMatrix. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
ComplexMatrixView& ComplexMatrixView::operator=(const ComplexMatrix& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a vector. This copies the data from a ComplexVectorView
    to this ComplexMatrixView. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this ComplexMatrixView by
    setting its range. */
ComplexMatrixView& ComplexMatrixView::operator=(
    const ConstComplexVectorView& v) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mrr.mextent == v.nelem());
  ARTS_ASSERT(mcr.mextent == 1);
  //  dummy = ConstComplexMatrixView(v.mdata,v.mrange,Range(v.mrange.mstart,1));;
  ConstComplexMatrixView dummy(v);
  copy(dummy.begin(), dummy.end(), begin());
  return *this;
}

/** Assigning a scalar to a MatrixView will set all elements to this
    value. */
ComplexMatrixView& ComplexMatrixView::operator=(Complex x) {
  copy(x, begin(), end());
  return *this;
}

/** Multiplication by scalar. */
ComplexMatrixView& ComplexMatrixView::operator*=(Complex x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c *= x;
  }
  return *this;
}

/** Multiplication by scalar. */
ComplexMatrixView& ComplexMatrixView::operator*=(Numeric x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c *= x;
  }
  return *this;
}

/** Division by scalar. */
ComplexMatrixView& ComplexMatrixView::operator/=(Complex x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c /= x;
  }
  return *this;
}

/** Division by scalar. */
ComplexMatrixView& ComplexMatrixView::operator/=(Numeric x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c /= x;
  }
  return *this;
}

/** Addition of scalar. */
ComplexMatrixView& ComplexMatrixView::operator+=(Complex x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c += x;
  }
  return *this;
}

/** Addition of scalar. */
ComplexMatrixView& ComplexMatrixView::operator+=(Numeric x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c += x;
  }
  return *this;
}

/** Subtraction of scalar. */
ComplexMatrixView& ComplexMatrixView::operator-=(Complex x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c -= x;
  }
  return *this;
}

/** Subtraction of scalar. */
ComplexMatrixView& ComplexMatrixView::operator-=(Numeric x) {
  const ComplexIterator2D er = end();
  for (ComplexIterator2D r = begin(); r != er; ++r) {
    const ComplexIterator1D ec = r->end();
    for (ComplexIterator1D c = r->begin(); c != ec; ++c) *c -= x;
  }
  return *this;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Complex* ComplexMatrixView::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT (not (mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1),
        "A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");
  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  MatrixView is not pointing to the beginning of a Matrix or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
Complex* ComplexMatrixView::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT (not (mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1),
        "A MatrixView can only be converted to a plain C-array if it's pointing to a continuous block of data");
  return mdata;
}

/** Element-vise multiplication by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator*=(
    const ConstComplexMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstComplexIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstComplexIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c *= *sc;
  }
  return *this;
}

/** Element-vise multiplication by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator*=(const ConstMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c *= *sc;
  }
  return *this;
}

/** Element-vise division by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator/=(
    const ConstComplexMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstComplexIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstComplexIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c /= *sc;
  }
  return *this;
}

/** Element-vise division by another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator/=(const ConstMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c /= *sc;
  }
  return *this;
}

/** Element-vise addition of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator+=(
    const ConstComplexMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstComplexIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstComplexIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c += *sc;
  }
  return *this;
}

/** Element-vise addition of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator+=(const ConstMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c += *sc;
  }
  return *this;
}

/** Element-vise subtraction of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator-=(
    const ConstComplexMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstComplexIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstComplexIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c -= *sc;
  }
  return *this;
}

/** Element-vise subtraction of another Matrix. */
ComplexMatrixView& ComplexMatrixView::operator-=(const ConstMatrixView& x) {
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator2D sr = x.begin();
  ComplexIterator2D r = begin();
  const ComplexIterator2D er = end();
  for (; r != er; ++r, ++sr) {
    ConstIterator1D sc = sr->begin();
    ComplexIterator1D c = r->begin();
    const ComplexIterator1D ec = r->end();
    for (; c != ec; ++c, ++sc) *c -= *sc;
  }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Matrix. */
ComplexMatrixView::ComplexMatrixView() : ConstComplexMatrixView() {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by ComplexMatrix to initialize its
    own ComplexMatrixView part. The row range rr must have a
    stride to account for the length of one row. */
ComplexMatrixView::ComplexMatrixView(Complex* data,
                                     const Range& rr,
                                     const Range& cr)
    : ConstComplexMatrixView(data, rr, cr) {
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
ComplexMatrixView::ComplexMatrixView(Complex* data,
                                     const Range& pr,
                                     const Range& pc,
                                     const Range& nr,
                                     const Range& nc)
    : ConstComplexMatrixView(data, pr, pc, nr, nc) {
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
void copy(ConstComplexIterator2D origin,
          const ConstComplexIterator2D& end,
          ComplexIterator2D target) {
  for (; origin != end; ++origin, ++target) {
    ConstComplexIterator1D o = origin->begin();
    const ConstComplexIterator1D e = origin->end();
    ComplexIterator1D t = target->begin();
    for (; o != e; ++o, ++t) *t = *o;
  }
}

/** Copy a scalar to all elements. */
void copy(Complex x, ComplexIterator2D target, const ComplexIterator2D& end) {
  for (; target != end; ++target) {
    ComplexIterator1D t = target->begin();
    const ComplexIterator1D e = target->end();
    for (; t != e; ++t) *t = x;
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, ComplexIterator2D target, const ComplexIterator2D& end) {
  for (; target != end; ++target) {
    ComplexIterator1D t = target->begin();
    const ComplexIterator1D e = target->end();
    for (; t != e; ++t) *t = x;
  }
}

// Functions for ComplexMatrix:
// ---------------------

/** Constructor setting size. This constructor has to set the stride
    in the row range correctly! */
ComplexMatrix::ComplexMatrix(Index r, Index c)
    : ComplexMatrixView(new Complex[r * c], Range(0, r, c), Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
ComplexMatrix::ComplexMatrix(Index r, Index c, Complex fill)
    : ComplexMatrixView(new Complex[r * c], Range(0, r, c), Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  const Complex* stop = mdata + r * c;
  for (Complex* x = mdata; x < stop; ++x) *x = fill;
}

/** Constructor setting size and filling with constant value. */
ComplexMatrix::ComplexMatrix(Index r, Index c, Numeric fill)
    : ComplexMatrixView(new Complex[r * c], Range(0, r, c), Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  const Complex* stop = mdata + r * c;
  for (Complex* x = mdata; x < stop; ++x) *x = fill;
}

/** Copy constructor from MatrixView. This automatically sets the size
    and copies the data. */
ComplexMatrix::ComplexMatrix(const ConstComplexMatrixView& m)
    : ComplexMatrixView(new Complex[m.nrows() * m.ncols()],
                        Range(0, m.nrows(), m.ncols()),
                        Range(0, m.ncols())) {
  copy(m.begin(), m.end(), begin());
}

/** Copy constructor from Matrix. This automatically sets the size
    and copies the data. */
ComplexMatrix::ComplexMatrix(const ComplexMatrix& m)
    : ComplexMatrixView(new Complex[m.nrows() * m.ncols()],
                        Range(0, m.nrows(), m.ncols()),
                        Range(0, m.ncols())) {
  // There is a catch here: If m is an empty matrix, then it will have
  // 0 colunns. But this is used to initialize the stride of the row
  // Range! Thus, this method has to be consistent with the behaviour
  // of Range::Range. For now, Range::Range allows also stride 0.
  copy(m.begin(), m.end(), begin());
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
ComplexMatrix& ComplexMatrix::operator=(ComplexMatrix m) {
  swap(*this, m);
  return *this;
}

/** Assignment operator from scalar. Assignment operators also seem to
    be not inherited. */
ComplexMatrix& ComplexMatrix::operator=(Complex x) {
  copy(x, begin(), end());
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
ComplexMatrix& ComplexMatrix::operator=(const ConstComplexVectorView& v) {
  resize(v.nelem(), 1);
  ConstComplexMatrixView dummy(v);
  copy(dummy.begin(), dummy.end(), begin());
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new Matrix is not
    initialized, so it will contain random values.*/
void ComplexMatrix::resize(Index r, Index c) {
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  if (mrr.mextent != r || mcr.mextent != c) {
    delete[] mdata;
    mdata = new Complex[r * c];

    mrr.mstart = 0;
    mrr.mextent = r;
    mrr.mstride = c;

    mcr.mstart = 0;
    mcr.mextent = c;
    mcr.mstride = 1;
  }
}

/** Swaps two objects. */
void swap(ComplexMatrix& m1, ComplexMatrix& m2) {
  std::swap(m1.mrr, m2.mrr);
  std::swap(m1.mcr, m2.mcr);
  std::swap(m1.mdata, m2.mdata);
}

/** Destructor for Matrix. This is important, since Matrix uses new to
    allocate storage. */
ComplexMatrix::~ComplexMatrix() {
  //   cout << "Destroying a Matrix:\n"
  //        << *this << "\n........................................\n";
  delete[] mdata;
}

ComplexMatrix& ComplexMatrix::inv(const lapack_help::Inverse<Complex>& help)
{
  ARTS_ASSERT(ncols() == nrows());
  
  // help's internal variables are all mutable, so const is just for default parameter and to not use copies
  help.resize_if_smaller(int(ncols()));
  
  // Compute LU decomposition using LAPACK dgetrf_.
  int info;
  lapack::zgetrf_(help.size(), help.size(), mdata, help.size(), help.ipivdata(), &info);
  lapack::zgetri_(help.size(), mdata, help.size(), help.ipivdata(), help.workdata(), help.lsize(), &info);
  
  // Check for success.
  ARTS_USER_ERROR_IF (info not_eq 0,
    "Error inverting matrix: Matrix not of full rank.");
  
  return *this;
}

/** Const version of transpose. */
ConstComplexMatrixView transpose(ConstComplexMatrixView m) {
  return ConstComplexMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
    transpose. The original is not changed! */
ComplexMatrixView transpose(ComplexMatrixView m) {
  return ComplexMatrixView(m.mdata, m.mcr, m.mrr);
}

/** Returns the transpose. This creates a special MatrixView for the
 transpose. The original is not changed! */
ComplexMatrixView transpose(ComplexVector v) {
  return transpose((ComplexMatrixView)v);
}

/** Assignment operator from Array<Complex>. This copies the data from
    an Array<Complex> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Complex> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy. */
ComplexVectorView& ComplexVectorView::operator=(const Array<Complex>& v) {
  //  cout << "Assigning VectorView from Array<Numeric>.\n";

  // Check that sizes are compatible:
  ARTS_ASSERT(mrange.mextent == v.nelem());

  // Iterators for Array:
  auto i = v.begin();
  const auto e = v.end();
  // Iterator for Vector:
  ComplexIterator1D target = begin();

  for (; i != e; ++i, ++target) *target = *i;

  return *this;
}

//Functions operating on the complex vectors

/** Scalar product. The two vectors may be identical. */
Complex operator*(const ConstComplexVectorView& a,
                  const ConstComplexVectorView& b) {
  // Check dimensions:
  ARTS_ASSERT(a.nelem() == b.nelem());

  const ConstComplexIterator1D ae = a.end();
  ConstComplexIterator1D ai = a.begin();
  ConstComplexIterator1D bi = b.begin();

  Complex res = 0;
  for (; ai != ae; ++ai, ++bi) res += (*ai) * (*bi);

  return res;
}

//! Matrix-Vector Multiplication
/*!
 * Uses the Eigen library.  Be carful to test the size of your input beforehand.
 * 
 * For left-hand multiplication, please use pure matrix-mult.
 * 
 * \param[out] y The length-m ComplexVectorView where the result is stored.
 * \param[in] M Reference to the m-times-n Const{Complex,}MatrixView holding the matrix M.
 * \param[in] x Reference to the length-n Const{Complex,}VectorView holding the vector x.
 */
void mult(ComplexVectorView y,
          const ConstComplexMatrixView& M,
          const ConstComplexVectorView& x) {
  ARTS_ASSERT(x.nelem() == M.nrows());
  ARTS_ASSERT(y.nelem() == M.ncols());

  ComplexMatrixViewMap eigen_y = MapToEigenRow(y);
  if (y.mdata == x.mdata)
    eigen_y = MapToEigen(M) * MapToEigenRow(x);
  else
    eigen_y.noalias() = MapToEigen(M) * MapToEigenRow(x);
}

//! Matrix-Matrix Multiplication
/*!
 * Uses the Eigen library.  Be carful to test the size of your input beforehand.
 * Note that to keep speed, the inputs should be different variables.  There is
 * memory duplication if this is not the case.  Note that it is mdata that is 
 * checked, so even if the matrices are at different parts of a tensor, there
 * is still a slowdown
 * 
 * \param[in,out] A The matrix A, that will hold the result of the multiplication.
 * \param[in] B The matrix B
 * \param[in] C The matrix C
 */
void mult(ComplexMatrixView A,
          const ConstComplexMatrixView& B,
          const ConstComplexMatrixView& C) {
  ARTS_ASSERT(B.nrows() == A.nrows());
  ARTS_ASSERT(C.ncols() == A.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  ComplexMatrixViewMap eigen_A = MapToEigen(A);
  if (A.mdata == B.mdata || A.mdata == C.mdata)
    eigen_A = MapToEigen(B) * MapToEigen(C);
  else
    eigen_A.noalias() = MapToEigen(B) * MapToEigen(C);
}
void mult(ComplexMatrixView A,
          const ConstComplexMatrixView& B,
          const ConstMatrixView& C) {
  ARTS_ASSERT(B.nrows() == A.nrows());
  ARTS_ASSERT(C.ncols() == A.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  ComplexMatrixViewMap eigen_A = MapToEigen(A);
  if (A.mdata == B.mdata)
    eigen_A = MapToEigen(B) * MapToEigen(C);
  else
    eigen_A.noalias() = MapToEigen(B) * MapToEigen(C);
}
void mult(ComplexMatrixView A,
          const ConstMatrixView& B,
          const ConstComplexMatrixView& C) {
  ARTS_ASSERT(B.nrows() == A.nrows());
  ARTS_ASSERT(C.ncols() == A.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  ComplexMatrixViewMap eigen_A = MapToEigen(A);
  if (A.mdata == C.mdata)
    eigen_A = MapToEigen(B) * MapToEigen(C);
  else
    eigen_A.noalias() = MapToEigen(B) * MapToEigen(C);
}
void mult(ComplexMatrixView A,
          const ConstMatrixView& B,
          const ConstMatrixView& C) {
  ARTS_ASSERT(B.nrows() == A.nrows());
  ARTS_ASSERT(C.ncols() == A.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  ComplexMatrixViewMap eigen_A = MapToEigen(A);
  eigen_A.noalias() = MapToEigen(B) * MapToEigen(C);
}

// Converts constant matrix to constant eigen map
ConstComplexMatrixViewMap MapToEigen(const ConstComplexMatrixView& A) {
  return ConstComplexMatrixViewMap(
      A.mdata + A.mrr.get_start() + A.mcr.get_start(),
      A.nrows(),
      A.ncols(),
      StrideType(A.mrr.get_stride(), A.mcr.get_stride()));
}

// Converts constant vector to constant eigen row-view
ConstComplexMatrixViewMap MapToEigen(const ConstComplexVectorView& A) {
  return ConstComplexMatrixViewMap(A.mdata + A.mrange.get_start(),
                                   A.nelem(),
                                   1,
                                   StrideType(A.mrange.get_stride(), 1));
}

// Converts constant vector to constant eigen row-view
ConstComplexMatrixViewMap MapToEigenRow(const ConstComplexVectorView& A) {
  return MapToEigen(A);
}

// Converts constant vector to constant eigen column-view
ConstComplexMatrixViewMap MapToEigenCol(const ConstComplexVectorView& A) {
  return ConstComplexMatrixViewMap(A.mdata + A.mrange.get_start(),
                                   1,
                                   A.nelem(),
                                   StrideType(1, A.mrange.get_stride()));
}

// Converts matrix to eigen map
ComplexMatrixViewMap MapToEigen(ComplexMatrixView& A) {
  return ComplexMatrixViewMap(
      A.mdata + A.mrr.get_start() + A.mcr.get_start(),
      A.nrows(),
      A.ncols(),
      StrideType(A.mrr.get_stride(), A.mcr.get_stride()));
}

// Converts vector to eigen map row-view
ComplexMatrixViewMap MapToEigen(ComplexVectorView& A) {
  return ComplexMatrixViewMap(A.mdata + A.mrange.get_start(),
                              A.nelem(),
                              1,
                              StrideType(A.mrange.get_stride(), 1));
}

// Converts vector to eigen map row-view
ComplexMatrixViewMap MapToEigenRow(ComplexVectorView& A) {
  return MapToEigen(A);
}

// Converts vector to eigen map column-view
ComplexMatrixViewMap MapToEigenCol(ComplexVectorView& A) {
  return ComplexMatrixViewMap(A.mdata + A.mrange.get_start(),
                              1,
                              A.nelem(),
                              StrideType(1, A.mrange.get_stride()));
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
Complex debug_matrixview_get_elem(ComplexMatrixView& mv, Index r, Index c) {
  return mv(r, c);
}

#endif
