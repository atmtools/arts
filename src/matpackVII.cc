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
   \file   matpackVII.cc

   \author Oliver Lemke
   \date   2002-11-21
*/

#include "matpackVII.h"
#include "exceptions.h"

// Functions for ConstTensor7View:
// ------------------------------

//! Check if variable is empty.
/*!
 \param[in]  x The variable to check.
 \return True if the size of any dimension of x is 0.
 */
bool ConstTensor7View::empty() const {
  return (nlibraries() == 0 || nvitrines() == 0 || nshelves() == 0 ||
          nbooks() == 0 || npages() == 0 || nrows() == 0 || ncols() == 0);
}

/** Returns the number of libraries. */
Index ConstTensor7View::nlibraries() const { return mlr.mextent; }

/** Returns the number of vitrines. */
Index ConstTensor7View::nvitrines() const { return mvr.mextent; }

/** Returns the number of shelves. */
Index ConstTensor7View::nshelves() const { return msr.mextent; }

/** Returns the number of books. */
Index ConstTensor7View::nbooks() const { return mbr.mextent; }

/** Returns the number of pages. */
Index ConstTensor7View::npages() const { return mpr.mextent; }

/** Returns the number of rows. */
Index ConstTensor7View::nrows() const { return mrr.mextent; }

/** Returns the number of columns. */
Index ConstTensor7View::ncols() const { return mcr.mextent; }

// Const index operators:

// -------
ConstTensor7View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  return ConstTensor7View(
      mdata, mlr, mvr, msr, mbr, mpr, mrr, mcr, l, v, s, b, p, r, c);
}
// |------
ConstTensor6View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  return ConstTensor6View(
      mdata + OFFSET(l), mvr, msr, mbr, mpr, mrr, mcr, v, s, b, p, r, c);
}

// ------|
ConstTensor6View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(c);
  return ConstTensor6View(
      mdata + OFFSET(c), mlr, mvr, msr, mbr, mpr, mrr, l, v, s, b, p, r);
}
// |-----|
ConstTensor5View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(l) + OFFSET(c), mvr, msr, mbr, mpr, mrr, v, s, b, p, r);
}

// -----|-
ConstTensor6View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(r);
  return ConstTensor6View(
      mdata + OFFSET(r), mlr, mvr, msr, mbr, mpr, mcr, l, v, s, b, p, c);
}
// |----|-
ConstTensor5View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(r);
  return ConstTensor5View(
      mdata + OFFSET(l) + OFFSET(r), mvr, msr, mbr, mpr, mcr, v, s, b, p, c);
}

// ----|--
ConstTensor6View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(p);
  return ConstTensor6View(
      mdata + OFFSET(p), mlr, mvr, msr, mbr, mrr, mcr, l, v, s, b, r, c);
}
// |---|--
ConstTensor5View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(p);
  return ConstTensor5View(
      mdata + OFFSET(l) + OFFSET(p), mvr, msr, mbr, mrr, mcr, v, s, b, r, c);
}

// ---|---
ConstTensor6View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(b);
  return ConstTensor6View(
      mdata + OFFSET(b), mlr, mvr, msr, mpr, mrr, mcr, l, v, s, p, r, c);
}
// |--|---
ConstTensor5View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(b);
  return ConstTensor5View(
      mdata + OFFSET(l) + OFFSET(b), mvr, msr, mpr, mrr, mcr, v, s, p, r, c);
}

// --|----
ConstTensor6View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  return ConstTensor6View(
      mdata + OFFSET(s), mlr, mvr, mbr, mpr, mrr, mcr, l, v, b, p, r, c);
}
// |-|----
ConstTensor5View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  return ConstTensor5View(
      mdata + OFFSET(l) + OFFSET(s), mvr, mbr, mpr, mrr, mcr, v, b, p, r, c);
}

// -|-----
ConstTensor6View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  return ConstTensor6View(
      mdata + OFFSET(v), mlr, msr, mbr, mpr, mrr, mcr, l, s, b, p, r, c);
}
// ||-----
ConstTensor5View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  return ConstTensor5View(
      mdata + OFFSET(l) + OFFSET(v), msr, mbr, mpr, mrr, mcr, s, b, p, r, c);
}

// -----||
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(r);
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(r) + OFFSET(c), mlr, mvr, msr, mbr, mpr, l, v, s, b, p);
}
// |----||
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(l);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(r) + OFFSET(c),
                          mvr,
                          msr,
                          mbr,
                          mpr,
                          v,
                          s,
                          b,
                          p);
}

// ----|-|
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(p);
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(p) + OFFSET(c), mlr, mvr, msr, mbr, mrr, l, v, s, b, r);
}
// |---|-|
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(c),
                          mvr,
                          msr,
                          mbr,
                          mrr,
                          v,
                          s,
                          b,
                          r);
}

// ---|--|
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(b);
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(b) + OFFSET(c), mlr, mvr, msr, mpr, mrr, l, v, s, p, r);
}
// |--|--|
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(c),
                          mvr,
                          msr,
                          mpr,
                          mrr,
                          v,
                          s,
                          p,
                          r);
}

// --|---|
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(s) + OFFSET(c), mlr, mvr, mbr, mpr, mrr, l, v, b, p, r);
}
// |-|---|
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(c),
                          mvr,
                          mbr,
                          mpr,
                          mrr,
                          v,
                          b,
                          p,
                          r);
}

// -|----|
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(v) + OFFSET(c), mlr, msr, mbr, mpr, mrr, l, s, b, p, r);
}
// ||----|
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(c),
                          msr,
                          mbr,
                          mpr,
                          mrr,
                          s,
                          b,
                          p,
                          r);
}

// ----||-
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(p);
  CHECK(r);
  return ConstTensor5View(
      mdata + OFFSET(p) + OFFSET(r), mlr, mvr, msr, mbr, mcr, l, v, s, b, c);
}
// |---||-
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r),
                          mvr,
                          msr,
                          mbr,
                          mcr,
                          v,
                          s,
                          b,
                          c);
}

// ---|-|-
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(b);
  CHECK(r);
  return ConstTensor5View(
      mdata + OFFSET(b) + OFFSET(r), mlr, mvr, msr, mpr, mcr, l, v, s, p, c);
}
// |--|-|-
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r),
                          mvr,
                          msr,
                          mpr,
                          mcr,
                          v,
                          s,
                          p,
                          c);
}

// --|--|-
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(r);
  return ConstTensor5View(
      mdata + OFFSET(s) + OFFSET(r), mlr, mvr, mbr, mpr, mcr, l, v, b, p, c);
}
// |-|--|-
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r),
                          mvr,
                          mbr,
                          mpr,
                          mcr,
                          v,
                          b,
                          p,
                          c);
}

// -|---|-
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(r);
  return ConstTensor5View(
      mdata + OFFSET(v) + OFFSET(r), mlr, msr, mbr, mpr, mcr, l, s, b, p, c);
}
// ||---|-
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r),
                          msr,
                          mbr,
                          mpr,
                          mcr,
                          s,
                          b,
                          p,
                          c);
}

// ---||--
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(b);
  CHECK(p);
  return ConstTensor5View(
      mdata + OFFSET(b) + OFFSET(p), mlr, mvr, msr, mrr, mcr, l, v, s, r, c);
}
// |--||--
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p),
                          mvr,
                          msr,
                          mrr,
                          mcr,
                          v,
                          s,
                          r,
                          c);
}

// --|-|--
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(p);
  return ConstTensor5View(
      mdata + OFFSET(s) + OFFSET(p), mlr, mvr, mbr, mrr, mcr, l, v, b, r, c);
}
// |-|-|--
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p),
                          mvr,
                          mbr,
                          mrr,
                          mcr,
                          v,
                          b,
                          r,
                          c);
}

// -|--|--
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(p);
  return ConstTensor5View(
      mdata + OFFSET(v) + OFFSET(p), mlr, msr, mbr, mrr, mcr, l, s, b, r, c);
}
// ||--|--
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p),
                          msr,
                          mbr,
                          mrr,
                          mcr,
                          s,
                          b,
                          r,
                          c);
}

// --||---
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  return ConstTensor5View(
      mdata + OFFSET(s) + OFFSET(b), mlr, mvr, mpr, mrr, mcr, l, v, p, r, c);
}
// |-||---
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b),
                          mvr,
                          mpr,
                          mrr,
                          mcr,
                          v,
                          p,
                          r,
                          c);
}

// -|-|---
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  return ConstTensor5View(
      mdata + OFFSET(v) + OFFSET(b), mlr, msr, mpr, mrr, mcr, l, s, p, r, c);
}
// ||-|---
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b),
                          msr,
                          mpr,
                          mrr,
                          mcr,
                          s,
                          p,
                          r,
                          c);
}

// -||----
ConstTensor5View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  return ConstTensor5View(
      mdata + OFFSET(v) + OFFSET(s), mlr, mbr, mpr, mrr, mcr, l, b, p, r, c);
}
// |||----
ConstTensor4View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  return ConstTensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s),
                          mbr,
                          mpr,
                          mrr,
                          mcr,
                          b,
                          p,
                          r,
                          c);
}

// ----|||
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              Index c) const {
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mvr,
                          msr,
                          mbr,
                          l,
                          v,
                          s,
                          b);
}
// |---|||
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              Index c) const {
  CHECK(l);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mvr,
                          msr,
                          mbr,
                          v,
                          s,
                          b);
}

// ---|-||
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mvr,
                          msr,
                          mpr,
                          l,
                          v,
                          s,
                          p);
}
// |--|-||
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(l);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mvr,
                          msr,
                          mpr,
                          v,
                          s,
                          p);
}

// --|--||
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mvr,
                          mbr,
                          mpr,
                          l,
                          v,
                          b,
                          p);
}
// |-|--||
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mvr,
                          mbr,
                          mpr,
                          v,
                          b,
                          p);
}

// -|---||
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(r) + OFFSET(c),
                          mlr,
                          msr,
                          mbr,
                          mpr,
                          l,
                          s,
                          b,
                          p);
}
// ||---||
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r) + OFFSET(c),
                          msr,
                          mbr,
                          mpr,
                          s,
                          b,
                          p);
}

// ---||-|
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mlr,
                          mvr,
                          msr,
                          mrr,
                          l,
                          v,
                          s,
                          r);
}
// |--||-|
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mvr,
                          msr,
                          mrr,
                          v,
                          s,
                          r);
}

// --|-|-|
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mlr,
                          mvr,
                          mbr,
                          mrr,
                          l,
                          v,
                          b,
                          r);
}
// |-|-|-|
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mvr,
                          mbr,
                          mrr,
                          v,
                          b,
                          r);
}

// -|--|-|
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(c),
                          mlr,
                          msr,
                          mbr,
                          mrr,
                          l,
                          s,
                          b,
                          r);
}
// ||--|-|
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(c),
                          msr,
                          mbr,
                          mrr,
                          s,
                          b,
                          r);
}

// --||--|
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mlr,
                          mvr,
                          mpr,
                          mrr,
                          l,
                          v,
                          p,
                          r);
}
// |-||--|
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mvr,
                          mpr,
                          mrr,
                          v,
                          p,
                          r);
}

// -|-|--|
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(c),
                          mlr,
                          msr,
                          mpr,
                          mrr,
                          l,
                          s,
                          p,
                          r);
}
// ||-|--|
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(c),
                          msr,
                          mpr,
                          mrr,
                          s,
                          p,
                          r);
}

// -||---|
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(c),
                          mlr,
                          mbr,
                          mpr,
                          mrr,
                          l,
                          b,
                          p,
                          r);
}
// |||---|
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(c),
                          mbr,
                          mpr,
                          mrr,
                          b,
                          p,
                          r);
}

// ---|||-
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mlr,
                          mvr,
                          msr,
                          mcr,
                          l,
                          v,
                          s,
                          c);
}
// |--|||-
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mvr,
                          msr,
                          mcr,
                          v,
                          s,
                          c);
}

// --|-||-
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mlr,
                          mvr,
                          mbr,
                          mcr,
                          l,
                          v,
                          b,
                          c);
}
// |-|-||-
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mvr,
                          mbr,
                          mcr,
                          v,
                          b,
                          c);
}

// -|--||-
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r),
                          mlr,
                          msr,
                          mbr,
                          mcr,
                          l,
                          s,
                          b,
                          c);
}
// ||--||-
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r),
                          msr,
                          mbr,
                          mcr,
                          s,
                          b,
                          c);
}

// --||-|-
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mlr,
                          mvr,
                          mpr,
                          mcr,
                          l,
                          v,
                          p,
                          c);
}
// |-||-|-
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mvr,
                          mpr,
                          mcr,
                          v,
                          p,
                          c);
}

// -|-|-|-
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r),
                          mlr,
                          msr,
                          mpr,
                          mcr,
                          l,
                          s,
                          p,
                          c);
}
// ||-|-|-
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r),
                          msr,
                          mpr,
                          mcr,
                          s,
                          p,
                          c);
}

// -||--|-
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r),
                          mlr,
                          mbr,
                          mpr,
                          mcr,
                          l,
                          b,
                          p,
                          c);
}
// |||--|-
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r),
                          mbr,
                          mpr,
                          mcr,
                          b,
                          p,
                          c);
}

// --|||--
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mlr,
                          mvr,
                          mrr,
                          mcr,
                          l,
                          v,
                          r,
                          c);
}
// |-|||--
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mvr,
                          mrr,
                          mcr,
                          v,
                          r,
                          c);
}

// -|-||--
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p),
                          mlr,
                          msr,
                          mrr,
                          mcr,
                          l,
                          s,
                          r,
                          c);
}
// ||-||--
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p),
                          msr,
                          mrr,
                          mcr,
                          s,
                          r,
                          c);
}

// -||-|--
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p),
                          mlr,
                          mbr,
                          mrr,
                          mcr,
                          l,
                          b,
                          r,
                          c);
}
// |||-|--
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p),
                          mbr,
                          mrr,
                          mcr,
                          b,
                          r,
                          c);
}

// -|||---
ConstTensor4View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return ConstTensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b),
                          mlr,
                          mpr,
                          mrr,
                          mcr,
                          l,
                          p,
                          r,
                          c);
}
// ||||---
ConstTensor3View ConstTensor7View::operator()(Index l,
                                              Index v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return ConstTensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b),
                          mpr,
                          mrr,
                          mcr,
                          p,
                          r,
                          c);
}

// -||||--
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                          mlr,
                          mrr,
                          mcr,
                          l,
                          r,
                          c);
}
// |||||--
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             const Range& r,
                                             const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
      mrr,
      mcr,
      r,
      c);
}

// -|||-|-
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                          mlr,
                          mpr,
                          mcr,
                          l,
                          p,
                          c);
}
// ||||-|-
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             const Range& p,
                                             Index r,
                                             const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
      mpr,
      mcr,
      p,
      c);
}

// -||-||-
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                          mlr,
                          mbr,
                          mcr,
                          l,
                          b,
                          c);
}
// |||-||-
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             const Range& b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
      mbr,
      mcr,
      b,
      c);
}

// -|-|||-
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mlr,
                          msr,
                          mcr,
                          l,
                          s,
                          c);
}
// ||-|||-
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             const Range& s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      msr,
      mcr,
      s,
      c);
}

// --||||-
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                          mlr,
                          mvr,
                          mcr,
                          l,
                          v,
                          c);
}
// |-||||-
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             const Range& v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      mvr,
      mcr,
      v,
      c);
}

// -|||--|
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                          mlr,
                          mpr,
                          mrr,
                          l,
                          p,
                          r);
}
// ||||--|
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             const Range& p,
                                             const Range& r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
      mpr,
      mrr,
      p,
      r);
}

// -||-|-|
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                          mlr,
                          mbr,
                          mrr,
                          l,
                          b,
                          r);
}
// |||-|-|
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             const Range& b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
      mbr,
      mrr,
      b,
      r);
}

// -|-||-|
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mlr,
                          msr,
                          mrr,
                          l,
                          s,
                          r);
}
// ||-||-|
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             const Range& s,
                                             Index b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      msr,
      mrr,
      s,
      r);
}

// --|||-|
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                          mlr,
                          mvr,
                          mrr,
                          l,
                          v,
                          r);
}
// |-|||-|
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             const Range& v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      mvr,
      mrr,
      v,
      r);
}

// -||--||
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mbr,
                          mpr,
                          l,
                          b,
                          p);
}
// |||--||
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             const Range& b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
      mbr,
      mpr,
      b,
      p);
}

// -|-|-||
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mlr,
                          msr,
                          mpr,
                          l,
                          s,
                          p);
}
// ||-|-||
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             const Range& s,
                                             Index b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      msr,
      mpr,
      s,
      p);
}

// --||-||
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mvr,
                          mpr,
                          l,
                          v,
                          p);
}
// |-||-||
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             const Range& v,
                                             Index s,
                                             Index b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      mvr,
      mpr,
      v,
      p);
}

// -|--|||
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              Index c) const {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr,
                          msr,
                          mbr,
                          l,
                          s,
                          b);
}
// ||--|||
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             const Range& s,
                                             const Range& b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      msr,
      mbr,
      s,
      b);
}

// --|-|||
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              Index c) const {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mvr,
                          mbr,
                          l,
                          v,
                          b);
}
// |-|-|||
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             const Range& v,
                                             Index s,
                                             const Range& b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mvr,
      mbr,
      v,
      b);
}

// ---||||
ConstTensor3View ConstTensor7View::operator()(const Range& l,
                                              const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              Index r,
                                              Index c) const {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                          mlr,
                          mvr,
                          msr,
                          l,
                          v,
                          s);
}
// |--||||
ConstMatrixView ConstTensor7View::operator()(Index l,
                                             const Range& v,
                                             const Range& s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mvr,
      msr,
      v,
      s);
}

// -|||||-
ConstMatrixView ConstTensor7View::operator()(const Range& l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      mlr,
      mcr,
      l,
      c);
}
// ||||||-
ConstVectorView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             const Range& c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) +
                             OFFSET(p) + OFFSET(r),
                         mcr,
                         c);
}

// -||||-|
ConstMatrixView ConstTensor7View::operator()(const Range& l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      mlr,
      mrr,
      l,
      r);
}
// |||||-|
ConstVectorView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             const Range& r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) +
                             OFFSET(p) + OFFSET(c),
                         mrr,
                         r);
}

// -|||-||
ConstMatrixView ConstTensor7View::operator()(const Range& l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      mlr,
      mpr,
      l,
      p);
}
// ||||-||
ConstVectorView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             const Range& p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) +
                             OFFSET(r) + OFFSET(c),
                         mpr,
                         p);
}

// -||-|||
ConstMatrixView ConstTensor7View::operator()(const Range& l,
                                             Index v,
                                             Index s,
                                             const Range& b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mlr,
      mbr,
      l,
      b);
}
// |||-|||
ConstVectorView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             Index s,
                                             const Range& b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) +
                             OFFSET(r) + OFFSET(c),
                         mbr,
                         b);
}

// -|-||||
ConstMatrixView ConstTensor7View::operator()(const Range& l,
                                             Index v,
                                             const Range& s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mlr,
      msr,
      l,
      s);
}
// ||-||||
ConstVectorView ConstTensor7View::operator()(Index l,
                                             Index v,
                                             const Range& s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) +
                             OFFSET(r) + OFFSET(c),
                         msr,
                         s);
}

// --|||||
ConstMatrixView ConstTensor7View::operator()(const Range& l,
                                             const Range& v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mlr,
      mvr,
      l,
      v);
}
// |-|||||
ConstVectorView ConstTensor7View::operator()(Index l,
                                             const Range& v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) +
                             OFFSET(r) + OFFSET(c),
                         mvr,
                         v);
}

// -||||||
ConstVectorView ConstTensor7View::operator()(const Range& l,
                                             Index v,
                                             Index s,
                                             Index b,
                                             Index p,
                                             Index r,
                                             Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) +
                             OFFSET(r) + OFFSET(c),
                         mlr,
                         l);
}

/** Return const iterator to first sub-tensor. */
ConstIterator7D ConstTensor7View::begin() const {
  return ConstIterator7D(
      ConstTensor6View(mdata + mlr.mstart, mvr, msr, mbr, mpr, mrr, mcr),
      mlr.mstride);
}

/** Return const iterator behind last sub-tensor. */
ConstIterator7D ConstTensor7View::end() const {
  return ConstIterator7D(
      ConstTensor6View(mdata + mlr.mstart + (mlr.mextent) * mlr.mstride,
                       mvr,
                       msr,
                       mbr,
                       mpr,
                       mrr,
                       mcr),
      mlr.mstride);
}

/** Special constructor to make a Tensor7 view of a Tensor6. */
ConstTensor7View::ConstTensor7View(const ConstTensor6View& a)
    : mlr(0,
          1,
          a.mvr.mextent * a.msr.mextent * a.mbr.mextent * a.mpr.mextent *
              a.mrr.mextent * a.mcr.mextent),
      mvr(a.mvr),
      msr(a.msr),
      mbr(a.mbr),
      mpr(a.mpr),
      mrr(a.mrr),
      mcr(a.mcr),
      mdata(a.mdata) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor7 to initialize
    its own Tensor7View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
ConstTensor7View::ConstTensor7View(Numeric* data,
                                   const Range& l,
                                   const Range& v,
                                   const Range& s,
                                   const Range& b,
                                   const Range& p,
                                   const Range& r,
                                   const Range& c)
    : mlr(l), mvr(v), msr(s), mbr(b), mpr(p), mrr(r), mcr(c), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
ConstTensor7View::ConstTensor7View(Numeric* data,
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
                                   const Range& nc)
    : mlr(pl, nl),
      mvr(pv, nv),
      msr(ps, ns),
      mbr(pb, nb),
      mpr(pp, np),
      mrr(pr, nr),
      mcr(pc, nc),
      mdata(data) {
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Tensor6 to print each page in turn. */
std::ostream& operator<<(std::ostream& os, const ConstTensor7View& v) {
  // Page iterators:
  ConstIterator7D ip = v.begin();
  const ConstIterator7D end_page = v.end();

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

// Functions for Tensor7View:
// -------------------------

// Non-const index operators:

// -------
Tensor7View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  return Tensor7View(
      mdata, mlr, mvr, msr, mbr, mpr, mrr, mcr, l, v, s, b, p, r, c);
}
// |------
Tensor6View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  return Tensor6View(
      mdata + OFFSET(l), mvr, msr, mbr, mpr, mrr, mcr, v, s, b, p, r, c);
}

// ------|
Tensor6View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(c);
  return Tensor6View(
      mdata + OFFSET(c), mlr, mvr, msr, mbr, mpr, mrr, l, v, s, b, p, r);
}
// |-----|
Tensor5View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(c);
  return Tensor5View(
      mdata + OFFSET(l) + OFFSET(c), mvr, msr, mbr, mpr, mrr, v, s, b, p, r);
}

// -----|-
Tensor6View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(r);
  return Tensor6View(
      mdata + OFFSET(r), mlr, mvr, msr, mbr, mpr, mcr, l, v, s, b, p, c);
}
// |----|-
Tensor5View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(r);
  return Tensor5View(
      mdata + OFFSET(l) + OFFSET(r), mvr, msr, mbr, mpr, mcr, v, s, b, p, c);
}

// ----|--
Tensor6View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(p);
  return Tensor6View(
      mdata + OFFSET(p), mlr, mvr, msr, mbr, mrr, mcr, l, v, s, b, r, c);
}
// |---|--
Tensor5View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(p);
  return Tensor5View(
      mdata + OFFSET(l) + OFFSET(p), mvr, msr, mbr, mrr, mcr, v, s, b, r, c);
}

// ---|---
Tensor6View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(b);
  return Tensor6View(
      mdata + OFFSET(b), mlr, mvr, msr, mpr, mrr, mcr, l, v, s, p, r, c);
}
// |--|---
Tensor5View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(b);
  return Tensor5View(
      mdata + OFFSET(l) + OFFSET(b), mvr, msr, mpr, mrr, mcr, v, s, p, r, c);
}

// --|----
Tensor6View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  return Tensor6View(
      mdata + OFFSET(s), mlr, mvr, mbr, mpr, mrr, mcr, l, v, b, p, r, c);
}
// |-|----
Tensor5View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  return Tensor5View(
      mdata + OFFSET(l) + OFFSET(s), mvr, mbr, mpr, mrr, mcr, v, b, p, r, c);
}

// -|-----
Tensor6View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  return Tensor6View(
      mdata + OFFSET(v), mlr, msr, mbr, mpr, mrr, mcr, l, s, b, p, r, c);
}
// ||-----
Tensor5View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  return Tensor5View(
      mdata + OFFSET(l) + OFFSET(v), msr, mbr, mpr, mrr, mcr, s, b, p, r, c);
}

// -----||
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(r);
  CHECK(c);
  return Tensor5View(
      mdata + OFFSET(r) + OFFSET(c), mlr, mvr, msr, mbr, mpr, l, v, s, b, p);
}
// |----||
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(l);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(r) + OFFSET(c),
                     mvr,
                     msr,
                     mbr,
                     mpr,
                     v,
                     s,
                     b,
                     p);
}

// ----|-|
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(p);
  CHECK(c);
  return Tensor5View(
      mdata + OFFSET(p) + OFFSET(c), mlr, mvr, msr, mbr, mrr, l, v, s, b, r);
}
// |---|-|
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(c),
                     mvr,
                     msr,
                     mbr,
                     mrr,
                     v,
                     s,
                     b,
                     r);
}

// ---|--|
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(b);
  CHECK(c);
  return Tensor5View(
      mdata + OFFSET(b) + OFFSET(c), mlr, mvr, msr, mpr, mrr, l, v, s, p, r);
}
// |--|--|
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(c),
                     mvr,
                     msr,
                     mpr,
                     mrr,
                     v,
                     s,
                     p,
                     r);
}

// --|---|
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(s);
  CHECK(c);
  return Tensor5View(
      mdata + OFFSET(s) + OFFSET(c), mlr, mvr, mbr, mpr, mrr, l, v, b, p, r);
}
// |-|---|
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(c),
                     mvr,
                     mbr,
                     mpr,
                     mrr,
                     v,
                     b,
                     p,
                     r);
}

// -|----|
Tensor5View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(c);
  return Tensor5View(
      mdata + OFFSET(v) + OFFSET(c), mlr, msr, mbr, mpr, mrr, l, s, b, p, r);
}
// ||----|
Tensor4View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(c),
                     msr,
                     mbr,
                     mpr,
                     mrr,
                     s,
                     b,
                     p,
                     r);
}

// ----||-
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(p);
  CHECK(r);
  return Tensor5View(
      mdata + OFFSET(p) + OFFSET(r), mlr, mvr, msr, mbr, mcr, l, v, s, b, c);
}
// |---||-
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r),
                     mvr,
                     msr,
                     mbr,
                     mcr,
                     v,
                     s,
                     b,
                     c);
}

// ---|-|-
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(b);
  CHECK(r);
  return Tensor5View(
      mdata + OFFSET(b) + OFFSET(r), mlr, mvr, msr, mpr, mcr, l, v, s, p, c);
}
// |--|-|-
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r),
                     mvr,
                     msr,
                     mpr,
                     mcr,
                     v,
                     s,
                     p,
                     c);
}

// --|--|-
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(s);
  CHECK(r);
  return Tensor5View(
      mdata + OFFSET(s) + OFFSET(r), mlr, mvr, mbr, mpr, mcr, l, v, b, p, c);
}
// |-|--|-
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r),
                     mvr,
                     mbr,
                     mpr,
                     mcr,
                     v,
                     b,
                     p,
                     c);
}

// -|---|-
Tensor5View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(r);
  return Tensor5View(
      mdata + OFFSET(v) + OFFSET(r), mlr, msr, mbr, mpr, mcr, l, s, b, p, c);
}
// ||---|-
Tensor4View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r),
                     msr,
                     mbr,
                     mpr,
                     mcr,
                     s,
                     b,
                     p,
                     c);
}

// ---||--
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(b);
  CHECK(p);
  return Tensor5View(
      mdata + OFFSET(b) + OFFSET(p), mlr, mvr, msr, mrr, mcr, l, v, s, r, c);
}
// |--||--
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p),
                     mvr,
                     msr,
                     mrr,
                     mcr,
                     v,
                     s,
                     r,
                     c);
}

// --|-|--
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  CHECK(p);
  return Tensor5View(
      mdata + OFFSET(s) + OFFSET(p), mlr, mvr, mbr, mrr, mcr, l, v, b, r, c);
}
// |-|-|--
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p),
                     mvr,
                     mbr,
                     mrr,
                     mcr,
                     v,
                     b,
                     r,
                     c);
}

// -|--|--
Tensor5View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(p);
  return Tensor5View(
      mdata + OFFSET(v) + OFFSET(p), mlr, msr, mbr, mrr, mcr, l, s, b, r, c);
}
// ||--|--
Tensor4View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p),
                     msr,
                     mbr,
                     mrr,
                     mcr,
                     s,
                     b,
                     r,
                     c);
}

// --||---
Tensor5View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  CHECK(b);
  return Tensor5View(
      mdata + OFFSET(s) + OFFSET(b), mlr, mvr, mpr, mrr, mcr, l, v, p, r, c);
}
// |-||---
Tensor4View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b),
                     mvr,
                     mpr,
                     mrr,
                     mcr,
                     v,
                     p,
                     r,
                     c);
}

// -|-|---
Tensor5View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(b);
  return Tensor5View(
      mdata + OFFSET(v) + OFFSET(b), mlr, msr, mpr, mrr, mcr, l, s, p, r, c);
}
// ||-|---
Tensor4View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b),
                     msr,
                     mpr,
                     mrr,
                     mcr,
                     s,
                     p,
                     r,
                     c);
}

// -||----
Tensor5View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  return Tensor5View(
      mdata + OFFSET(v) + OFFSET(s), mlr, mbr, mpr, mrr, mcr, l, b, p, r, c);
}
// |||----
Tensor4View Tensor7View::operator()(Index l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  return Tensor4View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s),
                     mbr,
                     mpr,
                     mrr,
                     mcr,
                     b,
                     p,
                     r,
                     c);
}

// ----|||
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    Index c) {
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(p) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mvr,
                     msr,
                     mbr,
                     l,
                     v,
                     s,
                     b);
}
// |---|||
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    Index c) {
  CHECK(l);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                     mvr,
                     msr,
                     mbr,
                     v,
                     s,
                     b);
}

// ---|-||
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mvr,
                     msr,
                     mpr,
                     l,
                     v,
                     s,
                     p);
}
// |--|-||
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(l);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                     mvr,
                     msr,
                     mpr,
                     v,
                     s,
                     p);
}

// --|--||
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mvr,
                     mbr,
                     mpr,
                     l,
                     v,
                     b,
                     p);
}
// |-|--||
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                     mvr,
                     mbr,
                     mpr,
                     v,
                     b,
                     p);
}

// -|---||
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(r) + OFFSET(c),
                     mlr,
                     msr,
                     mbr,
                     mpr,
                     l,
                     s,
                     b,
                     p);
}
// ||---||
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(r) + OFFSET(c),
                     msr,
                     mbr,
                     mpr,
                     s,
                     b,
                     p);
}

// ---||-|
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(c),
                     mlr,
                     mvr,
                     msr,
                     mrr,
                     l,
                     v,
                     s,
                     r);
}
// |--||-|
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                     mvr,
                     msr,
                     mrr,
                     v,
                     s,
                     r);
}

// --|-|-|
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(c),
                     mlr,
                     mvr,
                     mbr,
                     mrr,
                     l,
                     v,
                     b,
                     r);
}
// |-|-|-|
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                     mvr,
                     mbr,
                     mrr,
                     v,
                     b,
                     r);
}

// -|--|-|
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(c),
                     mlr,
                     msr,
                     mbr,
                     mrr,
                     l,
                     s,
                     b,
                     r);
}
// ||--|-|
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(c),
                     msr,
                     mbr,
                     mrr,
                     s,
                     b,
                     r);
}

// --||--|
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(c),
                     mlr,
                     mvr,
                     mpr,
                     mrr,
                     l,
                     v,
                     p,
                     r);
}
// |-||--|
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                     mvr,
                     mpr,
                     mrr,
                     v,
                     p,
                     r);
}

// -|-|--|
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(c),
                     mlr,
                     msr,
                     mpr,
                     mrr,
                     l,
                     s,
                     p,
                     r);
}
// ||-|--|
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(c),
                     msr,
                     mpr,
                     mrr,
                     s,
                     p,
                     r);
}

// -||---|
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(c),
                     mlr,
                     mbr,
                     mpr,
                     mrr,
                     l,
                     b,
                     p,
                     r);
}
// |||---|
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(c),
                     mbr,
                     mpr,
                     mrr,
                     b,
                     p,
                     r);
}

// ---|||-
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r),
                     mlr,
                     mvr,
                     msr,
                     mcr,
                     l,
                     v,
                     s,
                     c);
}
// |--|||-
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                     mvr,
                     msr,
                     mcr,
                     v,
                     s,
                     c);
}

// --|-||-
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r),
                     mlr,
                     mvr,
                     mbr,
                     mcr,
                     l,
                     v,
                     b,
                     c);
}
// |-|-||-
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                     mvr,
                     mbr,
                     mcr,
                     v,
                     b,
                     c);
}

// -|--||-
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r),
                     mlr,
                     msr,
                     mbr,
                     mcr,
                     l,
                     s,
                     b,
                     c);
}
// ||--||-
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r),
                     msr,
                     mbr,
                     mcr,
                     s,
                     b,
                     c);
}

// --||-|-
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r),
                     mlr,
                     mvr,
                     mpr,
                     mcr,
                     l,
                     v,
                     p,
                     c);
}
// |-||-|-
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                     mvr,
                     mpr,
                     mcr,
                     v,
                     p,
                     c);
}

// -|-|-|-
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r),
                     mlr,
                     msr,
                     mpr,
                     mcr,
                     l,
                     s,
                     p,
                     c);
}
// ||-|-|-
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r),
                     msr,
                     mpr,
                     mcr,
                     s,
                     p,
                     c);
}

// -||--|-
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r),
                     mlr,
                     mbr,
                     mpr,
                     mcr,
                     l,
                     b,
                     p,
                     c);
}
// |||--|-
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r),
                     mbr,
                     mpr,
                     mcr,
                     b,
                     p,
                     c);
}

// --|||--
Tensor4View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p),
                     mlr,
                     mvr,
                     mrr,
                     mcr,
                     l,
                     v,
                     r,
                     c);
}
// |-|||--
Tensor3View Tensor7View::operator()(Index l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                     mvr,
                     mrr,
                     mcr,
                     v,
                     r,
                     c);
}

// -|-||--
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p),
                     mlr,
                     msr,
                     mrr,
                     mcr,
                     l,
                     s,
                     r,
                     c);
}
// ||-||--
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p),
                     msr,
                     mrr,
                     mcr,
                     s,
                     r,
                     c);
}

// -||-|--
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p),
                     mlr,
                     mbr,
                     mrr,
                     mcr,
                     l,
                     b,
                     r,
                     c);
}
// |||-|--
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p),
                     mbr,
                     mrr,
                     mcr,
                     b,
                     r,
                     c);
}

// -|||---
Tensor4View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return Tensor4View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b),
                     mlr,
                     mpr,
                     mrr,
                     mcr,
                     l,
                     p,
                     r,
                     c);
}
// ||||---
Tensor3View Tensor7View::operator()(Index l,
                                    Index v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return Tensor3View(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b),
                     mpr,
                     mrr,
                     mcr,
                     p,
                     r,
                     c);
}

// -||||--
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
                     mlr,
                     mrr,
                     mcr,
                     l,
                     r,
                     c);
}
// |||||--
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   const Range& r,
                                   const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p),
      mrr,
      mcr,
      r,
      c);
}

// -|||-|-
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
                     mlr,
                     mpr,
                     mcr,
                     l,
                     p,
                     c);
}
// ||||-|-
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   const Range& p,
                                   Index r,
                                   const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r),
      mpr,
      mcr,
      p,
      c);
}

// -||-||-
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
                     mlr,
                     mbr,
                     mcr,
                     l,
                     b,
                     c);
}
// |||-||-
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   Index s,
                                   const Range& b,
                                   Index p,
                                   Index r,
                                   const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r),
      mbr,
      mcr,
      b,
      c);
}

// -|-|||-
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                     mlr,
                     msr,
                     mcr,
                     l,
                     s,
                     c);
}
// ||-|||-
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   const Range& s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      msr,
      mcr,
      s,
      c);
}

// --||||-
Tensor3View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
                     mlr,
                     mvr,
                     mcr,
                     l,
                     v,
                     c);
}
// |-||||-
MatrixView Tensor7View::operator()(Index l,
                                   const Range& v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   const Range& c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      mvr,
      mcr,
      v,
      c);
}

// -|||--|
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
                     mlr,
                     mpr,
                     mrr,
                     l,
                     p,
                     r);
}
// ||||--|
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   const Range& p,
                                   const Range& r,
                                   Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c),
      mpr,
      mrr,
      p,
      r);
}

// -||-|-|
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
                     mlr,
                     mbr,
                     mrr,
                     l,
                     b,
                     r);
}
// |||-|-|
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   Index s,
                                   const Range& b,
                                   Index p,
                                   const Range& r,
                                   Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c),
      mbr,
      mrr,
      b,
      r);
}

// -|-||-|
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                     mlr,
                     msr,
                     mrr,
                     l,
                     s,
                     r);
}
// ||-||-|
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   const Range& s,
                                   Index b,
                                   Index p,
                                   const Range& r,
                                   Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      msr,
      mrr,
      s,
      r);
}

// --|||-|
Tensor3View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
                     mlr,
                     mvr,
                     mrr,
                     l,
                     v,
                     r);
}
// |-|||-|
MatrixView Tensor7View::operator()(Index l,
                                   const Range& v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   const Range& r,
                                   Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      mvr,
      mrr,
      v,
      r);
}

// -||--||
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mbr,
                     mpr,
                     l,
                     b,
                     p);
}
// |||--||
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   Index s,
                                   const Range& b,
                                   const Range& p,
                                   Index r,
                                   Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c),
      mbr,
      mpr,
      b,
      p);
}

// -|-|-||
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                     mlr,
                     msr,
                     mpr,
                     l,
                     s,
                     p);
}
// ||-|-||
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   const Range& s,
                                   Index b,
                                   const Range& p,
                                   Index r,
                                   Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      msr,
      mpr,
      s,
      p);
}

// --||-||
Tensor3View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mvr,
                     mpr,
                     l,
                     v,
                     p);
}
// |-||-||
MatrixView Tensor7View::operator()(Index l,
                                   const Range& v,
                                   Index s,
                                   Index b,
                                   const Range& p,
                                   Index r,
                                   Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      mvr,
      mpr,
      v,
      p);
}

// -|--|||
Tensor3View Tensor7View::operator()(const Range& l,
                                    Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    Index c) {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                     mlr,
                     msr,
                     mbr,
                     l,
                     s,
                     b);
}
// ||--|||
MatrixView Tensor7View::operator()(Index l,
                                   Index v,
                                   const Range& s,
                                   const Range& b,
                                   Index p,
                                   Index r,
                                   Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      msr,
      mbr,
      s,
      b);
}

// --|-|||
Tensor3View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    Index c) {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mvr,
                     mbr,
                     l,
                     v,
                     b);
}
// |-|-|||
MatrixView Tensor7View::operator()(Index l,
                                   const Range& v,
                                   Index s,
                                   const Range& b,
                                   Index p,
                                   Index r,
                                   Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mvr,
      mbr,
      v,
      b);
}

// ---||||
Tensor3View Tensor7View::operator()(const Range& l,
                                    const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    Index r,
                                    Index c) {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
                     mlr,
                     mvr,
                     msr,
                     l,
                     v,
                     s);
}
// |--||||
MatrixView Tensor7View::operator()(Index l,
                                   const Range& v,
                                   const Range& s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   Index c) {
  CHECK(l);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(l) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mvr,
      msr,
      v,
      s);
}

// -|||||-
MatrixView Tensor7View::operator()(const Range& l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      mlr,
      mcr,
      l,
      c);
}
// ||||||-
VectorView Tensor7View::operator()(
    Index l, Index v, Index s, Index b, Index p, Index r, const Range& c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) +
                        OFFSET(p) + OFFSET(r),
                    mcr,
                    c);
}

// -||||-|
MatrixView Tensor7View::operator()(const Range& l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   const Range& r,
                                   Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      mlr,
      mrr,
      l,
      r);
}
// |||||-|
VectorView Tensor7View::operator()(
    Index l, Index v, Index s, Index b, Index p, const Range& r, Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) +
                        OFFSET(p) + OFFSET(c),
                    mrr,
                    r);
}

// -|||-||
MatrixView Tensor7View::operator()(const Range& l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   const Range& p,
                                   Index r,
                                   Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      mlr,
      mpr,
      l,
      p);
}
// ||||-||
VectorView Tensor7View::operator()(
    Index l, Index v, Index s, Index b, const Range& p, Index r, Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(b) +
                        OFFSET(r) + OFFSET(c),
                    mpr,
                    p);
}

// -||-|||
MatrixView Tensor7View::operator()(const Range& l,
                                   Index v,
                                   Index s,
                                   const Range& b,
                                   Index p,
                                   Index r,
                                   Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mlr,
      mbr,
      l,
      b);
}
// |||-|||
VectorView Tensor7View::operator()(
    Index l, Index v, Index s, const Range& b, Index p, Index r, Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(s) + OFFSET(p) +
                        OFFSET(r) + OFFSET(c),
                    mbr,
                    b);
}

// -|-||||
MatrixView Tensor7View::operator()(const Range& l,
                                   Index v,
                                   const Range& s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mlr,
      msr,
      l,
      s);
}
// ||-||||
VectorView Tensor7View::operator()(
    Index l, Index v, const Range& s, Index b, Index p, Index r, Index c) {
  CHECK(l);
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(mdata + OFFSET(l) + OFFSET(v) + OFFSET(b) + OFFSET(p) +
                        OFFSET(r) + OFFSET(c),
                    msr,
                    s);
}

// --|||||
MatrixView Tensor7View::operator()(const Range& l,
                                   const Range& v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mlr,
      mvr,
      l,
      v);
}
// |-|||||
VectorView Tensor7View::operator()(
    Index l, const Range& v, Index s, Index b, Index p, Index r, Index c) {
  CHECK(l);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(mdata + OFFSET(l) + OFFSET(s) + OFFSET(b) + OFFSET(p) +
                        OFFSET(r) + OFFSET(c),
                    mvr,
                    v);
}

// -||||||
VectorView Tensor7View::operator()(
    const Range& l, Index v, Index s, Index b, Index p, Index r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) +
                        OFFSET(r) + OFFSET(c),
                    mlr,
                    l);
}

/** Conversion to plain C-array.
    
This function returns a pointer to the raw data. It fails if the
Tensor5View is not pointing to the beginning of a Tensor5 or the stride
is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric* Tensor7View::get_c_array() {
  if (mlr.mstart != 0 ||
      mlr.mstride != mvr.mextent * msr.mextent * mbr.mextent * mpr.mextent *
                         mrr.mextent * mcr.mextent ||
      mvr.mstart != 0 ||
      mvr.mstride !=
          msr.mextent * mbr.mextent * mpr.mextent * mrr.mextent * mcr.mextent ||
      msr.mstart != 0 ||
      msr.mstride != mbr.mextent * mpr.mextent * mrr.mextent * mcr.mextent ||
      mbr.mstart != 0 ||
      mbr.mstride != mpr.mextent * mrr.mextent * mcr.mextent ||
      mpr.mstart != 0 || mpr.mstride != mrr.mextent * mcr.mextent ||
      mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1)
    throw std::runtime_error(
        "A Tensor7View can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** Conversion to plain C-array.

  This function returns a pointer to the raw data. It fails if the
  Tensor5View is not pointing to the beginning of a Tensor5 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric* Tensor7View::get_c_array() const {
  if (mlr.mstart != 0 ||
      mlr.mstride != mvr.mextent * msr.mextent * mbr.mextent * mpr.mextent *
                         mrr.mextent * mcr.mextent ||
      mvr.mstart != 0 ||
      mvr.mstride !=
          msr.mextent * mbr.mextent * mpr.mextent * mrr.mextent * mcr.mextent ||
      msr.mstart != 0 ||
      msr.mstride != mbr.mextent * mpr.mextent * mrr.mextent * mcr.mextent ||
      mbr.mstart != 0 ||
      mbr.mstride != mpr.mextent * mrr.mextent * mcr.mextent ||
      mpr.mstart != 0 || mpr.mstride != mrr.mextent * mcr.mextent ||
      mrr.mstart != 0 || mrr.mstride != mcr.mextent || mcr.mstart != 0 ||
      mcr.mstride != 1)
    throw std::runtime_error(
        "A Tensor7View can only be converted to a plain C-array if it's pointing to a continuous block of data");

  return mdata;
}

/** Return iterator to first sub-tensor. */
Iterator7D Tensor7View::begin() {
  return Iterator7D(
      Tensor6View(mdata + mlr.mstart, mvr, msr, mbr, mpr, mrr, mcr),
      mlr.mstride);
}

/** Return iterator behind last sub-tensor. */
Iterator7D Tensor7View::end() {
  return Iterator7D(
      Tensor6View(mdata + mlr.mstart + (mlr.mextent) * mlr.mstride,
                  mvr,
                  msr,
                  mbr,
                  mpr,
                  mrr,
                  mcr),
      mlr.mstride);
}

/** Assignment operator. This copies the data from another Tensor7View
    to this Tensor7View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor7View by
    setting its range. */
Tensor7View& Tensor7View::operator=(const ConstTensor7View& m) {
  // Check that sizes are compatible:
  assert(mlr.mextent == m.mlr.mextent);
  assert(mvr.mextent == m.mvr.mextent);
  assert(msr.mextent == m.msr.mextent);
  assert(mbr.mextent == m.mbr.mextent);
  assert(mpr.mextent == m.mpr.mextent);
  assert(mrr.mextent == m.mrr.mextent);
  assert(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from Tensor7View to Tensor7View. This is a tricky
    one. The problem is that since Tensor7View is derived from
    ConstTensor7View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
Tensor7View& Tensor7View::operator=(const Tensor7View& m) {
  // Check that sizes are compatible:
  assert(mlr.mextent == m.mlr.mextent);
  assert(mvr.mextent == m.mvr.mextent);
  assert(msr.mextent == m.msr.mextent);
  assert(mbr.mextent == m.mbr.mextent);
  assert(mpr.mextent == m.mpr.mextent);
  assert(mrr.mextent == m.mrr.mextent);
  assert(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a Tensor7. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
Tensor7View& Tensor7View::operator=(const Tensor7& m) {
  // Check that sizes are compatible:
  assert(mlr.mextent == m.mlr.mextent);
  assert(mvr.mextent == m.mvr.mextent);
  assert(msr.mextent == m.msr.mextent);
  assert(mbr.mextent == m.mbr.mextent);
  assert(mpr.mextent == m.mpr.mextent);
  assert(mrr.mextent == m.mrr.mextent);
  assert(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assigning a scalar to a Tensor7View will set all elements to this
    value. */
Tensor7View& Tensor7View::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
Tensor7View& Tensor7View::operator*=(Numeric x) {
  const Iterator7D ep = end();
  for (Iterator7D p = begin(); p != ep; ++p) {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
Tensor7View& Tensor7View::operator/=(Numeric x) {
  const Iterator7D ep = end();
  for (Iterator7D p = begin(); p != ep; ++p) {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
Tensor7View& Tensor7View::operator+=(Numeric x) {
  const Iterator7D ep = end();
  for (Iterator7D p = begin(); p != ep; ++p) {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
Tensor7View& Tensor7View::operator-=(Numeric x) {
  const Iterator7D ep = end();
  for (Iterator7D p = begin(); p != ep; ++p) {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor7. */
Tensor7View& Tensor7View::operator*=(const ConstTensor7View& x) {
  assert(nlibraries() == x.nlibraries());
  assert(nvitrines() == x.nvitrines());
  assert(nshelves() == x.nshelves());
  assert(nbooks() == x.nbooks());
  assert(npages() == x.npages());
  assert(nrows() == x.nrows());
  assert(ncols() == x.ncols());
  ConstIterator7D xp = x.begin();
  Iterator7D p = begin();
  const Iterator7D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p *= *xp;
  }
  return *this;
}

/** Element-vise division by another Tensor7. */
Tensor7View& Tensor7View::operator/=(const ConstTensor7View& x) {
  assert(nlibraries() == x.nlibraries());
  assert(nvitrines() == x.nvitrines());
  assert(nshelves() == x.nshelves());
  assert(nbooks() == x.nbooks());
  assert(npages() == x.npages());
  assert(nrows() == x.nrows());
  assert(ncols() == x.ncols());
  ConstIterator7D xp = x.begin();
  Iterator7D p = begin();
  const Iterator7D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p /= *xp;
  }
  return *this;
}

/** Element-vise addition of another Tensor7. */
Tensor7View& Tensor7View::operator+=(const ConstTensor7View& x) {
  assert(nlibraries() == x.nlibraries());
  assert(nvitrines() == x.nvitrines());
  assert(nshelves() == x.nshelves());
  assert(nbooks() == x.nbooks());
  assert(npages() == x.npages());
  assert(nrows() == x.nrows());
  assert(ncols() == x.ncols());
  ConstIterator7D xp = x.begin();
  Iterator7D p = begin();
  const Iterator7D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p += *xp;
  }
  return *this;
}

/** Element-vise subtraction of another Tensor7. */
Tensor7View& Tensor7View::operator-=(const ConstTensor7View& x) {
  assert(nlibraries() == x.nlibraries());
  assert(nvitrines() == x.nvitrines());
  assert(nshelves() == x.nshelves());
  assert(nbooks() == x.nbooks());
  assert(npages() == x.npages());
  assert(nrows() == x.nrows());
  assert(ncols() == x.ncols());
  ConstIterator7D xp = x.begin();
  Iterator7D p = begin();
  const Iterator7D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p -= *xp;
  }
  return *this;
}

/** Special constructor to make a Tensor7 view of a Tensor6. */
Tensor7View::Tensor7View(const Tensor6View& a)
    : ConstTensor7View(a.mdata,
                       Range(0,
                             1,
                             a.mvr.mextent * a.msr.mextent * a.mbr.mextent *
                                 a.mpr.mextent * a.mrr.mextent * a.mcr.mextent),
                       a.mvr,
                       a.msr,
                       a.mbr,
                       a.mpr,
                       a.mrr,
                       a.mcr) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor7 to initialize its
    own Tensor7View part. */
Tensor7View::Tensor7View(Numeric* data,
                         const Range& l,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c)
    : ConstTensor7View(data, l, v, s, b, p, r, c) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pl Previous range.
    \param pv Previous range.
    \param ps Previous range.
    \param pb Previous range.
    \param pp Previous range.
    \param pr Previous range.
    \param pc Previous range.
    \param nl New Range.
    \param nv New Range.
    \param ns New Range.
    \param nb New Range.
    \param np New Range.
    \param nr New Range.
    \param nc New Range.
  */
Tensor7View::Tensor7View(Numeric* data,
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
                         const Range& nc)
    : ConstTensor7View(
          data, pl, pv, ps, pb, pp, pr, pc, nl, nv, ns, nb, np, nr, nc) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
void copy(ConstIterator7D origin,
          const ConstIterator7D& end,
          Iterator7D target) {
  for (; origin != end; ++origin, ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(origin->begin(), origin->end(), target->begin());
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, Iterator7D target, const Iterator7D& end) {
  for (; target != end; ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(x, target->begin(), target->end());
  }
}

// Functions for Tensor7:
// ---------------------

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
Tensor7::Tensor7(Index l, Index v, Index s, Index b, Index p, Index r, Index c)
    : Tensor7View(new Numeric[l * v * s * b * p * r * c],
                  Range(0, l, v * s * b * p * r * c),
                  Range(0, v, s * b * p * r * c),
                  Range(0, s, b * p * r * c),
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Tensor7::Tensor7(
    Index l, Index v, Index s, Index b, Index p, Index r, Index c, Numeric fill)
    : Tensor7View(new Numeric[l * v * s * b * p * r * c],
                  Range(0, l, v * s * b * p * r * c),
                  Range(0, v, s * b * p * r * c),
                  Range(0, s, b * p * r * c),
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, l * v * s * b * p * r * c, fill);
}

/** Copy constructor from Tensor7View. This automatically sets the size
    and copies the data. */
Tensor7::Tensor7(const ConstTensor7View& m)
    : Tensor7View(
          new Numeric[m.nlibraries() * m.nvitrines() * m.nshelves() *
                      m.nbooks() * m.npages() * m.nrows() * m.ncols()],
          Range(0,
                m.nlibraries(),
                m.nvitrines() * m.nshelves() * m.nbooks() * m.npages() *
                    m.nrows() * m.ncols()),
          Range(0,
                m.nvitrines(),
                m.nshelves() * m.nbooks() * m.npages() * m.nrows() * m.ncols()),
          Range(
              0, m.nshelves(), m.nbooks() * m.npages() * m.nrows() * m.ncols()),
          Range(0, m.nbooks(), m.npages() * m.nrows() * m.ncols()),
          Range(0, m.npages(), m.nrows() * m.ncols()),
          Range(0, m.nrows(), m.ncols()),
          Range(0, m.ncols())) {
  copy(m.begin(), m.end(), begin());
}

/** Copy constructor from Tensor7. This automatically sets the size
    and copies the data. */
Tensor7::Tensor7(const Tensor7& m)
    : Tensor7View(
          new Numeric[m.nlibraries() * m.nvitrines() * m.nshelves() *
                      m.nbooks() * m.npages() * m.nrows() * m.ncols()],
          Range(0,
                m.nlibraries(),
                m.nvitrines() * m.nshelves() * m.nbooks() * m.npages() *
                    m.nrows() * m.ncols()),
          Range(0,
                m.nvitrines(),
                m.nshelves() * m.nbooks() * m.npages() * m.nrows() * m.ncols()),
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
  std::memcpy(mdata,
              m.mdata,
              nlibraries() * nvitrines() * nshelves() * nbooks() * npages() *
                  nrows() * ncols() * sizeof(Numeric));
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
Tensor7& Tensor7::operator=(const Tensor7& x) {
  if (this != &x) {
    resize(x.nlibraries(),
           x.nvitrines(),
           x.nshelves(),
           x.nbooks(),
           x.npages(),
           x.nrows(),
           x.ncols());
    std::memcpy(mdata,
                x.mdata,
                nlibraries() * nvitrines() * nshelves() * nbooks() * npages() *
                    nrows() * ncols() * sizeof(Numeric));
  }
  return *this;
}

//! Copy assignment operator from another tensor.
Tensor7& Tensor7::operator=(Tensor7&& x) noexcept {
  if (this != &x) {
    delete[] mdata;
    mdata = x.mdata;
    mlr = x.mlr;
    mvr = x.mvr;
    msr = x.msr;
    mbr = x.mbr;
    mpr = x.mpr;
    mrr = x.mrr;
    mcr = x.mcr;
    x.mlr = Range(0, 0);
    x.mvr = Range(0, 0);
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
Tensor7& Tensor7::operator=(Numeric x) {
  std::fill_n(mdata,
              nlibraries() * nvitrines() * nshelves() * nbooks() * npages() *
                  nrows() * ncols(),
              x);
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values. */
void Tensor7::resize(
    Index l, Index v, Index s, Index b, Index p, Index r, Index c) {
  assert(0 <= l);
  assert(0 <= v);
  assert(0 <= s);
  assert(0 <= b);
  assert(0 <= p);
  assert(0 <= r);
  assert(0 <= c);

  if (mlr.mextent != l || mvr.mextent != v || msr.mextent != s ||
      mbr.mextent != b || mpr.mextent != p || mrr.mextent != r ||
      mcr.mextent != c) {
    delete[] mdata;
    mdata = new Numeric[l * v * s * b * p * r * c];

    mlr.mstart = 0;
    mlr.mextent = l;
    mlr.mstride = v * s * b * p * r * c;

    mvr.mstart = 0;
    mvr.mextent = v;
    mvr.mstride = s * b * p * r * c;

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
void swap(Tensor7& t1, Tensor7& t2) {
  std::swap(t1.mlr, t2.mlr);
  std::swap(t1.mvr, t2.mvr);
  std::swap(t1.msr, t2.msr);
  std::swap(t1.mbr, t2.mbr);
  std::swap(t1.mpr, t2.mpr);
  std::swap(t1.mrr, t2.mrr);
  std::swap(t1.mcr, t2.mcr);
  std::swap(t1.mdata, t2.mdata);
}

/** Destructor for Tensor7. This is important, since Tensor7 uses new to
    allocate storage. */
Tensor7::~Tensor7() {
  //   cout << "Destroying a Tensor7:\n"
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
void transform(Tensor7View y, double (&my_func)(double), ConstTensor7View x) {
  // Check dimensions:
  assert(y.nlibraries() == x.nlibraries());
  assert(y.nvitrines() == x.nvitrines());
  assert(y.nshelves() == x.nshelves());
  assert(y.nbooks() == x.nbooks());
  assert(y.npages() == x.npages());
  assert(y.nrows() == x.nrows());
  assert(y.ncols() == x.ncols());

  const ConstIterator7D xe = x.end();
  ConstIterator7D xi = x.begin();
  Iterator7D yi = y.begin();
  for (; xi != xe; ++xi, ++yi) {
    // Use the transform function of lower dimensional tensors
    // recursively:
    transform(*yi, my_func, *xi);
  }
}

/** Max function, tensor version. */
Numeric max(const ConstTensor7View& x) {
  const ConstIterator7D xe = x.end();
  ConstIterator7D xi = x.begin();

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
Numeric min(const ConstTensor7View& x) {
  const ConstIterator7D xe = x.end();
  ConstIterator7D xi = x.begin();

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
    \param l  Library index
    \param v  Vitrine index
    \param s  Shelf index
    \param b  Book index
    \param p  Page index
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Numeric debug_tensor7view_get_elem(Tensor7View& tv,
                                   Index l,
                                   Index v,
                                   Index s,
                                   Index b,
                                   Index p,
                                   Index r,
                                   Index c) {
  return tv(l, v, s, b, p, r, c);
}

#endif
////////////////////////////////
