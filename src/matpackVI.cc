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
   \file   matpackVI.cc

   \author Oliver Lemke
   \date   2002-11-21
*/

#include "matpackVI.h"

#include "exceptions.h"

// Functions for ConstTensor6View:
// ------------------------------

// Const index operators:

// Result 6D (1 combination)
// ------
ConstTensor6View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  return ConstTensor6View(
      mdata, mvr, msr, mbr, mpr, mrr, mcr, v, s, b, p, r, c);
}

// Result 5D (6 combinations)
// -----|
ConstTensor5View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(c);
  return ConstTensor5View(
      mdata + OFFSET(c), mvr, msr, mbr, mpr, mrr, v, s, b, p, r);
}

// ----|-
ConstTensor5View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(r);
  return ConstTensor5View(
      mdata + OFFSET(r), mvr, msr, mbr, mpr, mcr, v, s, b, p, c);
}

// ---|--
ConstTensor5View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(p);
  return ConstTensor5View(
      mdata + OFFSET(p), mvr, msr, mbr, mrr, mcr, v, s, b, r, c);
}

// --|---
ConstTensor5View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(b);
  return ConstTensor5View(
      mdata + OFFSET(b), mvr, msr, mpr, mrr, mcr, v, s, p, r, c);
}

// -|----
ConstTensor5View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  return ConstTensor5View(
      mdata + OFFSET(s), mvr, mbr, mpr, mrr, mcr, v, b, p, r, c);
}

// |-----
ConstTensor5View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  return ConstTensor5View(
      mdata + OFFSET(v), msr, mbr, mpr, mrr, mcr, s, b, p, r, c);
}

// Result 4D (5+4+3+2+1 = 15 combinations)
// ----||
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(r);
  CHECK(c);
  return ConstTensor4View(
      mdata + OFFSET(r) + OFFSET(c), mvr, msr, mbr, mpr, v, s, b, p);
}

// ---|-|
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(p);
  CHECK(c);
  return ConstTensor4View(
      mdata + OFFSET(p) + OFFSET(c), mvr, msr, mbr, mrr, v, s, b, r);
}

// --|--|
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(b);
  CHECK(c);
  return ConstTensor4View(
      mdata + OFFSET(b) + OFFSET(c), mvr, msr, mpr, mrr, v, s, p, r);
}

// -|---|
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(c);
  return ConstTensor4View(
      mdata + OFFSET(s) + OFFSET(c), mvr, mbr, mpr, mrr, v, b, p, r);
}

// |----|
ConstTensor4View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(c);
  return ConstTensor4View(
      mdata + OFFSET(v) + OFFSET(c), msr, mbr, mpr, mrr, s, b, p, r);
}

// ---||-
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(p);
  CHECK(r);
  return ConstTensor4View(
      mdata + OFFSET(p) + OFFSET(r), mvr, msr, mbr, mcr, v, s, b, c);
}

// --|-|-
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(b);
  CHECK(r);
  return ConstTensor4View(
      mdata + OFFSET(b) + OFFSET(r), mvr, msr, mpr, mcr, v, s, p, c);
}

// -|--|-
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(r);
  return ConstTensor4View(
      mdata + OFFSET(s) + OFFSET(r), mvr, mbr, mpr, mcr, v, b, p, c);
}

// |---|-
ConstTensor4View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(r);
  return ConstTensor4View(
      mdata + OFFSET(v) + OFFSET(r), msr, mbr, mpr, mcr, s, b, p, c);
}

// --||--
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(b);
  CHECK(p);
  return ConstTensor4View(
      mdata + OFFSET(b) + OFFSET(p), mvr, msr, mrr, mcr, v, s, r, c);
}

// -|-|--
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(p);
  return ConstTensor4View(
      mdata + OFFSET(s) + OFFSET(p), mvr, mbr, mrr, mcr, v, b, r, c);
}

// |--|--
ConstTensor4View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(p);
  return ConstTensor4View(
      mdata + OFFSET(v) + OFFSET(p), msr, mbr, mrr, mcr, s, b, r, c);
}

// -||---
ConstTensor4View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  return ConstTensor4View(
      mdata + OFFSET(s) + OFFSET(b), mvr, mpr, mrr, mcr, v, p, r, c);
}

// |-|---
ConstTensor4View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  return ConstTensor4View(
      mdata + OFFSET(v) + OFFSET(b), msr, mpr, mrr, mcr, s, p, r, c);
}

// ||----
ConstTensor4View ConstTensor6View::operator()(Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  return ConstTensor4View(
      mdata + OFFSET(v) + OFFSET(s), mbr, mpr, mrr, mcr, b, p, r, c);
}

// Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
// ---|||
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              Index c) const {
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(p) + OFFSET(r) + OFFSET(c), mvr, msr, mbr, v, s, b);
}

// --|-||
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(b) + OFFSET(r) + OFFSET(c), mvr, msr, mpr, v, s, p);
}

// -|--||
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(s) + OFFSET(r) + OFFSET(c), mvr, mbr, mpr, v, b, p);
}

// |---||
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              Index c) const {
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(r) + OFFSET(c), msr, mbr, mpr, s, b, p);
}

// --||-|
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(b) + OFFSET(p) + OFFSET(c), mvr, msr, mrr, v, s, r);
}

// -|-|-|
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(s) + OFFSET(p) + OFFSET(c), mvr, mbr, mrr, v, b, r);
}

// |--|-|
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(p) + OFFSET(c), msr, mbr, mrr, s, b, r);
}

// -||--|
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(c), mvr, mpr, mrr, v, p, r);
}

// |-|--|
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(c), msr, mpr, mrr, s, p, r);
}

// ||---|
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              const Range& r,
                                              Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(c), mbr, mpr, mrr, b, p, r);
}

// --|||-
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(
      mdata + OFFSET(b) + OFFSET(p) + OFFSET(r), mvr, msr, mcr, v, s, c);
}

// -|-||-
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(
      mdata + OFFSET(s) + OFFSET(p) + OFFSET(r), mvr, mbr, mcr, v, b, c);
}

// |--||-
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              const Range& b,
                                              Index p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(p) + OFFSET(r), msr, mbr, mcr, s, b, c);
}

// -||-|-
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(r), mvr, mpr, mcr, v, p, c);
}

// |-|-|-
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              Index b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(r), msr, mpr, mcr, s, p, c);
}

// ||--|-
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              Index s,
                                              const Range& b,
                                              const Range& p,
                                              Index r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(r), mbr, mpr, mcr, b, p, c);
}

// -|||--
ConstTensor3View ConstTensor6View::operator()(const Range& v,
                                              Index s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p), mvr, mrr, mcr, v, r, c);
}

// |-||--
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              const Range& s,
                                              Index b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p), msr, mrr, mcr, s, r, c);
}

// ||-|--
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              Index s,
                                              const Range& b,
                                              Index p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p), mbr, mrr, mcr, b, r, c);
}

// |||---
ConstTensor3View ConstTensor6View::operator()(Index v,
                                              Index s,
                                              Index b,
                                              const Range& p,
                                              const Range& r,
                                              const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return ConstTensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b), mpr, mrr, mcr, p, r, c);
}

// Result 2D (15 combinations)
// IIII--
ConstMatrixView ConstTensor6View::operator()(
    Index v, Index s, Index b, Index p, const Range& r, const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p), mrr, mcr, r, c);
}

// III-I-
ConstMatrixView ConstTensor6View::operator()(
    Index v, Index s, Index b, const Range& p, Index r, const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r), mpr, mcr, p, c);
}

// II-II-
ConstMatrixView ConstTensor6View::operator()(
    Index v, Index s, const Range& b, Index p, Index r, const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r), mbr, mcr, b, c);
}

// I-III-
ConstMatrixView ConstTensor6View::operator()(
    Index v, const Range& s, Index b, Index p, Index r, const Range& c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r), msr, mcr, s, c);
}

// -IIII-
ConstMatrixView ConstTensor6View::operator()(
    const Range& v, Index s, Index b, Index p, Index r, const Range& c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstMatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r), mvr, mcr, v, c);
}

// III--I
ConstMatrixView ConstTensor6View::operator()(
    Index v, Index s, Index b, const Range& p, const Range& r, Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c), mpr, mrr, p, r);
}

// II-I-I
ConstMatrixView ConstTensor6View::operator()(
    Index v, Index s, const Range& b, Index p, const Range& r, Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c), mbr, mrr, b, r);
}

// I-II-I
ConstMatrixView ConstTensor6View::operator()(
    Index v, const Range& s, Index b, Index p, const Range& r, Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c), msr, mrr, s, r);
}

// -III-I
ConstMatrixView ConstTensor6View::operator()(
    const Range& v, Index s, Index b, Index p, const Range& r, Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c), mvr, mrr, v, r);
}

// II--II
ConstMatrixView ConstTensor6View::operator()(
    Index v, Index s, const Range& b, const Range& p, Index r, Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c), mbr, mpr, b, p);
}

// I-I-II
ConstMatrixView ConstTensor6View::operator()(
    Index v, const Range& s, Index b, const Range& p, Index r, Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c), msr, mpr, s, p);
}

// -II-II
ConstMatrixView ConstTensor6View::operator()(
    const Range& v, Index s, Index b, const Range& p, Index r, Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c), mvr, mpr, v, p);
}

// I--III
ConstMatrixView ConstTensor6View::operator()(
    Index v, const Range& s, const Range& b, Index p, Index r, Index c) const {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c), msr, mbr, s, b);
}

// -I-III
ConstMatrixView ConstTensor6View::operator()(
    const Range& v, Index s, const Range& b, Index p, Index r, Index c) const {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c), mvr, mbr, v, b);
}

// --IIII
ConstMatrixView ConstTensor6View::operator()(
    const Range& v, const Range& s, Index b, Index p, Index r, Index c) const {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstMatrixView(
      mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c), mvr, msr, v, s);
}

// Result 1D (6 combinations)
// IIIII-
ConstVectorView ConstTensor6View::operator()(
    Index v, Index s, Index b, Index p, Index r, const Range& c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return ConstVectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      mcr,
      c);
}

// IIII-I
ConstVectorView ConstTensor6View::operator()(
    Index v, Index s, Index b, Index p, const Range& r, Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return ConstVectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      mrr,
      r);
}

// III-II
ConstVectorView ConstTensor6View::operator()(
    Index v, Index s, Index b, const Range& p, Index r, Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      mpr,
      p);
}

// II-III
ConstVectorView ConstTensor6View::operator()(
    Index v, Index s, const Range& b, Index p, Index r, Index c) const {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mbr,
      b);
}

// I-IIII
ConstVectorView ConstTensor6View::operator()(
    Index v, const Range& s, Index b, Index p, Index r, Index c) const {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      msr,
      s);
}

// -IIIII
ConstVectorView ConstTensor6View::operator()(
    const Range& v, Index s, Index b, Index p, Index r, Index c) const {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return ConstVectorView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mvr,
      v);
}

/** Return const iterator to first sub-tensor. */
ConstIterator6D ConstTensor6View::begin() const {
  return ConstIterator6D(
      ConstTensor5View(mdata + mvr.mstart, msr, mbr, mpr, mrr, mcr),
      mvr.mstride);
}

/** Return const iterator behind last sub-tensor. */
ConstIterator6D ConstTensor6View::end() const {
  return ConstIterator6D(
      ConstTensor5View(mdata + mvr.mstart + (mvr.mextent) * mvr.mstride,
                       msr,
                       mbr,
                       mpr,
                       mrr,
                       mcr),
      mvr.mstride);
}

/** Special constructor to make a Tensor6 view of a Tensor5. */
ConstTensor6View::ConstTensor6View(const ConstTensor5View& a)
    : mvr(0,
          1,
          a.msr.mextent * a.mbr.mextent * a.mpr.mextent * a.mrr.mextent *
              a.mcr.mextent),
      msr(a.msr),
      mbr(a.mbr),
      mpr(a.mpr),
      mrr(a.mrr),
      mcr(a.mcr),
      mdata(a.mdata) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor6 to initialize
    its own Tensor6View part. The row range rr must have a stride to
    account for the length of one row. The page range pr must have a
    stride to account for the length of one page. */
ConstTensor6View::ConstTensor6View(Numeric* data,
                                   const Range& v,
                                   const Range& s,
                                   const Range& b,
                                   const Range& p,
                                   const Range& r,
                                   const Range& c)
    : mvr(v), msr(s), mbr(b), mpr(p), mrr(r), mcr(c), mdata(data) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
ConstTensor6View::ConstTensor6View(Numeric* data,
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
                                   const Range& nc)
    : mvr(pv, nv),
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
    Tensor5 to print each page in turn. */
std::ostream& operator<<(std::ostream& os, const ConstTensor6View& v) {
  // Page iterators:
  ConstIterator6D ip = v.begin();
  const ConstIterator6D end_page = v.end();

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

// Functions for Tensor6View:
// -------------------------

// Non-const index operators:

// Result 6D (1 combination)
// ------
Tensor6View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  return Tensor6View(mdata, mvr, msr, mbr, mpr, mrr, mcr, v, s, b, p, r, c);
}

// Result 5D (6 combinations)
// -----|
Tensor5View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(c);
  return Tensor5View(mdata + OFFSET(c), mvr, msr, mbr, mpr, mrr, v, s, b, p, r);
}

// ----|-
Tensor5View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(r);
  return Tensor5View(mdata + OFFSET(r), mvr, msr, mbr, mpr, mcr, v, s, b, p, c);
}

// ---|--
Tensor5View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(p);
  return Tensor5View(mdata + OFFSET(p), mvr, msr, mbr, mrr, mcr, v, s, b, r, c);
}

// --|---
Tensor5View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(b);
  return Tensor5View(mdata + OFFSET(b), mvr, msr, mpr, mrr, mcr, v, s, p, r, c);
}

// -|----
Tensor5View Tensor6View::operator()(const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  return Tensor5View(mdata + OFFSET(s), mvr, mbr, mpr, mrr, mcr, v, b, p, r, c);
}

// |-----
Tensor5View Tensor6View::operator()(Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  return Tensor5View(mdata + OFFSET(v), msr, mbr, mpr, mrr, mcr, s, b, p, r, c);
}

// Result 4D (5+4+3+2+1 = 15 combinations)
// ----||
Tensor4View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    Index c) {
  CHECK(r);
  CHECK(c);
  return Tensor4View(
      mdata + OFFSET(r) + OFFSET(c), mvr, msr, mbr, mpr, v, s, b, p);
}

// ---|-|
Tensor4View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    Index c) {
  CHECK(p);
  CHECK(c);
  return Tensor4View(
      mdata + OFFSET(p) + OFFSET(c), mvr, msr, mbr, mrr, v, s, b, r);
}

// --|--|
Tensor4View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(b);
  CHECK(c);
  return Tensor4View(
      mdata + OFFSET(b) + OFFSET(c), mvr, msr, mpr, mrr, v, s, p, r);
}

// -|---|
Tensor4View Tensor6View::operator()(const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(s);
  CHECK(c);
  return Tensor4View(
      mdata + OFFSET(s) + OFFSET(c), mvr, mbr, mpr, mrr, v, b, p, r);
}

// |----|
Tensor4View Tensor6View::operator()(Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    Index c) {
  CHECK(v);
  CHECK(c);
  return Tensor4View(
      mdata + OFFSET(v) + OFFSET(c), msr, mbr, mpr, mrr, s, b, p, r);
}

// ---||-
Tensor4View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    Index r,
                                    const Range& c) {
  CHECK(p);
  CHECK(r);
  return Tensor4View(
      mdata + OFFSET(p) + OFFSET(r), mvr, msr, mbr, mcr, v, s, b, c);
}

// --|-|-
Tensor4View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(b);
  CHECK(r);
  return Tensor4View(
      mdata + OFFSET(b) + OFFSET(r), mvr, msr, mpr, mcr, v, s, p, c);
}

// -|--|-
Tensor4View Tensor6View::operator()(const Range& v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(s);
  CHECK(r);
  return Tensor4View(
      mdata + OFFSET(s) + OFFSET(r), mvr, mbr, mpr, mcr, v, b, p, c);
}

// |---|-
Tensor4View Tensor6View::operator()(Index v,
                                    const Range& s,
                                    const Range& b,
                                    const Range& p,
                                    Index r,
                                    const Range& c) {
  CHECK(v);
  CHECK(r);
  return Tensor4View(
      mdata + OFFSET(v) + OFFSET(r), msr, mbr, mpr, mcr, s, b, p, c);
}

// --||--
Tensor4View Tensor6View::operator()(const Range& v,
                                    const Range& s,
                                    Index b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(b);
  CHECK(p);
  return Tensor4View(
      mdata + OFFSET(b) + OFFSET(p), mvr, msr, mrr, mcr, v, s, r, c);
}

// -|-|--
Tensor4View Tensor6View::operator()(const Range& v,
                                    Index s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  CHECK(p);
  return Tensor4View(
      mdata + OFFSET(s) + OFFSET(p), mvr, mbr, mrr, mcr, v, b, r, c);
}

// |--|--
Tensor4View Tensor6View::operator()(Index v,
                                    const Range& s,
                                    const Range& b,
                                    Index p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(p);
  return Tensor4View(
      mdata + OFFSET(v) + OFFSET(p), msr, mbr, mrr, mcr, s, b, r, c);
}

// -||---
Tensor4View Tensor6View::operator()(const Range& v,
                                    Index s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(s);
  CHECK(b);
  return Tensor4View(
      mdata + OFFSET(s) + OFFSET(b), mvr, mpr, mrr, mcr, v, p, r, c);
}

// |-|---
Tensor4View Tensor6View::operator()(Index v,
                                    const Range& s,
                                    Index b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(b);
  return Tensor4View(
      mdata + OFFSET(v) + OFFSET(b), msr, mpr, mrr, mcr, s, p, r, c);
}

// ||----
Tensor4View Tensor6View::operator()(Index v,
                                    Index s,
                                    const Range& b,
                                    const Range& p,
                                    const Range& r,
                                    const Range& c) {
  CHECK(v);
  CHECK(s);
  return Tensor4View(
      mdata + OFFSET(v) + OFFSET(s), mbr, mpr, mrr, mcr, b, p, r, c);
}

// Result 3D (4+3+2+1+ 3+2+1+ 2+1 +1 = 20 combinations)
// ---|||
Tensor3View Tensor6View::operator()(
    const Range& v, const Range& s, const Range& b, Index p, Index r, Index c) {
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(p) + OFFSET(r) + OFFSET(c), mvr, msr, mbr, v, s, b);
}

// --|-||
Tensor3View Tensor6View::operator()(
    const Range& v, const Range& s, Index b, const Range& p, Index r, Index c) {
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(b) + OFFSET(r) + OFFSET(c), mvr, msr, mpr, v, s, p);
}

// -|--||
Tensor3View Tensor6View::operator()(
    const Range& v, Index s, const Range& b, const Range& p, Index r, Index c) {
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(s) + OFFSET(r) + OFFSET(c), mvr, mbr, mpr, v, b, p);
}

// |---||
Tensor3View Tensor6View::operator()(
    Index v, const Range& s, const Range& b, const Range& p, Index r, Index c) {
  CHECK(v);
  CHECK(r);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(r) + OFFSET(c), msr, mbr, mpr, s, b, p);
}

// --||-|
Tensor3View Tensor6View::operator()(
    const Range& v, const Range& s, Index b, Index p, const Range& r, Index c) {
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(b) + OFFSET(p) + OFFSET(c), mvr, msr, mrr, v, s, r);
}

// -|-|-|
Tensor3View Tensor6View::operator()(
    const Range& v, Index s, const Range& b, Index p, const Range& r, Index c) {
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(s) + OFFSET(p) + OFFSET(c), mvr, mbr, mrr, v, b, r);
}

// |--|-|
Tensor3View Tensor6View::operator()(
    Index v, const Range& s, const Range& b, Index p, const Range& r, Index c) {
  CHECK(v);
  CHECK(p);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(p) + OFFSET(c), msr, mbr, mrr, s, b, r);
}

// -||--|
Tensor3View Tensor6View::operator()(
    const Range& v, Index s, Index b, const Range& p, const Range& r, Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(c), mvr, mpr, mrr, v, p, r);
}

// |-|--|
Tensor3View Tensor6View::operator()(
    Index v, const Range& s, Index b, const Range& p, const Range& r, Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(c), msr, mpr, mrr, s, p, r);
}

// ||---|
Tensor3View Tensor6View::operator()(
    Index v, Index s, const Range& b, const Range& p, const Range& r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(c);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(c), mbr, mpr, mrr, b, p, r);
}

// --|||-
Tensor3View Tensor6View::operator()(
    const Range& v, const Range& s, Index b, Index p, Index r, const Range& c) {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return Tensor3View(
      mdata + OFFSET(b) + OFFSET(p) + OFFSET(r), mvr, msr, mcr, v, s, c);
}

// -|-||-
Tensor3View Tensor6View::operator()(
    const Range& v, Index s, const Range& b, Index p, Index r, const Range& c) {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return Tensor3View(
      mdata + OFFSET(s) + OFFSET(p) + OFFSET(r), mvr, mbr, mcr, v, b, c);
}

// |--||-
Tensor3View Tensor6View::operator()(
    Index v, const Range& s, const Range& b, Index p, Index r, const Range& c) {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(p) + OFFSET(r), msr, mbr, mcr, s, b, c);
}

// -||-|-
Tensor3View Tensor6View::operator()(
    const Range& v, Index s, Index b, const Range& p, Index r, const Range& c) {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return Tensor3View(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(r), mvr, mpr, mcr, v, p, c);
}

// |-|-|-
Tensor3View Tensor6View::operator()(
    Index v, const Range& s, Index b, const Range& p, Index r, const Range& c) {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(r), msr, mpr, mcr, s, p, c);
}

// ||--|-
Tensor3View Tensor6View::operator()(
    Index v, Index s, const Range& b, const Range& p, Index r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(r), mbr, mpr, mcr, b, p, c);
}

// -|||--
Tensor3View Tensor6View::operator()(
    const Range& v, Index s, Index b, Index p, const Range& r, const Range& c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return Tensor3View(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p), mvr, mrr, mcr, v, r, c);
}

// |-||--
Tensor3View Tensor6View::operator()(
    Index v, const Range& s, Index b, Index p, const Range& r, const Range& c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p), msr, mrr, mcr, s, r, c);
}

// ||-|--
Tensor3View Tensor6View::operator()(
    Index v, Index s, const Range& b, Index p, const Range& r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p), mbr, mrr, mcr, b, r, c);
}

// |||---
Tensor3View Tensor6View::operator()(
    Index v, Index s, Index b, const Range& p, const Range& r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  return Tensor3View(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b), mpr, mrr, mcr, p, r, c);
}

// Result 2D (15 combinations)
// IIII--
MatrixView Tensor6View::operator()(
    Index v, Index s, Index b, Index p, const Range& r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p), mrr, mcr, r, c);
}

// III-I-
MatrixView Tensor6View::operator()(
    Index v, Index s, Index b, const Range& p, Index r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r), mpr, mcr, p, c);
}

// II-II-
MatrixView Tensor6View::operator()(
    Index v, Index s, const Range& b, Index p, Index r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r), mbr, mcr, b, c);
}

// I-III-
MatrixView Tensor6View::operator()(
    Index v, const Range& s, Index b, Index p, Index r, const Range& c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r), msr, mcr, s, c);
}

// -IIII-
MatrixView Tensor6View::operator()(
    const Range& v, Index s, Index b, Index p, Index r, const Range& c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return MatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r), mvr, mcr, v, c);
}

// III--I
MatrixView Tensor6View::operator()(
    Index v, Index s, Index b, const Range& p, const Range& r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(c), mpr, mrr, p, r);
}

// II-I-I
MatrixView Tensor6View::operator()(
    Index v, Index s, const Range& b, Index p, const Range& r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(c), mbr, mrr, b, r);
}

// I-II-I
MatrixView Tensor6View::operator()(
    Index v, const Range& s, Index b, Index p, const Range& r, Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(c), msr, mrr, s, r);
}

// -III-I
MatrixView Tensor6View::operator()(
    const Range& v, Index s, Index b, Index p, const Range& r, Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c), mvr, mrr, v, r);
}

// II--II
MatrixView Tensor6View::operator()(
    Index v, Index s, const Range& b, const Range& p, Index r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(r) + OFFSET(c), mbr, mpr, b, p);
}

// I-I-II
MatrixView Tensor6View::operator()(
    Index v, const Range& s, Index b, const Range& p, Index r, Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(r) + OFFSET(c), msr, mpr, s, p);
}

// -II-II
MatrixView Tensor6View::operator()(
    const Range& v, Index s, Index b, const Range& p, Index r, Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c), mvr, mpr, v, p);
}

// I--III
MatrixView Tensor6View::operator()(
    Index v, const Range& s, const Range& b, Index p, Index r, Index c) {
  CHECK(v);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(v) + OFFSET(p) + OFFSET(r) + OFFSET(c), msr, mbr, s, b);
}

// -I-III
MatrixView Tensor6View::operator()(
    const Range& v, Index s, const Range& b, Index p, Index r, Index c) {
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c), mvr, mbr, v, b);
}

// --IIII
MatrixView Tensor6View::operator()(
    const Range& v, const Range& s, Index b, Index p, Index r, Index c) {
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return MatrixView(
      mdata + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c), mvr, msr, v, s);
}

// Result 1D (6 combinations)
// IIIII-
VectorView Tensor6View::operator()(
    Index v, Index s, Index b, Index p, Index r, const Range& c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  return VectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r),
      mcr,
      c);
}

// IIII-I
VectorView Tensor6View::operator()(
    Index v, Index s, Index b, Index p, const Range& r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(c);
  return VectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(c),
      mrr,
      r);
}

// III-II
VectorView Tensor6View::operator()(
    Index v, Index s, Index b, const Range& p, Index r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(b);
  CHECK(r);
  CHECK(c);
  return VectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(b) + OFFSET(r) + OFFSET(c),
      mpr,
      p);
}

// II-III
VectorView Tensor6View::operator()(
    Index v, Index s, const Range& b, Index p, Index r, Index c) {
  CHECK(v);
  CHECK(s);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(
      mdata + OFFSET(v) + OFFSET(s) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mbr,
      b);
}

// I-IIII
VectorView Tensor6View::operator()(
    Index v, const Range& s, Index b, Index p, Index r, Index c) {
  CHECK(v);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(
      mdata + OFFSET(v) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      msr,
      s);
}

// -IIIII
VectorView Tensor6View::operator()(
    const Range& v, Index s, Index b, Index p, Index r, Index c) {
  CHECK(s);
  CHECK(b);
  CHECK(p);
  CHECK(r);
  CHECK(c);
  return VectorView(
      mdata + OFFSET(s) + OFFSET(b) + OFFSET(p) + OFFSET(r) + OFFSET(c),
      mvr,
      v);
}

/** Conversion to plain C-array.
    
This function returns a pointer to the raw data. It fails if the
Tensor5View is not pointing to the beginning of a Tensor5 or the stride
is not 1 because the caller expects to get a C array with continuous data.
*/
Numeric* Tensor6View::get_c_array() ARTS_NOEXCEPT {
  ARTS_ASSERT(mvr.mstart == 0 and
                  (mvr.mstride == mcr.mextent * mrr.mextent * mpr.mextent *
                                      mbr.mextent * msr.mextent or
                   size() == 0),
              "Vitrine ",
              mvr)
  ARTS_ASSERT(msr.mstart == 0 and
                  (msr.mstride ==
                       mcr.mextent * mrr.mextent * mpr.mextent * mbr.mextent or
                   size() == 0),
              "Shelve ",
              msr)
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
  Tensor5View is not pointing to the beginning of a Tensor5 or the stride
  is not 1 because the caller expects to get a C array with continuous data.
*/
const Numeric* Tensor6View::get_c_array() const ARTS_NOEXCEPT {
  ARTS_ASSERT(mvr.mstart == 0 and
                  (mvr.mstride == mcr.mextent * mrr.mextent * mpr.mextent *
                                      mbr.mextent * msr.mextent or
                   size() == 0),
              "Vitrine ",
              mvr)
  ARTS_ASSERT(msr.mstart == 0 and
                  (msr.mstride ==
                       mcr.mextent * mrr.mextent * mpr.mextent * mbr.mextent or
                   size() == 0),
              "Shelve ",
              msr)
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

/** Return iterator to first sub-tensor. */
Iterator6D Tensor6View::begin() {
  return Iterator6D(Tensor5View(mdata + mvr.mstart, msr, mbr, mpr, mrr, mcr),
                    mvr.mstride);
}

/** Return iterator behind last sub-tensor. */
Iterator6D Tensor6View::end() {
  return Iterator6D(
      Tensor5View(mdata + mvr.mstart + (mvr.mextent) * mvr.mstride,
                  msr,
                  mbr,
                  mpr,
                  mrr,
                  mcr),
      mvr.mstride);
}

/** Assignment operator. This copies the data from another Tensor6View
    to this Tensor6View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor6View by
    setting its range. */
Tensor6View& Tensor6View::operator=(const ConstTensor6View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mvr.mextent == m.mvr.mextent);
  ARTS_ASSERT(msr.mextent == m.msr.mextent);
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from Tensor6View to Tensor6View. This is a tricky
    one. The problem is that since Tensor6View is derived from
    ConstTensor6View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
Tensor6View& Tensor6View::operator=(const Tensor6View& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mvr.mextent == m.mvr.mextent);
  ARTS_ASSERT(msr.mextent == m.msr.mextent);
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assignment from a Tensor6. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
Tensor6View& Tensor6View::operator=(const Tensor6& m) {
  // Check that sizes are compatible:
  ARTS_ASSERT(mvr.mextent == m.mvr.mextent);
  ARTS_ASSERT(msr.mextent == m.msr.mextent);
  ARTS_ASSERT(mbr.mextent == m.mbr.mextent);
  ARTS_ASSERT(mpr.mextent == m.mpr.mextent);
  ARTS_ASSERT(mrr.mextent == m.mrr.mextent);
  ARTS_ASSERT(mcr.mextent == m.mcr.mextent);

  copy(m.begin(), m.end(), begin());
  return *this;
}

/** Assigning a scalar to a Tensor6View will set all elements to this
    value. */
Tensor6View& Tensor6View::operator=(Numeric x) {
  copy(x, begin(), end());
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
Tensor6View& Tensor6View::operator*=(Numeric x) {
  const Iterator6D ep = end();
  for (Iterator6D p = begin(); p != ep; ++p) {
    *p *= x;
  }
  return *this;
}

/** Division by scalar. */
Tensor6View& Tensor6View::operator/=(Numeric x) {
  const Iterator6D ep = end();
  for (Iterator6D p = begin(); p != ep; ++p) {
    *p /= x;
  }
  return *this;
}

/** Addition of scalar. */
Tensor6View& Tensor6View::operator+=(Numeric x) {
  const Iterator6D ep = end();
  for (Iterator6D p = begin(); p != ep; ++p) {
    *p += x;
  }
  return *this;
}

/** Subtraction of scalar. */
Tensor6View& Tensor6View::operator-=(Numeric x) {
  const Iterator6D ep = end();
  for (Iterator6D p = begin(); p != ep; ++p) {
    *p -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor6. */
Tensor6View& Tensor6View::operator*=(const ConstTensor6View& x) {
  ARTS_ASSERT(nvitrines() == x.nvitrines());
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator6D xp = x.begin();
  Iterator6D p = begin();
  const Iterator6D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p *= *xp;
  }
  return *this;
}

/** Element-vise division by another Tensor6. */
Tensor6View& Tensor6View::operator/=(const ConstTensor6View& x) {
  ARTS_ASSERT(nvitrines() == x.nvitrines());
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator6D xp = x.begin();
  Iterator6D p = begin();
  const Iterator6D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p /= *xp;
  }
  return *this;
}

/** Element-vise addition of another Tensor6. */
Tensor6View& Tensor6View::operator+=(const ConstTensor6View& x) {
  ARTS_ASSERT(nvitrines() == x.nvitrines());
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator6D xp = x.begin();
  Iterator6D p = begin();
  const Iterator6D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p += *xp;
  }
  return *this;
}

/** Element-vise subtraction of another Tensor6. */
Tensor6View& Tensor6View::operator-=(const ConstTensor6View& x) {
  ARTS_ASSERT(nvitrines() == x.nvitrines());
  ARTS_ASSERT(nshelves() == x.nshelves());
  ARTS_ASSERT(nbooks() == x.nbooks());
  ARTS_ASSERT(npages() == x.npages());
  ARTS_ASSERT(nrows() == x.nrows());
  ARTS_ASSERT(ncols() == x.ncols());
  ConstIterator6D xp = x.begin();
  Iterator6D p = begin();
  const Iterator6D ep = end();
  for (; p != ep; ++p, ++xp) {
    *p -= *xp;
  }
  return *this;
}

/** Special constructor to make a Tensor6 view of a Tensor5. */
Tensor6View::Tensor6View(const Tensor5View& a)
    : ConstTensor6View(a.mdata,
                       Range(0,
                             1,
                             a.msr.mextent * a.mbr.mextent * a.mpr.mextent *
                                 a.mrr.mextent * a.mcr.mextent),
                       a.msr,
                       a.mbr,
                       a.mpr,
                       a.mrr,
                       a.mcr) {
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor6 to initialize its
    own Tensor6View part. The row range rr must have a
    stride to account for the length of one row. */
Tensor6View::Tensor6View(Numeric* data,
                         const Range& v,
                         const Range& s,
                         const Range& b,
                         const Range& p,
                         const Range& r,
                         const Range& c)
    : ConstTensor6View(data, v, s, b, p, r, c) {
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges. 

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param pv Previous range.
    \param ps Previous range.
    \param pb Previous range.
    \param pp Previous range.
    \param pr Previous range.
    \param pc Previous range.
    \param nv New Range.
    \param ns New Range.
    \param nb New Range.
    \param np New Range.
    \param nr New Range.
    \param nc New Range.
  */
Tensor6View::Tensor6View(Numeric* data,
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
                         const Range& nc)
    : ConstTensor6View(data, pv, ps, pb, pp, pr, pc, nv, ns, nb, np, nr, nc) {
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
void copy(ConstIterator6D origin,
          const ConstIterator6D& end,
          Iterator6D target) {
  for (; origin != end; ++origin, ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(origin->begin(), origin->end(), target->begin());
  }
}

/** Copy a scalar to all elements. */
void copy(Numeric x, Iterator6D target, const Iterator6D& end) {
  for (; target != end; ++target) {
    // We use the copy function for the next smaller rank of tensor
    // recursively:
    copy(x, target->begin(), target->end());
  }
}

// Functions for Tensor6:
// ---------------------

/** Constructor setting size. This constructor has to set the strides
    in the page and row ranges correctly! */
Tensor6::Tensor6(Index v, Index s, Index b, Index p, Index r, Index c)
    : Tensor6View(new Numeric[v * s * b * p * r * c],
                  Range(0, v, s * b * p * r * c),
                  Range(0, s, b * p * r * c),
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
Tensor6::Tensor6(
    Index v, Index s, Index b, Index p, Index r, Index c, Numeric fill)
    : Tensor6View(new Numeric[v * s * b * p * r * c],
                  Range(0, v, s * b * p * r * c),
                  Range(0, s, b * p * r * c),
                  Range(0, b, p * r * c),
                  Range(0, p, r * c),
                  Range(0, r, c),
                  Range(0, c)) {
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  std::fill_n(mdata, v * s * b * p * r * c, fill);
}

/** Copy constructor from Tensor6View. This automatically sets the size
    and copies the data. */
Tensor6::Tensor6(const ConstTensor6View& m)
    : Tensor6View(
          new Numeric[m.nvitrines() * m.nshelves() * m.nbooks() * m.npages() *
                      m.nrows() * m.ncols()],
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

/** Copy constructor from Tensor6. This automatically sets the size
    and copies the data. */
Tensor6::Tensor6(const Tensor6& m)
    : Tensor6View(
          new Numeric[m.nvitrines() * m.nshelves() * m.nbooks() * m.npages() *
                      m.nrows() * m.ncols()],
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
              nvitrines() * nshelves() * nbooks() * npages() * nrows() *
                  ncols() * sizeof(Numeric));
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
Tensor6& Tensor6::operator=(const Tensor6& x) {
  if (this != &x) {
    resize(x.nvitrines(),
           x.nshelves(),
           x.nbooks(),
           x.npages(),
           x.nrows(),
           x.ncols());
    std::memcpy(mdata,
                x.mdata,
                nvitrines() * nshelves() * nbooks() * npages() * nrows() *
                    ncols() * sizeof(Numeric));
  }
  return *this;
}

//! Move assignment operator from another tensor.
Tensor6& Tensor6::operator=(Tensor6&& x) noexcept {
  if (this != &x) {
    delete[] mdata;
    mdata = x.mdata;
    mvr = x.mvr;
    msr = x.msr;
    mbr = x.mbr;
    mpr = x.mpr;
    mrr = x.mrr;
    mcr = x.mcr;
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
Tensor6& Tensor6::operator=(Numeric x) {
  std::fill_n(
      mdata,
      nvitrines() * nshelves() * nbooks() * npages() * nrows() * ncols(),
      x);
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values. */
void Tensor6::resize(Index v, Index s, Index b, Index p, Index r, Index c) {
  ARTS_ASSERT(0 <= v);
  ARTS_ASSERT(0 <= s);
  ARTS_ASSERT(0 <= b);
  ARTS_ASSERT(0 <= p);
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  if (mvr.mextent != v || msr.mextent != s || mbr.mextent != b ||
      mpr.mextent != p || mrr.mextent != r || mcr.mextent != c) {
    delete[] mdata;
    mdata = new Numeric[v * s * b * p * r * c];

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
void swap(Tensor6& t1, Tensor6& t2) noexcept {
  using std::swap;
  swap(t1.mvr, t2.mvr);
  swap(t1.msr, t2.msr);
  swap(t1.mbr, t2.mbr);
  swap(t1.mpr, t2.mpr);
  swap(t1.mrr, t2.mrr);
  swap(t1.mcr, t2.mcr);
  swap(t1.mdata, t2.mdata);
}

/** Destructor for Tensor6. This is important, since Tensor6 uses new to
    allocate storage. */
Tensor6::~Tensor6() noexcept {
  //   cout << "Destroying a Tensor6:\n"
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
void transform(Tensor6View y, double (&my_func)(double), ConstTensor6View x) {
  // Check dimensions:
  ARTS_ASSERT(y.nvitrines() == x.nvitrines());
  ARTS_ASSERT(y.nshelves() == x.nshelves());
  ARTS_ASSERT(y.nbooks() == x.nbooks());
  ARTS_ASSERT(y.npages() == x.npages());
  ARTS_ASSERT(y.nrows() == x.nrows());
  ARTS_ASSERT(y.ncols() == x.ncols());

  const ConstIterator6D xe = x.end();
  ConstIterator6D xi = x.begin();
  Iterator6D yi = y.begin();
  for (; xi != xe; ++xi, ++yi) {
    // Use the transform function of lower dimensional tensors
    // recursively:
    transform(*yi, my_func, *xi);
  }
}

/** Max function, tensor version. */
Numeric max(const ConstTensor6View& x) {
  const ConstIterator6D xe = x.end();
  ConstIterator6D xi = x.begin();

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
Numeric min(const ConstTensor6View& x) {
  const ConstIterator6D xe = x.end();
  ConstIterator6D xi = x.begin();

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
    \param v  Vitrine index
    \param s  Shelf index
    \param b  Book index
    \param p  Page index
    \param r  Row index
    \param c  Column index

    \author Oliver Lemke
    \date   2004-05-10
*/
Numeric debug_tensor6view_get_elem(
    Tensor6View& tv, Index v, Index s, Index b, Index p, Index r, Index c) {
  return tv(v, s, b, p, r, c);
}

#endif
////////////////////////////////
