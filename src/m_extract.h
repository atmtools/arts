/* Copyright (C) 2012 Stefan Buehler <sbuehler@ltu.se>
                      Oliver Lemke <olemke@core-dump.info>

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
  \file   m_extract.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-10-08
  
  \brief  Implementation of Extract.
  
  This file contains the implementation of the supergeneric method
  Extract.
*/

#ifndef m_extract_h
#define m_extract_h

#include "array.h"
#include "exceptions.h"
#include "gridded_fields.h"
#include "matpackV.h"
#include "quantum.h"

/* Workspace method: Doxygen documentation will be auto-generated */
/* This template function covers all implementations of extracting
   an element X from an ArrayOfX
 */
template <typename T>
void Extract(
    // WS Generic Output:
    T& e,
    // WS Input:
    // WS Generic Input:
    const Array<T>& arr,
    const Index& index,
    const Verbosity&) {
  if (index >= arr.nelem()) {
    ostringstream os;
    os << "The index " << index << " is outside the range of the array.";
    throw runtime_error(os.str());
  }

  e = arr[index];
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void ArrayOfIndexExtractFromArrayOfArrayOfIndex(
    // WS Generic Output:
    ArrayOfIndex& aoi,
    // WS Input:
    // WS Generic Input:
    const ArrayOfArrayOfIndex& aoaoi,
    const Index& index,
    const Verbosity&) {
  if (index >= aoaoi.nelem()) {
    ostringstream os;
    os << "The index " << index << " is outside the range of the Array.";
    throw runtime_error(os.str());
  }

  aoi.resize(aoaoi[index].nelem());
  aoi = aoaoi[index];
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Extract(
    // WS Generic Output:
    Numeric& n,
    // WS Input:
    // WS Generic Input:
    const Vector& v,
    const Index& index,
    const Verbosity&) {
  if (index >= v.nelem()) {
    ostringstream os;
    os << "The index " << index << " is outside the range of the Vector.";
    throw runtime_error(os.str());
  }

  n = v[index];
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Extract(
    // WS Generic Output:
    Matrix& m,
    // WS Input:
    // WS Generic Input:
    const Tensor3& t3,
    const Index& index,
    const Verbosity&) {
  if (index >= t3.npages()) {
    ostringstream os;
    os << "The index " << index << " is outside the page range of the Tensor3.";
    throw runtime_error(os.str());
  }

  m = t3(index, joker, joker);
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Extract(
    // WS Generic Output:
    Tensor3& t3,
    // WS Input:
    // WS Generic Input:
    const Tensor4& t4,
    const Index& index,
    const Verbosity&) {
  if (index >= t4.nbooks()) {
    ostringstream os;
    os << "The index " << index << " is outside the book range of the Tensor4.";
    throw runtime_error(os.str());
  }

  t3.resize(t4.npages(), t4.nrows(), t4.ncols());
  t3 = t4(index, joker, joker, joker);
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Extract(
    // WS Generic Output:
    Tensor4& t4,
    // WS Input:
    // WS Generic Input:
    const Tensor5& t5,
    const Index& index,
    const Verbosity&) {
  if (index >= t5.nshelves()) {
    ostringstream os;
    os << "The index " << index << "is outside the shelf range of the Tensor5.";
    throw runtime_error(os.str());
  }

  t4.resize(t5.nbooks(), t5.npages(), t5.nrows(), t5.ncols());
  t4 = t5(index, joker, joker, joker, joker);
}

/* Workspace method: Doxygen documentation will be auto-generated 

   Implementation largely copied from Patrick's MatrixExtractFromTensor3 method. 

   2007-10-26 Oliver Lemke */
inline void Extract(
    // WS Generic Output:
    ArrayOfGriddedField3& agf,
    // WS Input:
    // WS Generic Input:
    const ArrayOfArrayOfGriddedField3& aagf,
    const Index& index,
    const Verbosity&) {
  if (index >= aagf.nelem()) {
    ostringstream os;
    os << "The index " << index
       << " is outside the range of the ArrayOfArrayOfGriddedField3.";
    throw runtime_error(os.str());
  }

  agf.resize(aagf[index].nelem());
  agf = aagf[index];
}

/* Workspace method: Doxygen documentation will be auto-generated 

   Implementation largely copied from MatrixExtractFromArrayOfMatrix.

   2007-11-26 Stefan Buehler */
inline void Extract(
    // WS Generic Output:
    GriddedField4& m,
    // WS Input:
    // WS Generic Input:
    const ArrayOfGriddedField4& agf4,
    const Index& index,
    const Verbosity&) {
  if (index >= agf4.nelem()) {
    ostringstream os;
    os << "The index " << index
       << " is outside the range of The ArrayOfGriddedField4.";
    throw runtime_error(os.str());
  }

  // I simply use the copy operator here, since I'm too lazy to go
  // through all members of the structure to resize them. That is not
  // necessary, since sizes are adjusted automatically.
  m = agf4[index];
}

inline void Extract(
    // WS Generic Output:
    QuantumIdentifier& qi,
    // WS Input:
    // WS Generic Input:
    const ArrayOfQuantumIdentifier& aoqi,
    const Index& index,
    const Verbosity&) {
  if (index > aoqi.nelem() or index < 0) throw std::runtime_error("Bad index");
  qi = aoqi[index];
}

#endif /* m_extract_h */
