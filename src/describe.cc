/* Copyright (C) 2003-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   describe.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Feb 25 15:35:56 2003
  
  \brief  Describe type and dimensions of a tensor view.
  
  This file contains a set of helper functions called \c describe,
  which you can use to output the dimensions of a tensor. This is just
  for testing purposes.
*/

#include "describe.h"

using std::ostringstream;

//! Describe Tensor7.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstTensor7View x) {
  ostringstream os;
  os << "Tensor7 [" << x.nlibraries() << "," << x.nvitrines() << ","
     << x.nshelves() << "," << x.nbooks() << "," << x.npages() << ","
     << x.nrows() << "," << x.ncols() << "]";
  return os.str();
}

//! Describe Tensor6.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstTensor6View x) {
  ostringstream os;
  os << "Tensor6 [" << x.nvitrines() << "," << x.nshelves() << "," << x.nbooks()
     << "," << x.npages() << "," << x.nrows() << "," << x.ncols() << "]";
  return os.str();
}

//! Describe Tensor5.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstTensor5View x) {
  ostringstream os;
  os << "Tensor5 [" << x.nshelves() << "," << x.nbooks() << "," << x.npages()
     << "," << x.nrows() << "," << x.ncols() << "]";
  return os.str();
}

//! Describe Tensor4.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstTensor4View x) {
  ostringstream os;
  os << "Tensor4 [" << x.nbooks() << "," << x.npages() << "," << x.nrows()
     << "," << x.ncols() << "]";
  return os.str();
}

//! Describe Tensor3.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstTensor3View x) {
  ostringstream os;
  os << "Tensor3 [" << x.npages() << "," << x.nrows() << "," << x.ncols()
     << "]";
  return os.str();
}

//! Describe Matrix.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstMatrixView x) {
  ostringstream os;
  os << "Matrix [" << x.nrows() << "," << x.ncols() << "]";
  return os.str();
}

//! Describe Vector.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(ConstVectorView x) {
  ostringstream os;
  os << "Vector [" << x.nelem() << "]";
  return os.str();
}

//! Describe Scalar.
/*! 
  \param x  What to describe.
  \return   Output stream. */
String describe(const Numeric& x) {
  ostringstream os;
  os << "Scalar (" << x << ")";
  return os.str();
}
