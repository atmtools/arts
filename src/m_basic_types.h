/* Copyright (C) 2004-2012 Oliver Lemke <olemke@core-dump.info>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   m_basic_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2004-09-20
 
  \brief  Functions for straightforward operations on variables 
          of basic types.
*/

#ifndef M_BASIC_TYPES_H
#define M_BASIC_TYPES_H

#include <sstream>

#include "agenda_class.h"
#include "array.h"
#include "arts.h"
#include "exceptions.h"
#include "gridded_fields.h"
#include "matpackII.h"
#include "matpackVII.h"
#include "messages.h"
#include "mystring.h"
#include "workspace_ng.h"

/* To avoid redundant code preprocessor macros are used for code generation.
 *
 * Author: Oliver Lemke
 * Date:   2004-09-20
 */

////////////////////////////////////////////////////////////////////////
// The following template functions catch all calls to nelemGet,
// ncolsGet,etc. for data types (groups) which do not provide the
// requested attribute. A runtime error is thrown.

#define TMPL_NGET_GENERIC(what)                        \
  template <typename T>                                \
  void what##Get(Index&, const T&, const Verbosity&) { \
    ostringstream os;                                  \
    os << "The variable has no such attribute.\n";     \
    throw runtime_error(os.str());                     \
  }

TMPL_NGET_GENERIC(nelem)
TMPL_NGET_GENERIC(ncols)
TMPL_NGET_GENERIC(nrows)
TMPL_NGET_GENERIC(npages)
TMPL_NGET_GENERIC(nbooks)
TMPL_NGET_GENERIC(nshelves)
TMPL_NGET_GENERIC(nvitrines)
TMPL_NGET_GENERIC(nlibraries)

// Undefine the macro to make sure that it is never used anywhere else
#undef TMPL_NGET_GENERIC

#define TMPL_NGET_AGENDA(what)                                                 \
  void what##Get(Workspace& ws _U_, Index&, const Agenda&, const Verbosity&) { \
    ostringstream os;                                                          \
    os << "The variable has no such attribute.\n";                             \
    throw runtime_error(os.str());                                             \
  }

TMPL_NGET_AGENDA(nelem)
TMPL_NGET_AGENDA(ncols)
TMPL_NGET_AGENDA(nrows)
TMPL_NGET_AGENDA(npages)
TMPL_NGET_AGENDA(nbooks)
TMPL_NGET_AGENDA(nshelves)
TMPL_NGET_AGENDA(nvitrines)
TMPL_NGET_AGENDA(nlibraries)

// Undefine the macro to make sure that it is never used anywhere else
#undef TMPL_NGET_AGENDA

template <typename T>
void IndexSetToLast(Index&, const T&, const Verbosity&) {
  ostringstream os;
  os << "The variable has no such attribute.\n";
  throw runtime_error(os.str());
}

////////////////////////////////////////////////////////////////////////
// The following functions are special implementations of the template
// functions above. They set the corresponding workspace variable to the
// value of the requested attribute.

#define NGET_GENERIC(what, type)                                 \
  void what##Get(Index& what, const type& x, const Verbosity&) { \
    what = x.what();                                             \
  }

#define SET_TO_LAST_GENERIC(type)                                  \
  void IndexSetToLast(Index& i, const type& x, const Verbosity&) { \
    i = x.nelem() - 1;                                             \
  }

// If you add a group here, add it is also to SET_TO_LAST
NGET_GENERIC(nelem, Vector)
NGET_GENERIC(nelem, ArrayOfIndex)
NGET_GENERIC(nelem, ArrayOfArrayOfIndex)
NGET_GENERIC(nelem, ArrayOfString)
NGET_GENERIC(nelem, ArrayOfVector)
NGET_GENERIC(nelem, ArrayOfArrayOfVector)
NGET_GENERIC(nelem, ArrayOfMatrix)
NGET_GENERIC(nelem, ArrayOfArrayOfMatrix)
NGET_GENERIC(nelem, ArrayOfSparse)
NGET_GENERIC(nelem, ArrayOfTensor3)
NGET_GENERIC(nelem, ArrayOfArrayOfTensor3)
NGET_GENERIC(nelem, ArrayOfTensor4)
NGET_GENERIC(nelem, ArrayOfTensor5)
NGET_GENERIC(nelem, ArrayOfTensor6)
NGET_GENERIC(nelem, ArrayOfArrayOfTensor6)
NGET_GENERIC(nelem, ArrayOfTensor7)
NGET_GENERIC(nelem, ArrayOfArrayOfSpeciesTag)
NGET_GENERIC(nelem, ArrayOfSingleScatteringData)
NGET_GENERIC(nelem, ArrayOfScatteringMetaData)
NGET_GENERIC(nelem, ArrayOfArrayOfSingleScatteringData)
NGET_GENERIC(nelem, ArrayOfArrayOfScatteringMetaData)
NGET_GENERIC(nelem, ArrayOfGriddedField1)
NGET_GENERIC(nelem, ArrayOfGriddedField2)
NGET_GENERIC(nelem, ArrayOfGriddedField3)
NGET_GENERIC(nelem, ArrayOfGriddedField4)
NGET_GENERIC(nelem, ArrayOfArrayOfGriddedField1)
NGET_GENERIC(nelem, ArrayOfArrayOfGriddedField2)
NGET_GENERIC(nelem, ArrayOfArrayOfGriddedField3)
NGET_GENERIC(nelem, ArrayOfRetrievalQuantity)
NGET_GENERIC(nelem, ArrayOfAbsorptionLines)
NGET_GENERIC(nelem, ArrayOfArrayOfAbsorptionLines)
NGET_GENERIC(nelem, ArrayOfPropagationMatrix)
NGET_GENERIC(nelem, ArrayOfArrayOfPropagationMatrix)
NGET_GENERIC(nelem, ArrayOfTransmissionMatrix)
NGET_GENERIC(nelem, ArrayOfArrayOfTransmissionMatrix)

SET_TO_LAST_GENERIC(Vector)
SET_TO_LAST_GENERIC(ArrayOfIndex)
SET_TO_LAST_GENERIC(ArrayOfArrayOfIndex)
SET_TO_LAST_GENERIC(ArrayOfString)
SET_TO_LAST_GENERIC(ArrayOfVector)
SET_TO_LAST_GENERIC(ArrayOfArrayOfVector)
SET_TO_LAST_GENERIC(ArrayOfMatrix)
SET_TO_LAST_GENERIC(ArrayOfArrayOfMatrix)
SET_TO_LAST_GENERIC(ArrayOfSparse)
SET_TO_LAST_GENERIC(ArrayOfTensor3)
SET_TO_LAST_GENERIC(ArrayOfArrayOfTensor3)
SET_TO_LAST_GENERIC(ArrayOfTensor4)
SET_TO_LAST_GENERIC(ArrayOfTensor5)
SET_TO_LAST_GENERIC(ArrayOfTensor6)
SET_TO_LAST_GENERIC(ArrayOfArrayOfTensor6)
SET_TO_LAST_GENERIC(ArrayOfTensor7)
SET_TO_LAST_GENERIC(ArrayOfArrayOfSpeciesTag)
SET_TO_LAST_GENERIC(ArrayOfSingleScatteringData)
SET_TO_LAST_GENERIC(ArrayOfScatteringMetaData)
SET_TO_LAST_GENERIC(ArrayOfArrayOfSingleScatteringData)
SET_TO_LAST_GENERIC(ArrayOfArrayOfScatteringMetaData)
SET_TO_LAST_GENERIC(ArrayOfGriddedField1)
SET_TO_LAST_GENERIC(ArrayOfGriddedField2)
SET_TO_LAST_GENERIC(ArrayOfGriddedField3)
SET_TO_LAST_GENERIC(ArrayOfGriddedField4)
SET_TO_LAST_GENERIC(ArrayOfArrayOfGriddedField1)
SET_TO_LAST_GENERIC(ArrayOfArrayOfGriddedField2)
SET_TO_LAST_GENERIC(ArrayOfArrayOfGriddedField3)
SET_TO_LAST_GENERIC(ArrayOfRetrievalQuantity)
SET_TO_LAST_GENERIC(ArrayOfAbsorptionLines)
SET_TO_LAST_GENERIC(ArrayOfArrayOfAbsorptionLines)
SET_TO_LAST_GENERIC(ArrayOfPropagationMatrix)
SET_TO_LAST_GENERIC(ArrayOfArrayOfPropagationMatrix)
SET_TO_LAST_GENERIC(ArrayOfTransmissionMatrix)
SET_TO_LAST_GENERIC(ArrayOfArrayOfTransmissionMatrix)

NGET_GENERIC(ncols, Matrix)
NGET_GENERIC(ncols, Sparse)
NGET_GENERIC(ncols, Tensor3)
NGET_GENERIC(ncols, Tensor4)
NGET_GENERIC(ncols, Tensor5)
NGET_GENERIC(ncols, Tensor6)
NGET_GENERIC(ncols, Tensor7)

NGET_GENERIC(nrows, Matrix)
NGET_GENERIC(nrows, Sparse)
NGET_GENERIC(nrows, Tensor3)
NGET_GENERIC(nrows, Tensor4)
NGET_GENERIC(nrows, Tensor5)
NGET_GENERIC(nrows, Tensor6)
NGET_GENERIC(nrows, Tensor7)

NGET_GENERIC(npages, Tensor3)
NGET_GENERIC(npages, Tensor4)
NGET_GENERIC(npages, Tensor5)
NGET_GENERIC(npages, Tensor6)
NGET_GENERIC(npages, Tensor7)

NGET_GENERIC(nbooks, Tensor4)
NGET_GENERIC(nbooks, Tensor5)
NGET_GENERIC(nbooks, Tensor6)
NGET_GENERIC(nbooks, Tensor7)

NGET_GENERIC(nshelves, Tensor5)
NGET_GENERIC(nshelves, Tensor6)
NGET_GENERIC(nshelves, Tensor7)

NGET_GENERIC(nvitrines, Tensor6)
NGET_GENERIC(nvitrines, Tensor7)

NGET_GENERIC(nlibraries, Tensor7)

// Undefine the macros to make sure that it is never used anywhere else
#undef NGET_GENERIC
#undef SET_TO_LAST_GENERIC

void nelemGet(Workspace& /* ws */,
              Index& nelem,
              const ArrayOfAgenda& x,
              const Verbosity&) {
  nelem = x.nelem();
}

void IndexSetToLast(Workspace& /* ws */,
                    Index& nelem,
                    const ArrayOfAgenda& x,
                    const Verbosity&) {
  nelem = x.nelem() - 1;
}

#endif /* M_BASIC_TYPES_H */
