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

#include "array.h"
#include "arts.h"
#include "exceptions.h"
#include "matpackII.h"
#include "matpackVII.h"
#include "gridded_fields.h"
#include "mystring.h"
#include "messages.h"
#include "workspace_ng.h"
#include "agenda_class.h"

/* To avoid redundant code preprocessor macros are used for code generation.
 *
 * Author: Oliver Lemke
 * Date:   2004-09-20
 */

////////////////////////////////////////////////////////////////////////
// The following template functions catch all calls to nelemGet,
// ncolsGet,etc. for data types (groups) which do not provide the
// requested attribute. A runtime error is thrown.

#define TMPL_NGET_GENERIC(what) \
  template <typename T> \
  void what##Get( Index&, \
                  const T&, \
                  const Verbosity&) \
  { \
    ostringstream os; \
    os << "The variable has no such attribute.\n"; \
    throw runtime_error(os.str()); \
 }

TMPL_NGET_GENERIC (nelem)
TMPL_NGET_GENERIC (ncols)
TMPL_NGET_GENERIC (nrows)
TMPL_NGET_GENERIC (npages)
TMPL_NGET_GENERIC (nbooks)
TMPL_NGET_GENERIC (nshelves)
TMPL_NGET_GENERIC (nvitrines)
TMPL_NGET_GENERIC (nlibraries)

// Undefine the macro to make sure that it is never used anywhere else
#undef TMPL_NGET_GENERIC

#define TMPL_NGET_AGENDA(what) \
  void what##Get( Workspace& ws _U_, \
                  Index&, \
                  const Agenda&, \
                  const Verbosity&) \
  { \
    ostringstream os; \
    os << "The variable has no such attribute.\n"; \
    throw runtime_error(os.str()); \
 }

TMPL_NGET_AGENDA (nelem)
TMPL_NGET_AGENDA (ncols)
TMPL_NGET_AGENDA (nrows)
TMPL_NGET_AGENDA (npages)
TMPL_NGET_AGENDA (nbooks)
TMPL_NGET_AGENDA (nshelves)
TMPL_NGET_AGENDA (nvitrines)
TMPL_NGET_AGENDA (nlibraries)

// Undefine the macro to make sure that it is never used anywhere else
#undef TMPL_NGET_AGENDA


////////////////////////////////////////////////////////////////////////
// The following functions are special implementations of the template
// functions above. They set the corresponding workspace variable to the
// value of the requested attribute.

#define NGET_GENERIC(what, type) \
  void what##Get(Index&    what, \
                 const type&   x, \
                 const Verbosity&) \
  { \
    what = x.what (); \
  }

NGET_GENERIC (nelem, Vector)
NGET_GENERIC (nelem, ArrayOfIndex)
NGET_GENERIC (nelem, ArrayOfArrayOfIndex)
NGET_GENERIC (nelem, ArrayOfString)
NGET_GENERIC (nelem, ArrayOfVector)
NGET_GENERIC (nelem, ArrayOfArrayOfVector)
NGET_GENERIC (nelem, ArrayOfMatrix)
NGET_GENERIC (nelem, ArrayOfArrayOfMatrix)
NGET_GENERIC (nelem, ArrayOfSparse)
NGET_GENERIC (nelem, ArrayOfTensor3)
NGET_GENERIC (nelem, ArrayOfArrayOfTensor3)
NGET_GENERIC (nelem, ArrayOfTensor4)
NGET_GENERIC (nelem, ArrayOfTensor5)
NGET_GENERIC (nelem, ArrayOfTensor6)
NGET_GENERIC (nelem, ArrayOfArrayOfTensor6)
NGET_GENERIC (nelem, ArrayOfTensor7)
NGET_GENERIC (nelem, ArrayOfLineMixingRecord)
NGET_GENERIC (nelem, ArrayOfLineRecord)
NGET_GENERIC (nelem, ArrayOfArrayOfLineRecord)
NGET_GENERIC (nelem, ArrayOfLineshapeSpec)
NGET_GENERIC (nelem, ArrayOfArrayOfSpeciesTag)
NGET_GENERIC (nelem, ArrayOfSingleScatteringData)
NGET_GENERIC (nelem, ArrayOfScatteringMetaData)
NGET_GENERIC (nelem, ArrayOfArrayOfSingleScatteringData)
NGET_GENERIC (nelem, ArrayOfArrayOfScatteringMetaData)
NGET_GENERIC (nelem, ArrayOfGriddedField1)
NGET_GENERIC (nelem, ArrayOfGriddedField2)
NGET_GENERIC (nelem, ArrayOfGriddedField3)
NGET_GENERIC (nelem, ArrayOfGriddedField4)
NGET_GENERIC (nelem, ArrayOfArrayOfGriddedField1)
NGET_GENERIC (nelem, ArrayOfArrayOfGriddedField2)
NGET_GENERIC (nelem, ArrayOfArrayOfGriddedField3)
NGET_GENERIC (nelem, ArrayOfArrayOfLineMixingRecord)
NGET_GENERIC (nelem, ArrayOfRetrievalQuantity)

NGET_GENERIC (ncols, Matrix)
NGET_GENERIC (ncols, Sparse)
NGET_GENERIC (ncols, Tensor3)
NGET_GENERIC (ncols, Tensor4)
NGET_GENERIC (ncols, Tensor5)
NGET_GENERIC (ncols, Tensor6)
NGET_GENERIC (ncols, Tensor7)

NGET_GENERIC (nrows, Matrix)
NGET_GENERIC (nrows, Sparse)
NGET_GENERIC (nrows, Tensor3)
NGET_GENERIC (nrows, Tensor4)
NGET_GENERIC (nrows, Tensor5)
NGET_GENERIC (nrows, Tensor6)
NGET_GENERIC (nrows, Tensor7)

NGET_GENERIC (npages, Tensor3)
NGET_GENERIC (npages, Tensor4)
NGET_GENERIC (npages, Tensor5)
NGET_GENERIC (npages, Tensor6)
NGET_GENERIC (npages, Tensor7)

NGET_GENERIC (nbooks, Tensor4)
NGET_GENERIC (nbooks, Tensor5)
NGET_GENERIC (nbooks, Tensor6)
NGET_GENERIC (nbooks, Tensor7)

NGET_GENERIC (nshelves, Tensor5)
NGET_GENERIC (nshelves, Tensor6)
NGET_GENERIC (nshelves, Tensor7)

NGET_GENERIC (nvitrines, Tensor6)
NGET_GENERIC (nvitrines, Tensor7)

NGET_GENERIC (nlibraries, Tensor7)

// Undefine the macro to make sure that it is never used anywhere else
#undef NGET_GENERIC



#endif /* M_BASIC_TYPES_H */

