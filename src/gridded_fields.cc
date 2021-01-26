/* Copyright (C) 2008-2012 Oliver Lemke <olemke@core-dump.info>
     
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
   USA. 
*/

/*!
  \file   gridded_fields.cc
  \author Oliver Lemke  <olemke@core-dump.info>
  \date   2008-06-24
  
  \brief  Implementation of gridded fields.
  
  This file contains the implemenation gridded fields. Gridded fields are
  needed to store moredimesional data together with the corresponding grids 
  in the same variable.

  For further description see AUG.

*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "gridded_fields.h"
#include <iostream>
#include <stdexcept>
#include "exceptions.h"
#include "mystring.h"

using std::endl;
using std::ostream;
using std::ostringstream;
using std::runtime_error;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! Copy grids
/*!
  Copies the grids from the given GriddedField to the current one.

  \param[in] gf Source gridded field
*/
void GriddedField::copy_grids(const GriddedField& gf) {
  ARTS_ASSERT(gf.get_dim() == dim);

  for (Index i = 0; i < dim; i++) {
    switch (gf.get_grid_type(i)) {
      case GRID_TYPE_NUMERIC:
        mgridtypes[i] = GRID_TYPE_NUMERIC;
        mnumericgrids[i] = gf.get_numeric_grid(i);
        mstringgrids[i].resize(0);
        break;
      case GRID_TYPE_STRING:
        mgridtypes[i] = GRID_TYPE_STRING;
        mstringgrids[i] = gf.get_string_grid(i);
        mnumericgrids[i].resize(0);
        break;
    }
  }
}

//! Get a numeric grid.
/*!
  Returns the numeric grid with index i.

  Throws a runtime error if grid i is not of type Numeric.

  \param[in]  i  Grid index.
  \return        Numeric grid.
*/
const Vector& GriddedField::get_numeric_grid(Index i) const {
  ARTS_ASSERT(i < dim);
  if (mgridtypes[i] != GRID_TYPE_NUMERIC) {
    ostringstream os;

    if (mname.length()) os << mname << " ";

    os << "Grid ";
    if (mgridnames[i].length())
      os << mgridnames[i];
    else
      os << i;
    os << " is not a numeric grid.";

    throw runtime_error(os.str());
  }

  return mnumericgrids[i];
}

//! Get a numeric grid.
/*!
  Returns the numeric grid with index i.

  Throws a runtime error if grid i is not of type Numeric.

  \param[in]  i  Grid index.
  \return        Numeric grid.
*/
Vector& GriddedField::get_numeric_grid(Index i) {
  ARTS_ASSERT(i < dim);
  if (mgridtypes[i] != GRID_TYPE_NUMERIC) {
    ostringstream os;

    if (mname.length()) os << mname << " ";

    os << "Grid ";
    if (mgridnames[i].length())
      os << mgridnames[i];
    else
      os << i;
    os << " is not a numeric grid.";

    throw runtime_error(os.str());
  }

  return mnumericgrids[i];
}

//! Get a string grid.
/*!
  Returns the string grid with index i.

  Throws a runtime error if grid i is not of type String.

  \param[in]  i  Grid index.
  \return        String grid.
*/
const ArrayOfString& GriddedField::get_string_grid(Index i) const {
  ARTS_ASSERT(i < dim);
  if (mgridtypes[i] != GRID_TYPE_STRING) {
    ostringstream os;

    if (mname.length()) os << mname << " ";

    os << "Grid ";
    if (mgridnames[i].length())
      os << mgridnames[i];
    else
      os << i;
    os << " is not a string grid.";

    throw runtime_error(os.str());
  }

  return (mstringgrids[i]);
}

//! Get a string grid.
/*!
  Returns the string grid with index i.

  Throws a runtime error if grid i is not of type String.

  \param[in]  i  Grid index.
  \return        String grid.
*/
ArrayOfString& GriddedField::get_string_grid(Index i) {
  ARTS_ASSERT(i < dim);
  if (mgridtypes[i] != GRID_TYPE_STRING) {
    ostringstream os;

    if (mname.length()) os << mname << " ";

    os << "Grid ";
    if (mgridnames[i].length())
      os << mgridnames[i];
    else
      os << i;
    os << " is not a string grid.";

    throw runtime_error(os.str());
  }

  return (mstringgrids[i]);
}

//! Set a numeric grid.
/*!
  Sets grid i to the given grid.

  \param[in]  i  Grid index.
  \param[in]  g  New grid.
*/
void GriddedField::set_grid(Index i, const Vector& g) {
  ARTS_ASSERT(i < dim);
  mgridtypes[i] = GRID_TYPE_NUMERIC;
  mstringgrids[i].resize(0);
  mnumericgrids[i] = g;
}

//! Set a string grid.
/*!
  Sets grid i to the given grid.

  \param[in] i Grid index.
  \param[in] g New grid.
*/
void GriddedField::set_grid(Index i, const ArrayOfString& g) {
  ARTS_ASSERT(i < dim);
  mgridtypes[i] = GRID_TYPE_STRING;
  mnumericgrids[i].resize(0);
  mstringgrids[i] = g;
}

//! Output operator for GriddedField
/*!
 Outputs the grids for the given GriddedField.

 \param[in,out]  os  Output stream.
 \param[in]      gf  GriddedField.
 */
ostream& operator<<(ostream& os, const GriddedField& gf) {
  if (gf.mname.size()) os << gf.mname << ":" << endl;

  for (Index i = 0; i < gf.dim; i++) {
    if (gf.mgridnames[i].size())
      os << gf.mgridnames[i];
    else
      os << "Grid " << i;
    os << ": ";
    switch (gf.mgridtypes[i]) {
      case GRID_TYPE_STRING:
        os << gf.mstringgrids[i];
        break;
      case GRID_TYPE_NUMERIC:
        os << gf.mnumericgrids[i];
        break;
    }
    os << endl;
  }

  return os;
}

//! See GriddedField::operator<<
ostream& operator<<(ostream& os, const GriddedField1& gf) {
  return os << *((GriddedField*)&gf) << "Data:" << endl << gf.data << endl;
}

//! See GriddedField::operator<<
ostream& operator<<(ostream& os, const GriddedField2& gf) {
  return os << *((GriddedField*)&gf) << "Data:" << endl << gf.data << endl;
}

//! See GriddedField::operator<<
ostream& operator<<(ostream& os, const GriddedField3& gf) {
  return os << *((GriddedField*)&gf) << "Data:" << endl << gf.data << endl;
}

//! See GriddedField::operator<<
ostream& operator<<(ostream& os, const GriddedField4& gf) {
  return os << *((GriddedField*)&gf) << "Data:" << endl << gf.data << endl;
}

//! See GriddedField::operator<<
ostream& operator<<(ostream& os, const GriddedField5& gf) {
  return os << *((GriddedField*)&gf) << "Data:" << endl << gf.data << endl;
}

//! See GriddedField::operator<<
ostream& operator<<(ostream& os, const GriddedField6& gf) {
  return os << *((GriddedField*)&gf) << "Data:" << endl << gf.data << endl;
}
