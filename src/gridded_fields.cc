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
#include "artstime.h"
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
    set_grid_name(i, gf.get_grid_name(i));
    switch (gf.get_grid_type(i)) {
      case GRID_TYPE_NUMERIC:
        set_grid(i, gf.get_numeric_grid(i));
        break;
      case GRID_TYPE_STRING:
        set_grid(i, gf.get_string_grid(i));
        break;
      case GRID_TYPE_TIME:
        set_grid(i, gf.get_time_grid(i));
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
  ARTS_USER_ERROR_IF (mgridtypes[i] != GRID_TYPE_NUMERIC,
                      mname.length() ? var_string(mname, " Grid ") : var_string("Grid "),
                      mgridnames[i].length() ? var_string(mgridnames[i]) : var_string(i),
                      " is not a numeric grid.")
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
  ARTS_USER_ERROR_IF (mgridtypes[i] != GRID_TYPE_NUMERIC,
                      mname.length() ? var_string(mname, " Grid ") : var_string("Grid "),
                      mgridnames[i].length() ? var_string(mgridnames[i]) : var_string(i),
                      " is not a numeric grid.")
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
  ARTS_USER_ERROR_IF (mgridtypes[i] != GRID_TYPE_STRING,
                      mname.length() ? var_string(mname, " Grid ") : var_string("Grid "),
                      mgridnames[i].length() ? var_string(mgridnames[i]) : var_string(i),
                      " is not a string grid.")
  return mstringgrids[i];
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
  ARTS_USER_ERROR_IF (mgridtypes[i] != GRID_TYPE_STRING,
                      mname.length() ? var_string(mname, " Grid ") : var_string("Grid "),
                      mgridnames[i].length() ? var_string(mgridnames[i]) : var_string(i),
                      " is not a string grid.")
  return mstringgrids[i];
}

//! Get a time grid.
/*!
  Returns the time grid with index i.

  Throws a runtime error if grid i is not of type ArrayOfTime.

  \param[in]  i  Grid index.
  \return        Time grid.
*/
const ArrayOfTime& GriddedField::get_time_grid(Index i) const {
  ARTS_ASSERT(i < dim);
  ARTS_USER_ERROR_IF (mgridtypes[i] != GRID_TYPE_TIME,
                      mname.length() ? var_string(mname, " Grid ") : var_string("Grid "),
                      mgridnames[i].length() ? var_string(mgridnames[i]) : var_string(i),
                      " is not a time grid.")
  return mtimegrids[i];
}

//! Get a time grid.
/*!
  Returns the time grid with index i.

  Throws a runtime error if grid i is not of type ArrayOfTime.

  \param[in]  i  Grid index.
  \return        Time grid.
*/
ArrayOfTime& GriddedField::get_time_grid(Index i) {
  ARTS_ASSERT(i < dim);
  ARTS_USER_ERROR_IF (mgridtypes[i] != GRID_TYPE_TIME,
                      mname.length() ? var_string(mname, " Grid ") : var_string("Grid "),
                      mgridnames[i].length() ? var_string(mgridnames[i]) : var_string(i),
                      " is not a time grid.")
  return mtimegrids[i];
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
  mtimegrids[i].resize(0);
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
  mtimegrids[i].resize(0);
  mstringgrids[i] = g;
}

//! Set a time grid.
/*!
  Sets grid i to the given grid.

  \param[in] i Grid index.
  \param[in] g New grid.
*/
void GriddedField::set_grid(Index i, const ArrayOfTime& g) {
  ARTS_ASSERT(i < dim);
  mgridtypes[i] = GRID_TYPE_TIME;
  mnumericgrids[i].resize(0);
  mstringgrids[i].resize(0);
  mtimegrids[i] = g;
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
      case GRID_TYPE_TIME:
        os << gf.mtimegrids[i];
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

std::ostream& operator<<(std::ostream& os, const ArrayOfGriddedField1& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfGriddedField2& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfGriddedField3& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfGriddedField4& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfGriddedField5& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfGriddedField1& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfGriddedField2& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfGriddedField3& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
