/* Copyright (C) 2002-2008
   Claudia Emde <claudia.emde@dlr.de>
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


#include <stdexcept>
#include <iostream>
#include "mystring.h"
#include "exceptions.h"
#include "gridded_fields.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! Copy grids
/*!
  Copies the grids from the given GField to the current one.

  \param[in] gf Source gridded field
*/
void GField::copy_grids (const GField& gf)
{
  assert(gf.get_dim() == dim);

  for (Index i = 0; i < dim; i++)
    {
      switch (gf.get_grid_type(i))
        {
        case GRIDTYPE_NUMERIC:
          mgridtypes[i] = GRIDTYPE_NUMERIC;
          mnumericgrids[i] = gf.get_numeric_grid(i);
          mstringgrids[i].resize(0);
          break;
        case GRIDTYPE_STRING:
          mgridtypes[i] = GRIDTYPE_STRING;
          mstringgrids[i] = gf.get_string_grid(i);
          mnumericgrids[i].resize(0);
          break;
        }
    }
}


//! Get the size of a grid.
/*!
  Returns the size of grid i.

  \param[in]  i  Grid index.
  \return        Grid size.
*/
Index GField::get_grid_size (Index i) const
{
  Index ret = 0;
  assert (i < dim);
  switch (mgridtypes[i])
    {
    case GRIDTYPE_NUMERIC: ret = mnumericgrids[i].nelem(); break;
    case GRIDTYPE_STRING:  ret = mstringgrids[i].nelem(); break;
    }

  return ret;
}


//! Get a numeric grid.
/*!
  Returns the numeric grid with index i.

  Throws a runtime error if grid i is not of type Numeric.

  \param[in]  i  Grid index.
  \return        Numeric grid.
*/
ConstVectorView GField::get_numeric_grid (Index i) const
{
  assert (i < dim);
  if (mgridtypes[i] != GRIDTYPE_NUMERIC)
    {
      ostringstream os;
      if (mgridnames[i].length())
        os << mgridnames[i];
      else
        os << "Grid " << i;
      os << " is not a numeric grid.";
      throw runtime_error(os.str());
    }

  return (mnumericgrids[i]);
}


//! Get a numeric grid.
/*!
  Returns the numeric grid with index i.

  Throws a runtime error if grid i is not of type Numeric.

  \param[in]  i  Grid index.
  \return        Numeric grid.
*/
VectorView GField::get_numeric_grid (Index i)
{
  assert (i < dim);
  if (mgridtypes[i] != GRIDTYPE_NUMERIC)
    {
      ostringstream os;
      if (mgridnames[i].length())
        os << mgridnames[i];
      else
        os << "Grid " << i;
      os << " is not a numeric grid.";
      throw runtime_error(os.str());
    }

  return (mnumericgrids[i]);
}


//! Get a string grid.
/*!
  Returns the string grid with index i.

  Throws a runtime error if grid i is not of type String.

  \param[in]  i  Grid index.
  \return        String grid.
*/
const ArrayOfString& GField::get_string_grid (Index i) const
{
  assert (i < dim);
  if (mgridtypes[i] != GRIDTYPE_STRING)
    {
      ostringstream os;
      if (mgridnames[i].length())
        os << mgridnames[i];
      else
        os << "Grid " << i;
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
ArrayOfString& GField::get_string_grid (Index i)
{
  assert (i < dim);
  if (mgridtypes[i] != GRIDTYPE_STRING)
    {
      ostringstream os;
      if (mgridnames[i].length())
        os << mgridnames[i];
      else
        os << "Grid " << i;
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
void GField::set_grid (Index i, const Vector& g)
{
  assert (i < dim);
  mgridtypes[i] = GRIDTYPE_NUMERIC;
  mstringgrids[i].resize(0);
  mnumericgrids[i] = g;
}


//! Set a string grid.
/*!
  Sets grid i to the given grid.

  \param  i[in]  Grid index.
  \param  g[in]  New grid.
*/
void GField::set_grid (Index i, const ArrayOfString& g)
{
  assert (i < dim);
  mgridtypes[i] = GRIDTYPE_STRING;
  mnumericgrids[i].resize(0);
  mstringgrids[i] = g;
}


//! Output operator for GField
/*!
  Outputs the grids for the given GField.

  \param[in,out]  os  Output stream.
  \param[in]      gf  GField.
*/
ostream& operator<<(ostream& os, const GField& gf)
{
  if (gf.mname.size()) os << gf.mname << ":" << endl;

  for( Index i = 0; i < gf.dim; i++)
    {
      if (gf.mgridnames[i].size()) os << gf.mgridnames[i];
      else os << "Grid " << i;
      os << ": ";
      switch (gf.mgridtypes[i])
        {
        case GRIDTYPE_STRING:
          os << gf.mstringgrids[i];
          break;
        case GRIDTYPE_NUMERIC:
          os << gf.mnumericgrids[i];
          break;
        }
      os << endl;
    }

  return os;
}


//! Output operator for GField1
/*!
  Outputs the given GField1.

  \param[in,out]  os  Output stream.
  \param[in]      gf  GField1.
*/
ostream& operator<<(ostream& os, const GField1& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Vector&)gf << endl;
}


//! Output operator for GField2
/*!
  Outputs the given GField2.

  \param[in,out]  os  Output stream.
  \param[in]      gf  GField2.
*/
ostream& operator<<(ostream& os, const GField2& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Matrix&)gf;
}


//! Output operator for GField3
/*!
  Outputs the given GField3.

  \param[in,out]  os  Output stream.
  \param[in]      gf  GField3.
*/
ostream& operator<<(ostream& os, const GField3& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Tensor3&)gf;
}


//! Output operator for GField4
/*!
  Outputs the given GField4.

  \param[in,out]  os  Output stream.
  \param[in]      gf  GField4.
*/
ostream& operator<<(ostream& os, const GField4& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Tensor4&)gf;
}

