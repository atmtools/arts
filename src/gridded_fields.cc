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
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Wed Jun 26 17:48:29 2002
  
  \brief  Reading routines for gridded fields.
  
  This file contains reading routines for gridded fields. Gridded fields are
  needed to store moredimesional data together with the corresponding grids 
  in the same variable. The datatype og a gridded field is always an array
  of tensors.

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


ConstVectorView GField::get_numeric_grid (Index i) const
{
  assert (i < dim);
  assert (mgridtypes[i] == GRIDTYPE_NUMERIC);

  return (mnumericgrids[i]);
}


const ArrayOfString& GField::get_string_grid (Index i) const
{
  assert (i < dim);
  assert (mgridtypes[i] == GRIDTYPE_STRING);

  return (mstringgrids[i]);
}


void GField::set_grid (Index i, const Vector& g)
{
  assert (i < dim);
  mgridtypes[i] = GRIDTYPE_NUMERIC;
  mstringgrids[i].resize(0);
  mnumericgrids[i] = g;
}


void GField::set_grid (Index i, const ArrayOfString& g)
{
  assert (i < dim);
  mgridtypes[i] = GRIDTYPE_STRING;
  mnumericgrids[i].resize(0);
  mstringgrids[i] = g;
}


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


ostream& operator<<(ostream& os, const GField1& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Vector&)gf << endl;
}


ostream& operator<<(ostream& os, const GField2& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Matrix&)gf;
}


ostream& operator<<(ostream& os, const GField3& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Tensor3&)gf;
}


ostream& operator<<(ostream& os, const GField4& gf)
{
  os << (GField&)gf;
  return os << "Data:" << endl << (Tensor4&)gf;
}


ostream& operator<< (ostream& os, const GriddedField3& /*gfield3*/)
{
  os << "GriddedField3: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream& os, const ArrayOfGriddedField3& /*agfield3*/)
{
  os << "ArrayOfGriddedField3: Output operator not implemented";
  return os;
}

ostream& operator<< (ostream& os, const GriddedField4& /*gfield3*/)
{
  os << "GriddedField4: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream& os, const ArrayOfGriddedField4& /*agfield3*/)
{
  os << "ArrayOfGriddedField4: Output operator not implemented";
  return os;
}

