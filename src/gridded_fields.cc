/* Copyright (C) 2002-2007
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
#include "xml_io.h"
#include "array.h"
#include "matpackIII.h"
#include "matpackVI.h"
#include "gridded_fields.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

ostream& operator<< (ostream &os, const GriddedField3 &/*gfield3*/)
{
  os << "GriddedField3: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream &os, const ArrayOfGriddedField3 &/*agfield3*/)
{
  os << "ArrayOfGriddedField3: Output operator not implemented";
  return os;
}

ostream& operator<< (ostream &os, const GriddedField4 &/*gfield3*/)
{
  os << "GriddedField4: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream &os, const ArrayOfGriddedField4 &/*agfield3*/)
{
  os << "ArrayOfGriddedField4: Output operator not implemented";
  return os;
}

