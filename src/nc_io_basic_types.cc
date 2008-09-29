/* Copyright (C) 2003-2008 Oliver Lemke <olemke@core-dump.info>

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


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-26

  \brief This file contains basic functions to handle NetCDF data files.

*/

#include <netcdfcpp.h>
#include "arts.h"
#include "nc_io.h"
//#include "nc_io_private.h"
#include "nc_io_types.h"

//=== Vector ==========================================================

//! Writes Vector to NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param vector  Vector
*/
void
nc_write_to_file (NcFile &ncf _U_,
                  const Vector& v _U_)
{
  NcDim *dim = ncf.add_dim ("nelem", v.nelem());

  NcVar *data = ncf.add_var ("Vector", ncDouble, dim);

  for (Index i = 0; i < v.nelem(); i++)
    {
      double n = v[i];
      data->put (&n, 1);
    }

}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

