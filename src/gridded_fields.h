/* Copyright (C) 2002-2007 Claudia Emde <claudia.emde@dlr.de>

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
  \file   gridded_fields.h
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Wed Jun 26 17:48:29 2002

  \brief  Reading routines for gridded fields.

  This file contains reading routines for gridded fields. Gridded fields are
  needed to store moredimesional data together with the corresponding grids
  in the same variable. The datatype og a gridded field is always an array
  of tensors.

  For further description see AUG.

*/

#ifndef gridded_fields_h
#define gridded_fields_h

//! Contains a GriddedField3.
/*!
 A gridded field consists of a pressure grid vector, a latitude vector, a
 longitude vector and a Tensor3 for the data itself.

 \author Oliver Lemke
 \date   2003-06-19
 */
typedef struct {
      //! Pressure grid.
      Vector p_grid;
      //! Latitude grid.
      Vector lat_grid;
      //! Longitude grid.
      Vector lon_grid;
      //! Data.
      Tensor3 data;
} GriddedField3;

typedef Array<GriddedField3> ArrayOfGriddedField3;

ostream& operator<< (ostream &os, const GriddedField3 &gfield3);
ostream& operator<< (ostream &os, const ArrayOfGriddedField3 &agfield3);

void
check_gridded_tensor3 (ArrayOfTensor3 &gridded_tensor);

void
check_gridded_tensor6 (ArrayOfTensor6 &gridded_tensor);

void
read_gridded_tensor3 (const String& filename,
                      ArrayOfTensor3& gridded_tensor3
                     );

void
read_gridded_tensor6(
                     const String& filename,
                     ArrayOfTensor6& gridded_tensor6
                     );

#endif
