/* Copyright (C) 2002-2008 Claudia Emde <claudia.emde@dlr.de>

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

#include "matpackIV.h"
#include "mystring.h"


//! Contains a GriddedField3.
/*!
 A GriddedField3 consists of a pressure grid vector, a latitude vector, a
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

//! Contains a GriddedField4.
/*!
  A GriddedField4 can be used to store fields of several scalar
  variables on the same p-lat-lon grid. A typical example is to store
  several atmospheric variables from a circulation model, or a
  radiosonde profile.

  A GriddedField4 consists of a field name array, a pressure grid
  vector, a latitude vector, a longitude vector and a Tensor4 for the
  data itself. The dimensions in the data tensor are also sorted in
  the above order.
  
  \author Stefan Buehler
  \date   2007-07-25
 */
typedef struct {
      //! Field names
      ArrayOfString field_names;
      //! Pressure grid.
      Vector p_grid;
      //! Latitude grid.
      Vector lat_grid;
      //! Longitude grid.
      Vector lon_grid;
      //! Data.
      Tensor4 data;
} GriddedField4;

typedef Array< Array<GriddedField3> > ArrayOfArrayOfGriddedField3;
typedef Array<GriddedField3> ArrayOfGriddedField3;
typedef Array<GriddedField4> ArrayOfGriddedField4;

ostream& operator<< (ostream& os, const GriddedField3& gfield3);
ostream& operator<< (ostream& os, const ArrayOfGriddedField3& agfield3);

ostream& operator<< (ostream& os, const GriddedField4& gfield4);
ostream& operator<< (ostream& os, const ArrayOfGriddedField4& agfield4);

#endif
