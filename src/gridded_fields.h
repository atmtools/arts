/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
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

void
check_gridded_tensor3 (ArrayOfTensor3 &gridded_tensor);

void
check_gridded_tensor6 (ArrayOfTensor6 &gridded_tensor);

void 
read_gridded_tensor3(
                     const String& filename,
                     ArrayOfTensor3& gridded_tensor3
                     );

void 
read_gridded_tensor6(
                     const String& filename,
                     ArrayOfTensor6& gridded_tensor6
                     );

#endif
