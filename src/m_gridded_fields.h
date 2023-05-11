/* Copyright (C) 2002-2012 Oliver Lemke <olemke@core-dump.info>

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

/*!
  \file   m_gridded_fields.h
  \author Lukas Kluft  <lukas.kluft@gmail.com>
  \date   2017-02-01

  \brief  Implementation of GriddedField workspace methods.

  This file contains the implementation of the supergeneric methods.
*/

#ifndef m_gridded_fields_h
#define m_gridded_fields_h

#include "gridded_fields.h"
#include "mystring.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void GriddedFieldGetName(  // WS Generic Output:
    String& name,
    // WS Generic Input:
    const T& gf) {
  // Return the name of the given GriddedField.
  name = gf.get_name();
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ArrayOfGriddedFieldGetNames(  // WS Generic Output:
    ArrayOfString& names,
    // WS Generic Input:
    const Array<T>& aogf) {
  // Return the name of the given GriddedField.
  names.resize(aogf.nelem());
  for (Index i = 0; i < aogf.nelem(); i++) {
    names[i] = aogf[i].get_name();
  }
}

#endif  // m_gridded_fields_h
