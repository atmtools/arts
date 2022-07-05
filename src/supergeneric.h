/* Copyright (C) 2002-2012
   Stefan Buehler <sbuehler@ltu.se>

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
  \file   supergeneric.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Thu Jun 13 14:44:02 2002
  
  \brief  Declarations for supergeneric methods.

  The class Any can be used to mark supergenerio IO.
*/

#ifndef supergeneric_h
#define supergeneric_h

#include <ostream>

//! A placeholder for any type.
/*!  Used to mark supergeneric methods in file
  methods.cc.
*/
class Any {
  // Nothing to do here.
  friend std::ostream& operator<<(std::ostream&, const Any&) {}
};

#endif  // supergeneric_h
