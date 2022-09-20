/* Copyright (C) 2006-2019 Oliver Lemke <olemke@core-dump.info>

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

/**
    \file   matpack.h
    \author Stefan Buehler
*/

#ifndef matpack_h
#define matpack_h

#include "config.h"
#include "debug.h"

//--------------------< Set floating point type >--------------------
/** The type to use for all floating point numbers. You should never
    use float or double explicitly, unless you have a very good
    reason. Always use this type instead.  */
using Numeric = NUMERIC;

//--------------------< Set integer type >--------------------
/** The type to use for all integer numbers and indices. You should never
    use int, long, or size_t explicitly, unless you have a very good
    reason. Always use this type instead.  */
using Index = INDEX;

#endif  // matpackI_h
