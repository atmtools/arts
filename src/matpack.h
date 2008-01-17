/* Copyright (C) 2006-2007 Oliver Lemke <olemke@core-dump.info>

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

#ifndef matpack_h
#define matpack_h

#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Please run ./configure in the top arts directory before compiling."
#endif

#ifdef HAVE_NAMESPACES
  // We need those to support ansi-compliant compilers (gcc-3x)
  namespace std {}
  using namespace std;
#endif

//--------------------< Set floating point type >--------------------
/** The type to use for all floating point numbers. You should never
    use float or double explicitly, unless you have a very good
    reason. Always use this type instead.  */
typedef NUMERIC Numeric;

//--------------------< Set integer type >--------------------
/** The type to use for all integer numbers and indices. You should never
    use int, long, or size_t explicitly, unless you have a very good
    reason. Always use this type instead.  */
typedef INDEX Index;

#endif    // matpackI_h
