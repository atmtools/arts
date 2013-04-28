/* Copyright (C) 2013
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
   USA. */

/**
  \file  debug.h

  Helper macros for debugging.

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef debug_h
#define debug_h

#ifndef NDEBUG

#include <iostream>


// Use this macro around function parameter names and variable definitions
// which are only used in assertions
#define DEBUG_ONLY(...) __VA_ARGS__

// Use this macro to output a counter value everytime a
// certain place is reached
#define DEBUG_COUNTER(n) \
{ \
  static Index n = 0; \
  std::cerr << "DBG: " << #n << ": " << ++n << std::endl; \
}

// Print expression for debug
#define DEBUG_PRINT(e) \
{ \
  std::cerr << "DBG: " << #e << ": " << (e) << std::endl; \
}

#else

#define DEBUG_ONLY(...)

#define DEBUG_COUNTER(n)

#define DEBUG_PRINT(e)

#endif /* NDEBUG */

#endif /* debug_h */


