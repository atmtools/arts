/* Copyright (C) 2003 Oliver Lemke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

/*!
  \file   m_general.h
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2003-07-24
 
  \brief  Template functions for general supergeneric ws methods.
 
*/
 
//! Print
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Oliver Lemke
   \date   2003-07-24
*/

#ifndef m_general_h
#define m_general_h

#include <iostream>
#include <sys/times.h>
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "mystring.h"


#define SWITCH_OUTPUT(x,y) switch (x) { \
case 0: out0 << y;break; \
case 1: out1 << y;break; \
case 2: out2 << y;break; \
case 3: out3 << y;break; \
default: throw runtime_error ("Output level must have value from 0-3"); \
}


typedef struct {
  struct tms cputime;
  clock_t realtime;
} Timer;

//! Print
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Oliver Lemke
   \date   2003-07-24
*/
template<typename T> void
Print(
      // WS Generic Input:
      const T&        x,
      // WS Generic Input Names:
      const String&   x_name,
      // Keywords:
      const Index    level )
{
  SWITCH_OUTPUT (level, "  *" << x_name << "*:\n" << x << '\n')
}

void
Print(
      // WS Generic Input:
      const ArrayOfGridPos&   x,
      // WS Generic Input Names:
      const String&           x_name,
      // Keywords:
      const Index              level );

void
Print(
        // WS Generic Input:
        const ArrayOfIndex&   x,
        // WS Generic Input Names:
        const String&         x_name,
        // Keywords:
        const Index           level );

void
Print(
        // WS Generic Input:
        const ArrayOfString&   x,
        // WS Generic Input Names:
        const String&          x_name,
        // Keywords:
        const Index            level );

void
Print(
      // WS Generic Input:
      const Ppath&    ppath,
      // WS Generic Input Names:
      const String&   x_name,
      // Keywords:
      const Index     level );

void
Print(
      // WS Generic Input:
      const Timer&   x,
      // WS Generic Input Names:
      const String&  x_name,
      // Keywords:
      const Index    level);

#endif /* m_general_h */

