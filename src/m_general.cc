/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_general.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions of a general and overall character.

  This file is for general functions that do not fit in any other "m_"-file.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>

#include "arts.h"

#include "m_general.h"
#include "array.h"
#include "check_input.h"
#include "messages.h"
#include "mystring.h"

#include "math_funcs.h"
#include "make_vector.h"
#include "sensor.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

void Print(
        // WS Generic Input:
        const ArrayOfIndex&   x,
        // WS Generic Input Names:
        const String&         x_name,
        // Keywords:
        const Index&    level )
{
  ostringstream os;
  cout << "  *" << x_name <<"* =";
  SWITCH_OUTPUT (level, "  *" << x_name << "*:\n")
  for( Index i=0; i<x.nelem(); i++ )
    os << x[i];
  SWITCH_OUTPUT (level, os.str () << '\n')
}


void Print(
        // WS Generic Input:
        const ArrayOfString&   x,
        // WS Generic Input Names:
        const String&          x_name,
        // Keywords:
        const Index&    level )
{
  ostringstream os;
  cout << "  *" << x_name <<"* =";
  SWITCH_OUTPUT (level, "  *" << x_name << "*:\n")
  for( Index i=0; i<x.nelem(); i++ )
    os << x[i] << '\n';
  SWITCH_OUTPUT (level, os.str ())
}


void
Print(
      // WS Generic Input:
      const Ppath&    x,
      // WS Generic Input Names:
      const String&   x_name,
      // Keywords:
      const Index&    level )
{
  SWITCH_OUTPUT (level, "  The fields of *" << x_name << "*:\n")
  Print( x.dim, "dim", level );
  Print( x.np, "np", level );
  Print( x.refraction, "refraction", level );
  Print( x.method, "method", level );
  Print( x.constant, "constant", level );
  Print( x.pos, "pos", level );
  Print( x.z, "z", level );
  Print( x.l_step, "l_step", level );
  ArrayOfGridPosPrint( x.gp_p, "gp_p" );
  if( x.dim >= 2 )
    ArrayOfGridPosPrint( x.gp_lat, "gp_lat" );
  if( x.dim == 3 )
    ArrayOfGridPosPrint( x.gp_lon, "gp_lon" );
  Print( x.los, "los", level );
  Print( x.background, "background", level );
  if( x.tan_pos.nelem() )
    Print( x.tan_pos, "tan_pos", level );
  if( x.geom_tan_pos.nelem() )
    Print( x.geom_tan_pos, "geom_tan_pos", level );
}


//! Exit
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void Exit()
{
  out1 << "  Forced exit.\n";
  arts_exit (0);
}


//! Test
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-05-15
*/
void Test( )
{
  // This function can be used to test stuff.

}



