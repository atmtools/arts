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
#include <unistd.h>

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
        const ArrayOfGridPos&   x,
        // WS Generic Input Names:
        const String&           x_name,
        // Keywords:
        const Index             level )
{
  ostringstream os;
  SWITCH_OUTPUT (level, "  *" << x_name << "*:\n")
  for( Index i=0; i<x.nelem(); i++ )
    os << "     " << x[i].idx << "  " << x[i].fd[0] << "  " << x[i].fd[1]
         << "\n";
  SWITCH_OUTPUT (level, os.str ())
}


void Print(
        // WS Generic Input:
        const ArrayOfIndex&   x,
        // WS Generic Input Names:
        const String&         x_name,
        // Keywords:
        const Index           level )
{
  ostringstream os;
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
        const Index            level )
{
  ostringstream os;
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
      const Index     level )
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
  Print( x.gp_p, "gp_p", level );
  if( x.dim >= 2 )
    Print( x.gp_lat, "gp_lat", level );
  if( x.dim == 3 )
    Print( x.gp_lon, "gp_lon", level );
  Print( x.los, "los", level );
  Print( x.background, "background", level );
  if( x.tan_pos.nelem() )
    Print( x.tan_pos, "tan_pos", level );
  if( x.geom_tan_pos.nelem() )
    Print( x.geom_tan_pos, "geom_tan_pos", level );
}


void Print(
        // WS Generic Input:
        const Timer&   /* x */,
        // WS Generic Input Names:
        const String&  /* x_name */,
        // Keywords:
        const Index    /* level */ )
{
/*  ostringstream os;
  cout << "  *" << x_name <<"* =";
  SWITCH_OUTPUT (level, "  *" << x_name << "*:\n")
  for( Index i=0; i<x.nelem(); i++ )
    os << x[i] << '\n';
  SWITCH_OUTPUT (level, os.str ()) */
}


void
timerStart (// WS Output
            Timer& starttime)
{
  if ((starttime.realtime = times (&starttime.cputime)) == -1)
    throw runtime_error ("Timer error: Unable to get current CPU time");
}


void
timerStop (// WS Input
           const Timer& starttime)
{
  Timer endtime;
  static long clktck = 0;

  if (clktck == 0)
    if ((clktck = sysconf (_SC_CLK_TCK)) < 0)
      throw runtime_error ("Timer error: Unable to determine CPU clock ticks");

  if ((endtime.realtime = times (&endtime.cputime)) == -1)
    throw runtime_error ("Timer error: Unable to get current CPU time");

  // FIXME: out1 does not have setf, but we need to set it here
  cout.setf (ios::showpoint | ios::fixed);

  out1 << "  * CPU time  total: " << setprecision (2)
    << ((endtime.cputime.tms_stime - starttime.cputime.tms_stime)
        + (endtime.cputime.tms_utime - starttime.cputime.tms_utime))
    / (Numeric)clktck;

  out1 << "  user: " << setprecision (2)
    << (endtime.cputime.tms_utime - starttime.cputime.tms_utime)
    / (Numeric)clktck;

  out1 << "  system: " << setprecision (2)
    << (endtime.cputime.tms_stime - starttime.cputime.tms_stime)
    / (Numeric)clktck;

  out1 << "\n               real: " << setprecision (2)
    << (endtime.realtime - starttime.realtime) / (Numeric)clktck;

  out1 << "  " << setprecision (2)
    << ((endtime.cputime.tms_stime - starttime.cputime.tms_stime)
        + (endtime.cputime.tms_utime - starttime.cputime.tms_utime))
    / (Numeric)(endtime.realtime - starttime.realtime) * 100.
    << "%CPU\n";
}


//! Error
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void Error(
           const String& msg )
{
  out0 << msg << "\n";
  arts_exit();
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



