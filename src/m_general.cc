/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
#ifdef TIME_SUPPORT
#include <unistd.h>
#endif

#include "arts.h"

#include "m_general.h"
#include "array.h"
#include "check_input.h"
#include "messages.h"
#include "mystring.h"

#include "math_funcs.h"
#include "make_vector.h"
#include "wsv_aux.h"

#include "workspace_ng.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void INCLUDE()
{
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        Workspace& ws _U_,
        // WS Generic Input:
        const Agenda&   x,
        // Keywords:
        const Index&            level )
{
  ostringstream os;
  os << "     " << x << "\n";
  SWITCH_OUTPUT (level, os.str ())
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const ArrayOfGridPos&   x,
        // Keywords:
        const Index&            level )
{
  ostringstream os;
  for( Index i=0; i<x.nelem(); i++ )
    os << "     " << x[i].idx << "  " << x[i].fd[0] << "  " << x[i].fd[1]
         << "\n";
  SWITCH_OUTPUT (level, os.str ())
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const ArrayOfIndex&   x,
        // Keywords:
        const Index&          level )
{
  ostringstream os;
  for( Index i=0; i<x.nelem(); i++ )
    os << x[i] << " ";
  SWITCH_OUTPUT (level, os.str () << '\n')
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const ArrayOfString&   x,
        // Keywords:
        const Index&           level )
{
  ostringstream os;
  for( Index i=0; i<x.nelem(); i++ )
    os << x[i] << '\n';
  SWITCH_OUTPUT (level, os.str ())
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
Print(
      // WS Generic Input:
      const Ppath&    x,
      // Keywords:
      const Index&    level )
{
  SWITCH_OUTPUT (level, "dim: ")
  Print( x.dim, level );
  SWITCH_OUTPUT (level, "np: ")
  Print( x.np, level );
  SWITCH_OUTPUT (level, "refraction: ")
  Print( x.refraction, level );
  SWITCH_OUTPUT (level, "method: ")
  Print( x.method, level );
  SWITCH_OUTPUT (level, "constant: ")
  Print( x.constant, level );
  SWITCH_OUTPUT (level, "pos: ")
  Print( x.pos, level );
  SWITCH_OUTPUT (level, "z: ")
  Print( x.z, level );
  SWITCH_OUTPUT (level, "l_step: ")
  Print( x.l_step, level );
  SWITCH_OUTPUT (level, "gp_p: ")
  Print( x.gp_p, level );
  if( x.dim >= 2 )
    {
      SWITCH_OUTPUT (level, "gp_lat: ")
      Print( x.gp_lat, level );
    }
  if( x.dim == 3 )
    {
      SWITCH_OUTPUT (level, "gp_lon: ")
      Print( x.gp_lon, level );
    }
  SWITCH_OUTPUT (level, "los: ")
  Print( x.los, level );
  SWITCH_OUTPUT (level, "background: ")
  Print( x.background, level );
  if( x.tan_pos.nelem() )
    {
      SWITCH_OUTPUT (level, "tan_pos: ")
      Print( x.tan_pos, level );
    }
  if( x.geom_tan_pos.nelem() )
    {
      SWITCH_OUTPUT (level, "geom_tan_pos: ")
      Print( x.geom_tan_pos, level );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const ArrayOfPpath&   x,
        // Keywords:
        const Index&            level )
{
  for( Index i=0; i<x.nelem(); i++ )
    {
      ostringstream os;
      os << "Ppath element " << i << ": ";
      SWITCH_OUTPUT (level, os.str ())
      Print( x[i], level );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const Timer&   /* x */,
        // Keywords:
        const Index&   /* level */ )
{
/*  ostringstream os;
  cout << "  *" << x_name <<"* =";
  SWITCH_OUTPUT (level, "  *" << x_name << "*:\n")
  for( Index i=0; i<x.nelem(); i++ )
    os << x[i] << '\n';
  SWITCH_OUTPUT (level, os.str ()) */
}


/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(
        // Workspace reference
        Workspace& ws,
        // Keywords:
        const Index&   only_allocated,
        const Index&   level)
{
  ostringstream os;

  if (only_allocated)
    os << "  Allocated workspace variables: \n";
  else
    os << "  Workspace variables: \n";
  for (Index i = 0; i < ws.nelem(); i++)
    {
      if (!only_allocated)
        {
          os << "    ";
          PrintWsvName (os, i);
          if (ws.is_initialized(i)) os << "+";
          os << "\n";
        }
      else if (ws.is_initialized(i))
        {
          os << "    ";
          PrintWsvName (os, i);
          os << "\n";
        }
    }
  SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef TIME_SUPPORT
void
timerStart (// WS Output
            Timer& starttime)
{
  if ((starttime.realtime = times (&starttime.cputime)) == (clock_t)-1)
    throw runtime_error ("Timer error: Unable to get current CPU time");
}
#else
void
timerStart (// WS Output
            Timer& /*starttime*/)
{
  throw runtime_error ("Timer error: ARTS was compiled without POSIX support, thus timer\nfunctions are not available.");
}
#endif


/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef TIME_SUPPORT
void
timerStop (// WS Input
           const Timer& starttime)
{
  Timer endtime;
  static long clktck = 0;

  if (clktck == 0)
    if ((clktck = sysconf (_SC_CLK_TCK)) < 0)
      throw runtime_error ("Timer error: Unable to determine CPU clock ticks");

  if ((endtime.realtime = times (&endtime.cputime)) == (clock_t)-1)
    throw runtime_error ("Timer error: Unable to get current CPU time");

  // FIXME: out1 does not have setf, but we need to set it here
  cout.setf (ios::showpoint | ios::fixed);

  out1 << "  * CPU time  total: " << setprecision (2)
    << (Numeric)((endtime.cputime.tms_stime - starttime.cputime.tms_stime)
        + (endtime.cputime.tms_utime - starttime.cputime.tms_utime))
    / (Numeric)clktck;

  out1 << "  user: " << setprecision (2)
    << (Numeric)(endtime.cputime.tms_utime - starttime.cputime.tms_utime)
    / (Numeric)clktck;

  out1 << "  system: " << setprecision (2)
    << (Numeric)(endtime.cputime.tms_stime - starttime.cputime.tms_stime)
    / (Numeric)clktck;

  out1 << "\n               real: " << setprecision (2)
    << (Numeric)(endtime.realtime - starttime.realtime) / (Numeric)clktck;

  out1 << "  " << setprecision (2)
    << (Numeric)((endtime.cputime.tms_stime - starttime.cputime.tms_stime)
        + (endtime.cputime.tms_utime - starttime.cputime.tms_utime))
    / (Numeric)(endtime.realtime - starttime.realtime) * 100.
    << "%CPU\n";
}
#else
void
timerStop (// WS Input
           const Timer&)
{
  throw runtime_error ("Timer error: ARTS was compiled without POSIX support, thus timer\nfunctions are not available.");
}
#endif

/* Workspace method: Doxygen documentation will be auto-generated */
void Error(
           const String& msg )
{
  out0 << msg << "\n";
  arts_exit();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Exit()
{
  out1 << "  Forced exit.\n";
  arts_exit (EXIT_SUCCESS);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Test( )
{
  // This function can be used to test stuff.

}



