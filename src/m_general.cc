/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions of a general and overall character.

  This file is for general functions that do not fit in any other "m_"-file.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"

#include <stdexcept>
#include <cstdlib>
#ifdef TIME_SUPPORT
#include <unistd.h>
#endif

#include "m_general.h"
#include "array.h"
#include "check_input.h"
#include "messages.h"
#include "mystring.h"

#include "math_funcs.h"
#include "wsv_aux.h"

#include "auto_md.h"
#include "workspace_ng.h"

#include "sensor.h"

#include "fastem.h"
#include "tessem.h"

extern const Numeric SPEED_OF_LIGHT;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void INCLUDE(const Verbosity&)
{
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws _U_,
           // WS Generic Input:
           const Agenda& x,
           // Keywords:
           const Index& level,
           const Verbosity& verbosity)
{
  ostringstream os;
  os << "    " << x.name() << " {\n";
  x.print(os, "        ");
  os << "    " << "}";
  CREATE_OUTS;
  SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws _U_,
           // WS Generic Input:
           const ArrayOfAgenda&   x,
           // Keywords:
           const Index&    level,
           const Verbosity& verbosity)
{
  ostringstream os;
  os << "    " << x.nelem() << " agendas: {\n";
  for (Index i = 0; i < x.nelem(); i++)
  {
      os << "      " << x[i].name() << ": {\n";
      x[i].print(os, "          ");
      os << "      " << "}\n";
  }
  os << "    " << "}";
  CREATE_OUTS;
  SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(// WS Generic Input:
           const ArrayOfGridPos& x,
           // Keywords:
           const Index&          level,
           const Verbosity&      verbosity)
{
  ostringstream os;
  for( Index i=0; i<x.nelem(); i++ )
  {
      if (i) os << '\n';
      os << "  " << x[i].idx << "  " << x[i].fd[0] << "  " << x[i].fd[1];
  }
  CREATE_OUTS;
  SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(// WS Generic Input:
           const ArrayOfCIARecord& cia_data,
           // Keywords:
           const Index&         level,
           const Verbosity&     verbosity)
{
    CREATE_OUTS;

    ostringstream os;
    os << "  CIA tag; Spectral range [cm-1]; Temp range [K]; # of sets\n";
    for (Index i = 0; i < cia_data.nelem(); i++)
        for (Index j = 0; j < cia_data[i].DatasetCount(); j++)
        {
            Vector temp_grid = cia_data[i].TemperatureGrid(j);
            Vector freq_grid = cia_data[i].FrequencyGrid(j);

            os << setprecision(2) << std::fixed << "  "
            << cia_data[i].MoleculeName(0) << "-CIA-" << cia_data[i].MoleculeName(1)
            << "-" << j
            << "; " << freq_grid[0] / 100. / SPEED_OF_LIGHT
            << " - " << freq_grid[freq_grid.nelem()-1] / 100. / SPEED_OF_LIGHT
            << std::fixed
            << "; " << temp_grid[0] << " - " << temp_grid[temp_grid.nelem()-1]
            << "; " << temp_grid.nelem()
            << "\n";
        }
    SWITCH_OUTPUT(level, os.str());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(// WS Generic Input:
           const ArrayOfString& x,
           // Keywords:
           const Index&         level,
           const Verbosity&     verbosity)
{
  ostringstream os;
  for( Index i=0; i<x.nelem(); i++ )
  {
      if (i) os << '\n';
      os << "  " << x[i];
  }
  CREATE_OUTS;
  SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
Print(// WS Generic Input:
      const Ppath&     x,
      // Keywords:
      const Index&     level,
      const Verbosity& verbosity)
{
  CREATE_OUTS;
  SWITCH_OUTPUT (level, "dim: ");
  Print( x.dim, level, verbosity );
  SWITCH_OUTPUT (level, "np: ");
  Print( x.np, level, verbosity );
  SWITCH_OUTPUT (level, "constant: ");
  Print( x.constant, level, verbosity );
  SWITCH_OUTPUT (level, "background: ");
  Print( x.background, level, verbosity );
  SWITCH_OUTPUT (level, "start_pos: ");
  Print( x.start_pos, level, verbosity );
  SWITCH_OUTPUT (level, "start_los: ");
  Print( x.start_los, level, verbosity );
  SWITCH_OUTPUT (level, "start_lstep: ");
  Print( x.start_lstep, level, verbosity );
  SWITCH_OUTPUT (level, "pos: ");
  Print( x.pos, level, verbosity );
  SWITCH_OUTPUT (level, "los: ");
  Print( x.los, level, verbosity );
  SWITCH_OUTPUT (level, "r: ");
  Print( x.r, level, verbosity );
  SWITCH_OUTPUT (level, "lstep: ");
  Print( x.lstep, level, verbosity );
  SWITCH_OUTPUT (level, "end_pos: ");
  Print( x.end_pos, level, verbosity );
  SWITCH_OUTPUT (level, "end_los: ");
  Print( x.end_los, level, verbosity );
  SWITCH_OUTPUT (level, "end_lstep: ");
  Print( x.end_lstep, level, verbosity );
  SWITCH_OUTPUT (level, "nreal: ");
  Print( x.nreal, level, verbosity );
  SWITCH_OUTPUT (level, "ngroup: ");
  Print( x.ngroup, level, verbosity );
  SWITCH_OUTPUT (level, "gp_p: ");
  Print( x.gp_p, level, verbosity );
  if( x.dim >= 2 )
    {
      SWITCH_OUTPUT (level, "gp_lat: ");
      Print( x.gp_lat, level, verbosity );
    }
  if( x.dim == 3 )
    {
      SWITCH_OUTPUT (level, "gp_lon: ");
      Print( x.gp_lon, level, verbosity );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(// WS Generic Input:
           const ArrayOfPpath& x,
           // Keywords:
           const Index&        level,
           const Verbosity&    verbosity)
{
  CREATE_OUTS;
  for( Index i=0; i<x.nelem(); i++ )
    {
      ostringstream os;
      os << "Ppath element " << i << ": ";
      SWITCH_OUTPUT (level, os.str ());
      Print( x[i], level, verbosity );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef TIME_SUPPORT
void Print(// WS Generic Input:
           const Timer& timer,
           // Keywords:
           const Index& level,
           const Verbosity& verbosity)
{
    CREATE_OUTS;

    if (!timer.finished)
    {
        SWITCH_OUTPUT(level, "Timer error: Nothing to output. Use timerStart/timerStop first.");
        return;
    }

    ostringstream os;
    os.setf (ios::showpoint | ios::fixed);

    static long clktck = 0;

    if (clktck == 0)
        if ((clktck = sysconf (_SC_CLK_TCK)) < 0)
            throw runtime_error ("Timer error: Unable to determine CPU clock ticks");

    os << "  * CPU time  total: " << setprecision (2)
    << (Numeric)((timer.cputime_end.tms_stime - timer.cputime_start.tms_stime)
                 + (timer.cputime_end.tms_utime - timer.cputime_start.tms_utime))
    / (Numeric)clktck;

    os << "  user: " << setprecision (2)
    << (Numeric)(timer.cputime_end.tms_utime - timer.cputime_start.tms_utime)
    / (Numeric)clktck;

    os << "  system: " << setprecision (2)
    << (Numeric)(timer.cputime_end.tms_stime - timer.cputime_start.tms_stime)
    / (Numeric)clktck;

    os << "\n               real: " << setprecision (2)
    << (Numeric)(timer.realtime_end - timer.realtime_start) / (Numeric)clktck;

    os << "  " << setprecision (2)
    << (Numeric)((timer.cputime_end.tms_stime - timer.cputime_start.tms_stime)
                 + (timer.cputime_end.tms_utime - timer.cputime_start.tms_utime))
    / (Numeric)(timer.realtime_end - timer.realtime_start) * 100.
    << "%CPU\n";

    SWITCH_OUTPUT(level, os.str());
}
#else
void Print(// WS Generic Input:
           const Timer&,
           // Keywords:
           const Index& level,
           const Verbosity& verbosity)
{
    SWITCH_OUTPUT (level, "Timer error: ARTS was compiled without timer support");
}
#endif


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(// WS Generic Input:
        const TessemNN& x,
        // Keywords:
        const Index& level,
        const Verbosity& verbosity)
{
    CREATE_OUTS;
    ostringstream os;
    os << "TessemNN size: Inputs = " << x.nb_inputs << ", Outputs = " << x.nb_outputs << ", Cache = " << x.nb_cache;
    SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(// Workspace reference
                    Workspace& ws,
                    // Keywords:
                    const Index&     only_allocated,
                    const Index&     level,
                    const Verbosity& verbosity)
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
          if (ws.is_initialized(i)) os << ws.depth(i);
          os << "\n";
        }
      else if (ws.is_initialized(i))
        {
          os << "    ";
          PrintWsvName (os, i);
          os << ws.depth(i) << "\n";
        }
    }
  CREATE_OUTS;
  SWITCH_OUTPUT (level, os.str ());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
StringJoin(String& out,
              const String& in1,
              const String& in2,
              const String& in3,
              const String& in4,
              const String& in5,
              const String& in6,
              const String& in7,
              const String& in8,
              const String& in9,
              const String& in10,
              const Verbosity&
              )
{
    out = in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10;
}


/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef TIME_SUPPORT
void
timerStart (// WS Output
            Timer& timer,
            const Verbosity&)
{
    if ((timer.realtime_start = times (&timer.cputime_start)) == (clock_t)-1)
        throw runtime_error ("Timer error: Unable to get current CPU time");

    timer.running = true;
    timer.finished = false;
}
#else
void
timerStart (// WS Output
            Timer& /*starttime*/,
            const Verbosity&)
{
  throw runtime_error ("Timer error: ARTS was compiled without POSIX support, thus timer\nfunctions are not available.");
}
#endif


/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef TIME_SUPPORT
void
timerStop (// WS Input
           Timer& timer,
           const Verbosity&)
{
    if (!timer.running)
        throw runtime_error("Timer error: Unable to stop timer that's not running.");

    if ((timer.realtime_end = times (&(timer.cputime_end))) == (clock_t)-1)
        throw runtime_error ("Timer error: Unable to get current CPU time");

    timer.running = false;
    timer.finished = true;
}
#else
void
timerStop (// WS Input
           const Timer&,
           const Verbosity&)
{
  throw runtime_error ("Timer error: ARTS was compiled without POSIX support, thus timer\nfunctions are not available.");
}
#endif

/* Workspace method: Doxygen documentation will be auto-generated */
void Error(const String& msg, const Verbosity& verbosity)
{
  CREATE_OUT0;
  throw runtime_error(msg);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void Exit(const Verbosity& verbosity)
{
  CREATE_OUT1;
  out1 << "  Forced exit.\n";
  arts_exit (EXIT_SUCCESS);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void TestArrayOfAgenda(Workspace& ws,
                       const ArrayOfAgenda& test_agenda_array,
                       const Index& index,
                       const Verbosity&)
{
    ostringstream os;
    os << "  Local value of iy_unit, agenda #"
       << index << " of " << test_agenda_array.nelem();
    test_agenda_arrayExecute(ws, index, os.str(), test_agenda_array);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Test(const ArrayOfArrayOfLineRecord& ast, const ArrayOfQuantumIdentifier& aqi, const Verbosity& verbosity)
{
    CREATE_OUT1;

    ArrayOfIndex matches;
    ArrayOfQuantumMatchInfo match_info;

    for (Index qi = 0; qi < aqi.nelem(); qi++)
    {
        out1 << "QI: " << aqi[qi] << "\n";
        for (Index s = 0; s < ast.nelem(); s++)
        {
            const ArrayOfLineRecord& species_lines = ast[s];

            match_lines_by_quantum_identifier(matches, match_info, aqi[qi], species_lines);

            // Output info about matched lines for debugging purposes
            for (Index i = 0; i < matches.nelem(); i++)
            {
                const LineRecord& lr = ast[s][matches[i]];
                out1 << "  Line: " << lr.Name() << " ";
                switch (match_info[i].Upper())
                {
                    case QMI_NONE:    out1 << "Upper None    "; break;
                    case QMI_FULL:    out1 << "Upper full    "; break;
                    case QMI_PARTIAL: out1 << "Upper partial "; break;
                }
                switch (match_info[i].Lower())
                {
                    case QMI_NONE:    out1 << "Lower None    "; break;
                    case QMI_FULL:    out1 << "Lower full    "; break;
                    case QMI_PARTIAL: out1 << "Lower partial "; break;
                }
                out1 << " abs_lines_per_species[" << s << "][" << matches[i] << "]\n";
//                out1 << "  UP: " << lr.QuantumNumbers().Upper()
//                << " LO: " << lr.QuantumNumbers().Lower() << "\n";
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void verbosityInit(// WS Output:
                   Verbosity& verbosity)
{
  extern Verbosity verbosity_at_launch;
  
  verbosity.set_screen_verbosity(verbosity_at_launch.get_screen_verbosity());
  verbosity.set_agenda_verbosity(verbosity_at_launch.get_agenda_verbosity());
  verbosity.set_file_verbosity(verbosity_at_launch.get_file_verbosity());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void verbositySet(// WS Output:
                  Verbosity& verbosity,
                  // WS Generic Input:
                  const Index& agenda,
                  const Index& screen,
                  const Index& file)
{
  verbosity.set_agenda_verbosity(agenda);
  verbosity.set_screen_verbosity(screen);
  verbosity.set_file_verbosity(file);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void verbositySetAgenda(// WS Output:
                        Verbosity& verbosity,
                        // WS Generic Input:
                        const Index& level)
{
  verbosity.set_agenda_verbosity(level);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void verbositySetFile(// WS Output:
                      Verbosity& verbosity,
                      // WS Generic Input:
                      const Index& level)
{
  verbosity.set_file_verbosity(level);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void verbositySetScreen(// WS Output:
                        Verbosity& verbosity,
                        // WS Generic Input:
                        const Index& level)
{
  verbosity.set_screen_verbosity(level);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GetEnvironmentVariable(// WS Generic Output:
                            String& str,
                            // WS Generic Input:
                            const String& envvar,
                            const Verbosity& /* verbosity */)
{
    char *cstr;
    cstr = std::getenv (envvar.c_str());
    if (cstr == NULL)
    {
        std::ostringstream os;
        os << "Environment variable " << envvar << " does not exist.";
        throw std::runtime_error(os.str());
    }
    str = cstr != NULL?String(cstr):"";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetEnvironmentVariable(// WS Generic Output:
                            Index& i,
                            // WS Generic Input:
                            const String& envvar,
                            const Verbosity& /* verbosity */)
{
    char *cstr;
    cstr = std::getenv (envvar.c_str());
    if (cstr == NULL || std::strlen(cstr) == 0)
    {
        std::ostringstream os;
        os << "Environment variable " << envvar << " is empty or does not exist.";
        throw std::runtime_error(os.str());
    }
    std::istringstream is(cstr);
    is >> i;
    if (!is.eof())
    {
        ostringstream os;
        os << "Cannot convert environment variable " << envvar << " to Index: " << cstr;
        throw std::runtime_error(os.str());
    }
}
