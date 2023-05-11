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

#include <chrono>
#include <cstdlib>
#include <ratio>
#include <stdexcept>

#include "array.h"
#include "arts_constants.h"
#include "check_input.h"
#include "m_general.h"
#include "mystring.h"

#include "math_funcs.h"
#include "wsv_aux.h"

#include "auto_md.h"
#include "workspace_ng.h"

#include "sensor.h"

#include "fastem.h"
#include "tessem.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void INCLUDE() {}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws _U_,
           // WS Generic Input:
           const Agenda& x,
           // Keywords:
           const Index& level) {
  ostringstream os;
  os << "    " << x.name() << " {\n";
  x.print(os, "        ");
  os << "    "
     << "}";
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws _U_,
           // WS Generic Input:
           const ArrayOfAgenda& x,
           // Keywords:
           const Index& level) {
  ostringstream os;
  os << "    " << x.nelem() << " agendas: {\n";
  for (Index i = 0; i < x.nelem(); i++) {
    os << "      " << x[i].name() << ": {\n";
    x[i].print(os, "          ");
    os << "      "
       << "}\n";
  }
  os << "    "
     << "}";
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfGridPos& x,
    // Keywords:
    const Index& level) {
  ostringstream os;
  for (Index i = 0; i < x.nelem(); i++) {
    if (i) os << '\n';
    os << "  " << x[i].idx << "  " << x[i].fd[0] << "  " << x[i].fd[1];
  }
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfCIARecord& cia_data,
    // Keywords:
    const Index& level) {
  ostringstream os;
  os << "  CIA tag; Spectral range [cm-1]; Temp range [K]; # of sets\n";
  for (Index i = 0; i < cia_data.nelem(); i++)
    for (Index j = 0; j < cia_data[i].DatasetCount(); j++) {
      Vector temp_grid{cia_data[i].TemperatureGrid(j)};
      Vector freq_grid{cia_data[i].FrequencyGrid(j)};

      os << setprecision(2) << std::fixed << "  " << cia_data[i].MoleculeName(0)
         << "-CIA-" << cia_data[i].MoleculeName(1) << "-" << j << "; "
         << freq_grid[0] / 100. / SPEED_OF_LIGHT << " - "
         << freq_grid[freq_grid.nelem() - 1] / 100. / SPEED_OF_LIGHT
         << std::fixed << "; " << temp_grid[0] << " - "
         << temp_grid[temp_grid.nelem() - 1] << "; " << temp_grid.nelem()
         << "\n";
    }
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfString& x,
    // Keywords:
    const Index& level) {
  ostringstream os;
  for (Index i = 0; i < x.nelem(); i++) {
    if (i) os << '\n';
    os << "  " << x[i];
  }
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfPpath& x,
    // Keywords:
    const Index& level) {
  for (Index i = 0; i < x.nelem(); i++) {
    ostringstream os;
    os << "Ppath element " << i << ": ";
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const Timer& timer,
    // Keywords:
    const Index& level) {
  if (!timer.finished) {
  if (level) std::cerr << "Timer error: Nothing to output. Use timerStart/timerStop first." << '\n';
  else std::cout << "Timer error: Nothing to output. Use timerStart/timerStop first." << '\n';
    return;
  }

  ostringstream os;

  const auto cputime =
      static_cast<double>(timer.cputime_end - timer.cputime_start) /
      CLOCKS_PER_SEC;
  const auto walltime = std::chrono::duration<double, std::ratio<1>>(
                            timer.realtime_end - timer.realtime_start)
                            .count();
  os << std::fixed << setprecision(2) << "  * Timing: CPU " << cputime << "s, "
     << "Wall " << walltime << "s, " << 100. * cputime / walltime << "%CPU\n";

  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const TessemNN& x,
    // Keywords:
    const Index& level) {
  ostringstream os;
  os << "TessemNN size: Inputs = " << x.nb_inputs
     << ", Outputs = " << x.nb_outputs << ", Cache = " << x.nb_cache;
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(  // Workspace reference
    Workspace& ws,
    // Keywords:
    const Index& only_allocated,
    const Index& level) {
  ostringstream os;

  if (only_allocated)
    os << "  Allocated workspace variables: \n";
  else
    os << "  Workspace variables: \n";
  for (Index i = 0; i < ws.nelem(); i++) {
    if (!only_allocated) {
      os << "    ";
      ws.PrintWsvName(os, i);
      if (ws.is_initialized(i)) os << ws.depth(i);
      os << "\n";
    } else if (ws.is_initialized(i)) {
      os << "    ";
      ws.PrintWsvName(os, i);
      os << ws.depth(i) << "\n";
    }
  }
  if (level) std::cerr << os.str() << '\n';
  else std::cout << os.str() << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void StringJoin(String& out,
                const String& in1,
                const String& in2,
                const String& in3,
                const String& in4,
                const String& in5,
                const String& in6,
                const String& in7,
                const String& in8,
                const String& in9,
                const String& in10) {
  out = in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void timerStart(  // WS Output
    Timer& timer) {
  timer.cputime_start = std::clock();
  timer.realtime_start = std::chrono::high_resolution_clock::now();

  timer.running = true;
  timer.finished = false;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void timerStop(  // WS Input
    Timer& timer) {
  ARTS_USER_ERROR_IF(!timer.running,
                     "Timer error: Unable to stop timer that's not running.");

  timer.realtime_end = std::chrono::high_resolution_clock::now();
  timer.cputime_end = std::clock();

  timer.running = false;
  timer.finished = true;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Error(const String& msg) {
  ARTS_USER_ERROR ( msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Exit() {
  arts_exit(EXIT_SUCCESS);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void TestArrayOfAgenda(Workspace& ws,
                       const ArrayOfAgenda& test_agenda_array,
                       const Index& index) {
  ostringstream os;
  os << "  Local value of iy_unit, agenda #" << index << " of "
     << test_agenda_array.nelem();
  test_agenda_arrayExecute(ws, index, os.str(), test_agenda_array);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetEnvironmentVariable(  // WS Generic Output:
    String& str,
    // WS Generic Input:
    const String& envvar /* verbosity */) {
  char* cstr;
  cstr = std::getenv(envvar.c_str());
  ARTS_USER_ERROR_IF (cstr == NULL,
    "Environment variable ", envvar, " does not exist.")
  str = cstr != NULL ? String(cstr) : "";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetNumberOfThreads(Index& nthreads) {
  nthreads = arts_omp_get_max_threads();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetEnvironmentVariable(  // WS Generic Output:
    Index& i,
    // WS Generic Input:
    const String& envvar /* verbosity */) {
  char* cstr;
  cstr = std::getenv(envvar.c_str());
  ARTS_USER_ERROR_IF (cstr == NULL || std::strlen(cstr) == 0,
    "Environment variable ", envvar, " "
    "is empty or does not exist.")
  std::istringstream is(cstr);
  is >> i;
  ARTS_USER_ERROR_IF (!is.eof(),
      "Cannot convert environment variable ", envvar, " "
      "to Index: ", cstr)
}

/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef _OPENMP
void SetNumberOfThreads(const Index& nthreads) {
  omp_set_num_threads((int)nthreads);
}
#else
void SetNumberOfThreads(const Index& /* nthreads */,
                    ) {
  CREATE_OUT1;
  out1 << "No OpenMP support. Can't change number of threads.\n";
}
#endif
