/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   m_general.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-07-24
 
  \brief  Template functions for general supergeneric ws methods.
 
*/

#ifndef m_general_h
#define m_general_h

#include "arts.h"

#include <iostream>
#ifdef TIME_SUPPORT
#include <sys/times.h>
#endif
#include <stdexcept>
#include "cia.h"
#include "messages.h"
#include "mystring.h"
#include "ppath.h"
#include "special_interp.h"
#include "tessem.h"

class Workspace;

#define SWITCH_OUTPUT(x, y)                                           \
  {                                                                   \
    ostringstream so_os;                                              \
    so_os << y << '\n';                                               \
    switch (x) {                                                      \
      case 0:                                                         \
        out0 << so_os.str();                                          \
        break;                                                        \
      case 1:                                                         \
        out1 << so_os.str();                                          \
        break;                                                        \
      case 2:                                                         \
        out2 << so_os.str();                                          \
        break;                                                        \
      case 3:                                                         \
        out3 << so_os.str();                                          \
        break;                                                        \
      default:                                                        \
        throw runtime_error("Output level must have value from 0-3"); \
    }                                                                 \
  }

class Timer {
 public:
  Timer() : running(false), finished(false) {}
  bool running;
  bool finished;
#ifdef TIME_SUPPORT
  struct tms cputime_start;
  clock_t realtime_start;
  struct tms cputime_end;
  clock_t realtime_end;
#endif
};

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void Print(
    // WS Generic Input:
    const T& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity) {
  CREATE_OUTS;
  SWITCH_OUTPUT(level, x)
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws,
           // WS Generic Input:
           const Agenda& x,
           // Keywords:
           const Index& level,
           const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws,
           // WS Generic Input:
           const ArrayOfAgenda& x,
           // Keywords:
           const Index& level,
           const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfGridPos& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfCIARecord& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfString& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const Ppath& ppath,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfPpath& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const Timer& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const TessemNN& x,
    // Keywords:
    const Index& level,
    const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(  // Workspace reference
    Workspace& ws,
    // Keywords:
    const Index& only_allocated,
    const Index& level,
    const Verbosity& verbosity);

#endif /* m_general_h */
