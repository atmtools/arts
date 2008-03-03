/* Copyright (C) 2003-2007 Oliver Lemke <olemke@core-dump.info>
  
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

#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Please run ./configure in the top arts directory before compiling."
#endif

#include <iostream>
#if HAVE_UNISTD_H
# include <sys/types.h>
# include <unistd.h>
#endif
#ifdef _POSIX_VERSION
#include <sys/times.h>
#endif
#include <stdexcept>
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "mystring.h"

class Workspace;

#define SWITCH_OUTPUT(x,y) switch (x) { \
case 0: out0 << y;break; \
case 1: out1 << y;break; \
case 2: out2 << y;break; \
case 3: out3 << y;break; \
default: throw runtime_error ("Output level must have value from 0-3"); \
}


typedef struct {
#ifdef _POSIX_VERSION
  struct tms cputime;
  clock_t realtime;
#endif
} Timer;

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
Print(
      // WS Generic Input:
      const T&        x,
      // Keywords:
      const Index&   level )
{
  SWITCH_OUTPUT (level, x);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
Print(
      // WS Generic Input:
      const ArrayOfGridPos&   x,
      // Keywords:
      const Index&             level );


/* Workspace method: Doxygen documentation will be auto-generated */
void
Print(
        // WS Generic Input:
        const ArrayOfIndex&   x,
        // Keywords:
        const Index&          level );


/* Workspace method: Doxygen documentation will be auto-generated */
void
Print(
        // WS Generic Input:
        const ArrayOfString&   x,
        // Keywords:
        const Index&           level );


/* Workspace method: Doxygen documentation will be auto-generated */
void
Print(
        // WS Generic Input:
        const Ppath&    ppath,
        // Keywords:
        const Index&    level );


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const ArrayOfPpath&   x,
        // Keywords:
        const Index&            level );


/* Workspace method: Doxygen documentation will be auto-generated */
void Print(
        // WS Generic Input:
        const Timer&   x,
        // Keywords:
        const Index&   level);


/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(
        // Workspace reference
        Workspace& ws,
        // Keywords:
        const Index&   only_allocated,
        const Index&   level);

#endif /* m_general_h */

