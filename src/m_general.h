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

#include <iostream>
#include <stdexcept>

#include "agenda_class.h"
#include "arts.h"
#include "cia.h"
#include "mystring.h"
#include "ppath_struct.h"
#include "special_interp.h"
#include "tessem.h"
#include "timer_struct.h"

class Workspace;

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void Print(
    // WS Generic Input:
    const T& x,
    // Keywords:
    const Index& level) {
  if (level) std::cerr << x << '\n';
  else std::cout << x << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws,
           // WS Generic Input:
           const Agenda& x,
           // Keywords:
           const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(Workspace& ws,
           // WS Generic Input:
           const ArrayOfAgenda& x,
           // Keywords:
           const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfGridPos& x,
    // Keywords:
    const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfCIARecord& x,
    // Keywords:
    const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfString& x,
    // Keywords:
    const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const ArrayOfPpath& x,
    // Keywords:
    const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const Timer& x,
    // Keywords:
    const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void Print(  // WS Generic Input:
    const TessemNN& x,
    // Keywords:
    const Index& level);

/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(  // Workspace reference
    Workspace& ws,
    // Keywords:
    const Index& only_allocated,
    const Index& level);

#endif /* m_general_h */
