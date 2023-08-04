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

#include <workspace.h>
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
