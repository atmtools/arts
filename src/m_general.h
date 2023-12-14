/*!
  \file   m_general.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-07-24
 
  \brief  Template functions for general supergeneric ws methods.
 
*/

#ifndef m_general_h
#define m_general_h

#include <workspace.h>

#include <iostream>

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Print(
    // WS Generic Input:
    const T& x,
    // Keywords:
    const Index& level) {
  if (level) std::cerr << x << '\n';
  else std::cout << x << '\n';
}

#endif /* m_general_h */
