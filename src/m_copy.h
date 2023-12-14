/*!
  \file   m_copy.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Copy.
  
  This file contains the implementation of the supergeneric method
  Copy.
*/

#ifndef m_copy_h
#define m_copy_h

#include <workspace.h>

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Copy(  // WS Generic Output:
    T& out,
    // WS Generic Input:
    const T& in) {
  out = in;
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Set(  // WS Generic Output:
    T& out,
    const T& in) {
  out = in;
}

#endif  // m_copy_h
