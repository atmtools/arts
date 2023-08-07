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

#include "debug.h"
#include "mystring.h"
#include <workspace.h>

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Copy(  // WS Generic Output:
    T& out,
    const String& /* out_name */,
    // WS Generic Input:
    const T& in,
    const String& /* in_name */) {
  out = in;
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <>
inline void Copy<Agenda>(// WS Generic Output:
                         Agenda& out,
                         const String& out_name,
                         // WS Generic Input:
                         const Agenda& in,
                         const String& /* in_name */) {
  out = in;
  out.set_name(out_name);
  out.finalize();
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <>
inline void Copy<ArrayOfAgenda>(// WS Generic Output:
                                ArrayOfAgenda& out,
                                const String& out_name,
                                // WS Generic Input:
                                const ArrayOfAgenda& in,
                                const String& /* in_name */) {
  out = in;
  for (auto& it : out) {
    it.set_name(out_name);
    it.finalize();
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Set(  // WS Generic Output:
    T& out,
    const T& in) {
  out = in;
}

#endif  // m_copy_h
