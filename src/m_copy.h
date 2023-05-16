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

#include "agenda_class.h"
#include "debug.h"
#include "messages.h"
#include "mystring.h"
#include "workspace_ng.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <class T>
void Copy(  // WS Generic Output:
    T& out,
    const String& /* out_name */,
    // WS Generic Input:
    const T& in,
    const String& /* in_name */,
    const Verbosity&) {
  out = in;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Copy(Workspace& ws,
          // WS Generic Output:
          Agenda& out,
          const String& out_name,
          // WS Generic Input:
          const Agenda& in,
          const String& /* in_name */,
          const Verbosity& verbosity) {
  out = in;
  out.set_name(out_name);
  out.check(ws, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Copy(Workspace& ws,
          // WS Generic Output:
          ArrayOfAgenda& out,
          const String& out_name,
          // WS Generic Input:
          const ArrayOfAgenda& in,
          const String& /* in_name */,
          const Verbosity& verbosity) {
  out = in;
  for (auto & it : out) {
    it.set_name(out_name);
    it.check(ws, verbosity);
  }
}

#endif  // m_copy_h
