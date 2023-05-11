/*!
  \file   m_ignore.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Ignore.
  
  This file contains the implementation of the supergeneric method
  Ignore.
*/

#ifndef m_ignore_h
#define m_ignore_h

#include "agenda_class.h"
#include "exceptions.h"
#include "mystring.h"
#include "workspace_ng.h"

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Ignore(Workspace&,
            // WS Generic Input:
            const Agenda&) {}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Ignore(Workspace&,
            // WS Generic Input:
            const ArrayOfAgenda&) {}

/* Workspace method: Doxygen documentation will be auto-generated */
template <class T>
void Ignore(  // WS Generic Input:
    const T&) {}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void Touch(Workspace& ws,
                  // WS Generic Output:
                  Agenda& out) {
  out = Agenda{ws};
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <class T>
void Touch(  // WS Generic Output:
    T&) {}

#endif  // m_ignore_h
