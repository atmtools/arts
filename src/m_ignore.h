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

#include <workspace.h>

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Ignore(const T&) {
  ARTS_TIME_REPORT
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Touch(T&) {
  ARTS_TIME_REPORT
}

#endif  // m_ignore_h
