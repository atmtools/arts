/*!
  \file   m_delete.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2007-11-26
  
  \brief  Implementation of Delete.
  
  This file contains the implementation of the supergeneric method
  Delete.
*/

#ifndef m_delete_h
#define m_delete_h

#include <workspace.h>

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void Delete(T& x) {
  x = T{};
}

#endif  // m_ignore_h
