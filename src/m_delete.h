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

#include <map>

#include "m_general.h"
#include "mystring.h"
#include "workspace_ng.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void Delete(  // Workspace reference
    Workspace& ws,
    // WS Generic Input:
    const T& x _U_,
    // WS Generic Input Names:
    const String& x_name,
    const Verbosity&) {
  ws.set_empty(ws.WsvMap_ptr->at(x_name));
}

#endif  // m_ignore_h
