/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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
  ws.del(ws.WsvMap.find(x_name)->second);
}

#endif  // m_ignore_h
