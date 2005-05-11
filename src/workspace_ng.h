/* Copyright (C) 2004 Oliver Lemke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   workspace_ng.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2004-11-05

  \brief This file contains the declaration and partly the implementation
         of the workspace class.

*/


#ifndef WORKSPACE_NG_INCLUDED
#define WORKSPACE_NG_INCLUDED

#include <stack>

#include "array.h"
#include "wsv_aux.h"
#include "auto_wsv.h"

//! Workspace class
/*!
  Manages the workspace variables.
*/
class Workspace {
private:
  void *EMPTY_WSV;

  //! Workspace variable container.
  Array< stack< void * > > ws;

  //! Memory handler for allocation and deallocation of WSVs.
  WorkspaceMemoryHandler wsmh;

public:
  Workspace ();
  virtual ~Workspace ();

  //! Checks existence of the given WSV.
  bool is_occupied(Index i) { return (ws[i].top() != NULL); }

  void *operator[](Index i);
};

#endif /* WORKSPACE_NG_INCLUDED */

