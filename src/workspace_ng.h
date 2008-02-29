/* Copyright (C) 2004-2007 Oliver Lemke <olemke@core-dump.info>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

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

//! Workspace class
/*!
  Manages the workspace variables.
*/
class Workspace {
private:
  typedef struct {
    void *wsv;
    bool initialized;
    bool auto_allocated;
  } WsvStruct;

  //! Workspace variable container.
  Array< stack<WsvStruct *> > ws;

  //! Memory handler for allocation and deallocation of WSVs.

public:
  Workspace ();
  Workspace (const Workspace& workspace);
  virtual ~Workspace ();

  void initialize ();

  //! Checks existence of the given WSV.
  bool is_initialized (Index i) {
    return ((ws[i].size () != 0)
            && (ws[i].top()->initialized == true)); }

  void duplicate (Index i);

  void *pop (Index i);

  void pop_free (Index i);

  void push (Index i, void *wsv);

  void push_uninitialized (Index i, void *wsv);

  Index nelem () {return ws.nelem ();}

  void *operator[](Index i);
};

#endif /* WORKSPACE_NG_INCLUDED */

