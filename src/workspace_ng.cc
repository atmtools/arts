/* Copyright (C) 2004 Oliver Lemke <olemke@core-dump.info>
  
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
  \file   workspace_ng.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2004-11-05

  \brief This file contains the implementation of the workspace
         member functions.

*/

#include "workspace_ng.h"
#include "wsv_aux.h"
#include "auto_wsv.h"

WorkspaceMemoryHandler wsmh;

//! Construct a new workspace
/*!
  Create the stacks for the WSVs.
*/
Workspace::Workspace () : ws(N_WSV)
{
}

//! Workspace copy constructor
/*!
  Make a copy of a workspace. The copy constructor will only copy the topmost
  layer of the workspace variable stacks.

  \param[in] workspace The workspace to be copied

  \author Oliver Lemke
  \date   2007-11-28
*/
Workspace::Workspace (const Workspace& workspace) : ws(N_WSV)
{
  for (Index i=0; i < workspace.ws.nelem(); i++)
    {
      WsvStruct *wsvs = new WsvStruct;
      wsvs->auto_allocated = false;
      if (workspace.ws[i].size() && workspace.ws[i].top()->wsv)
        {
          wsvs->wsv = workspace.ws[i].top()->wsv;
          wsvs->initialized = workspace.ws[i].top()->initialized;
        }
      else
        {
          wsvs->wsv = NULL;
          wsvs->initialized = false;
        }
      ws[i].push (wsvs);
    }
}

//! Destruct the workspace
/*!
  Frees all WSVs.
*/
Workspace::~Workspace ()
{
  extern const Array<WsvRecord> wsv_data;

  for (int i = 0; i < ws.nelem (); i++)
    {
      WsvStruct *wsvs;

      while (ws[i].size ())
        {
          wsvs = ws[i].top ();
          if (wsvs->auto_allocated && wsvs->wsv)
            {
              wsmh.deallocate (wsv_data[i].Group(), wsvs->wsv);
            }
          delete (wsvs);
          ws[i].pop ();
        }
    }
  ws.empty ();

}

//! Duplicate WSV.
/*!
  Copies the topmost WSV and puts it back on the WSV stack.

  \param i WSV index.
 */
void Workspace::duplicate (Index i)
{
  extern const Array<WsvRecord> wsv_data;
  WsvStruct *wsvs = new WsvStruct;

  wsvs->auto_allocated = true;
  if (ws[i].size())
    {
      wsvs->wsv = wsmh.duplicate (wsv_data[i].Group(), ws[i].top()->wsv);
      wsvs->initialized = true;
    }
  else
    {
      wsvs->wsv = NULL;
      wsvs->initialized = false;
    }
  ws[i].push (wsvs);
}

//! Pop the topmost wsv from its stack.
/*!
  Removes the topmost element from the wsv's stack.
  If necessary, the calling function has to free the wsv's memory.

  \param i WSV index.
 */
void *Workspace::pop (Index i)
{
  WsvStruct *wsvs = ws[i].top ();
  void *vp = NULL;
  if (wsvs)
    {
      vp = wsvs->wsv;
      delete wsvs;
      ws[i].pop ();
    }
  return vp;
}

//! Pop the topmost wsv from its stack and free its memory.
/*!
  Removes the topmost element from the wsv's stack and frees memory.

  \param i WSV index.
 */
void Workspace::pop_free (Index i)
{
  extern const Array<WsvRecord> wsv_data;
  WsvStruct *wsvs = ws[i].top ();

  if (wsvs)
    {
      if (wsvs->wsv)
        wsmh.deallocate (wsv_data[i].Group(), wsvs->wsv);

      delete wsvs;
      ws[i].pop ();
    }
}

//! Push a new wsv onto its stack.
/*!
  Adds the pointer to the variable to the stack of the WSV with index i.

  \param i WSV index.
  \param wsv Void pointer to variable that should be put on the stack.
  */
void Workspace::push (Index i, void *wsv)
{
  WsvStruct *wsvs = new WsvStruct;
  wsvs->auto_allocated = false;
  wsvs->initialized = true;
  wsvs->wsv = wsv;
  ws[i].push (wsvs);
}

//! Push a new wsv onto its stack but mark it as uninitialized.
/*!
  Adds the pointer to the variable to the stack of the WSV with index i.
  The variable is flagged as uninitialized. This is used for agenda
  output-only variables.

  \param i WSV index.
  \param wsv Void pointer to variable that should be put on the stack.
  */
void Workspace::push_uninitialized (Index i, void *wsv)
{
  WsvStruct *wsvs = new WsvStruct;
  wsvs->auto_allocated = false;
  wsvs->initialized = false;
  wsvs->wsv = wsv;
  ws[i].push (wsvs);
}

//! Retrieve pointer to the given WSV.
/*!
  This method returns a void pointer to the topmost instance of the
  given workspace variable.

  \param i WSV index.
*/
void *Workspace::operator[](Index i)
{
  if (!ws[i].size ())
    push (i, NULL);

  if (!ws[i].top ()->wsv)
    {
      extern const Array<WsvRecord> wsv_data;

      ws[i].top ()->auto_allocated = true;
      ws[i].top ()->wsv = wsmh.allocate (wsv_data[i].Group());
    }

  ws[i].top ()->initialized = true;

  return (ws[i].top()->wsv);
}

