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
  \file   workspace_ng.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2004-11-05

  \brief This file contains the implementation of the workspace
         member functions.

*/

#include "workspace_ng.h"
#include "auto_wsv.h"

//! Construct a new workspace
/*!
  Initialize the stacks for the WSVs. A NULL element is added to each
  stack.
*/
Workspace::Workspace () : EMPTY_WSV (NULL)
{
  ws.resize (N_WSV);

  for (int i = 0; i < N_WSV; i++)
    {
      ws[i].push (EMPTY_WSV);
    }
}

//! Destruct the workspace
/*!
  Frees all WSVs.
*/
Workspace::~Workspace ()
{
  for (int i = 0; i < ws.nelem (); i++)
    {
      void *vp;

      // Store top of stack in vp
      while ((vp = ws[i].top ()))
        {
          ws[i].pop ();
          // And free memory if vp is not NULL.
          wsmh.deallocate (i, vp);
        }
    }
  ws.empty ();

}

//! Pop the topmost wsv from its stack.
/*!
  Removes the topmost element from the wsv's stack.
  If necessary, the calling function has to free the wsv's memory.
 */
void *Workspace::pop (Index i)
{
  void *vp = ws[i].top ();
  ws[i].pop ();
  return vp;
}

//! Push a new wsv onto its stack.
/*! Adds the pointer to the variable to the wsv stack i. */
void Workspace::push (Index i, void *wsv)
{
  ws[i].push (wsv);
}

//! Retrieve pointer to the given WSV.
/*!
  This method returns a void pointer to the topmost instance of the
  given workspace variable.
*/
void *Workspace::operator[](Index i)
{
  if (!ws[i].top ())
    {
      ws[i].push (wsmh.allocate (i));
    }

  return (ws[i].top());
}


