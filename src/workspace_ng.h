/* Copyright (C) 2004-2012 Oliver Lemke <olemke@core-dump.info>
  
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
#include <map>

class Workspace;

#include "array.h"
#include "wsv_aux.h"

//! Workspace class
/*!
  Manages the workspace variables.
*/
class Workspace {
protected:
  struct WsvStruct {
    void *wsv;
    bool initialized;
    bool auto_allocated;
  };

  //! Workspace variable container.
  Array< stack<WsvStruct *> > ws;

public:
#ifndef NDEBUG
  //! Only for debugging
  String context;
#endif

  static Array<WsvRecord> wsv_data;

  /*! The map associated with wsv_data. */
  static map<String, Index> WsvMap;

  Workspace ();
  Workspace (const Workspace& workspace);
  virtual ~Workspace ();

  static void  define_wsv_data();
  static void  define_wsv_map();
  static Index add_wsv (const WsvRecord& wsv);

  void del (Index i);

  void duplicate (Index i);

  void initialize ();

  //! Checks existence of the given WSV.
  bool is_initialized (Index i) {
    return ((ws[i].size () != 0)
            && (ws[i].top()->initialized == true)); }

  //! Return scoping level of the given WSV.
  Index depth (Index i) {
      return (Index)ws[i].size();
  }

  void *pop (Index i);

  void pop_free (Index i);

  void push (Index i, void *wsv);

  void push_uninitialized (Index i, void *wsv);

  Index nelem () {return ws.nelem ();}

  void *operator[](Index i);
};


//! Print WSV name to output stream.
/** Looks up the name of the WSV with index i and
 prints it to the given output stream.
 
 \param outstream OutputStream
 \param i Index of WSV
 */
template <typename OutputStream> void
PrintWsvName (OutputStream& outstream, Index i)
{
  outstream << Workspace::wsv_data[i].Name () << "(" << i << ") ";
}


//! Print list of WSV names to output stream.
/** Runs through the list of WSV indexes and print all names
 to the given output stream. The list of indexes can be any
 STL container such as Array, vector...
 
 \param outstream OutputStream
 \param container List of WSV indexes
 */
template <typename OutputStream, typename Container> void
PrintWsvNames (OutputStream& outstream, const Container& container)
{
  for (typename Container::const_iterator it = container.begin ();
       it != container.end (); it++ )
  {
    PrintWsvName (outstream, *it);
  }
  
}


#endif /* WORKSPACE_NG_INCLUDED */

