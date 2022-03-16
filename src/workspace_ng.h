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

/** This file contains the Workspace class.
 *
 * @file   workspace_ng.h
 * @author Oliver Lemke <olemke@core-dump.info>
 * @date   2004-11-05
 */

#ifndef WORKSPACE_NG_INCLUDED
#define WORKSPACE_NG_INCLUDED

#include <map>
#include <stack>

class Workspace;

#include "array.h"
#include "wsv_aux.h"

/** Workspace class.
 *
 * Manages the workspace variables.
 */
class Workspace {
 protected:
  struct WsvStruct {
    void *wsv;
    bool initialized;
    bool auto_allocated;
  };

  /** Workspace variable container. */
  Array<stack<WsvStruct *> > ws;

 public:
#ifndef NDEBUG
  /** Debugging context. */
  String context;
#endif

  /** Global WSV data. */
  static Array<WsvRecord> wsv_data;

  /** Global map associated with wsv_data. */
  static map<String, Index> WsvMap;

  /** Construct a new workspace
   *
   * Create the stacks for the WSVs.
   */
  Workspace();

  /** Workspace copy constructor.
   *
   * Make a copy of a workspace. The copy constructor will only copy the topmost
   * layer of the workspace variable stacks.
   *
   * @param[in] workspace The workspace to be copied
   */
  Workspace(const Workspace &workspace);

  /** Destruct the workspace and free all WSVs. */
  virtual ~Workspace();

  /** Define workspace variables. */
  static void define_wsv_data();

  /** Map WSV names to indices. */
  static void define_wsv_map();

  /** Append a new WSV to the workspace. */
  static Index add_wsv(const WsvRecord &wsv);

  /** Delete WSV.
   *
   * Frees the memory of the topmost WSV on the stack.
   *
   * @param[in] i WSV index.
   */
  void del(Index i);

  /** Duplicate WSV.
   *
   * Create another level of scope by duplicating the top element on the WSV
   * stack.
   *
   * @param[in] i
   */
  void duplicate(Index i);

  /** Reset the size of the workspace.
   *
   * Resize the workspace to match the number of WSVs in wsv_data.
   */
  void initialize() { ws.resize(wsv_data.nelem()); }

  /** Checks existence of the given WSV.
   *
   * @param[in] i WSV index.
   * @return true if the WSV exists, otherwise false.
   */
  bool is_initialized(Index i) const {
    return ((ws[i].size() != 0) && (ws[i].top()->initialized == true));
  }

  /** Return scoping level of the given WSV. */
  Index depth(Index i) { return (Index)ws[i].size(); }

  /** Remove the topmost WSV from its stack.
   *
   * Memory is not freed.
   *
   * @see pop_free
   *
   * @param[in] i WSV index.
   */
  void *pop(Index i);

  /** Remove the topmost WSV from its stack and free its memory.
   *
   * @see pop
   *
   *  @param[in] i WSV index.
   */
  void pop_free(Index i);

  /** Push a new WSV onto its stack.
   *
   * @see push_uninitialized
   *
   * @param[in] i WSV index.
   * @param[in] Pointer to variable that should be put on the stack.
   */
  void push(Index i, void *wsv);

  /** Put a new WSV onto its stack.
   *
   * @see push
   *
   * @param[in] i WSV index.
   * @param[in] Pointer to the WSV to put on the stack.
   */
  void push_uninitialized(Index i, void *wsv);

  /** Get the number of workspace variables. */
  Index nelem() const { return ws.nelem(); }
  
  /** Add a new variable to existing workspace and to the static maps */
  Index add_wsv_inplace(const WsvRecord &wsv);

  /** Retrieve a pointer to the given WSV. */
  void *operator[](Index i);

  /** Swap with another workspace */
  void swap(Workspace& other);
};

/** Print WSV name to output stream.
 *
 * Looks up the name of the WSV with index i and
 * prints it to the given output stream.
 *
 * @param[in,out] outstream OutputStream
 * @param[in] i Index of WSV
 */
template <typename OutputStream>
void PrintWsvName(OutputStream &outstream, Index i) {
  outstream << Workspace::wsv_data[i].Name() << "(" << i << ") ";
}

/** Print list of WSV names to output stream.
  *
  * Runs through the list of WSV indexes and print all names to the given output
  * stream. The list of indexes can be any STL container such as Array,
  * vector...
  * @param outstream OutputStream
  * @param container List of WSV indexes
  */
template <typename OutputStream, typename Container>
void PrintWsvNames(OutputStream &outstream, const Container &container) {
  for (typename Container::const_iterator it = container.begin();
       it != container.end();
       it++) {
    PrintWsvName(outstream, *it);
  }
}

#endif /* WORKSPACE_NG_INCLUDED */
