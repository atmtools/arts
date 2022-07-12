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
#include <memory>
#include <stack>
#include <vector>

#include "array.h"
#include "wsv_aux.h"

struct WorkspaceVariableStruct final {
  std::shared_ptr<void> wsv;
  bool initialized;
};

using WorkspaceVariable = stack<WorkspaceVariableStruct, std::vector<WorkspaceVariableStruct>>;

struct WorkspaceBorrowGuard final {
  WorkspaceVariable &wsv;
  ~WorkspaceBorrowGuard() noexcept { wsv.pop(); };
};

/** Workspace class.
 *
 * Manages the workspace variables.
 */
class Workspace final {
 protected:

  /** Workspace variable container. */
  Array<WorkspaceVariable> ws{};

 public:
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
  void set_empty(Index i);

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
  [[nodiscard]] bool is_initialized(Index i) const {
    return ws[i].size() and ws[i].top().initialized;
  }

  /** Return scoping level of the given WSV. */
  Index depth(Index i) { return static_cast<Index>(ws[i].size()); }

  /** Remove the topmost WSV from its stack.
   *
   * Memory is freed depending on if it was borrowed or not.
   *
   * @param[in] i WSV index.
   */
  void pop(Index i);

  void emplace(Index);

  /** Push a new WSV onto its stack.
   *
   * @see push_uninitialized
   *
   * @param[in] i WSV index.
   * @param[in] Pointer to variable that should be put on the stack.
   */
  template <typename T>
  WorkspaceBorrowGuard borrow(Index i, T &wsv) {
    using U = std::remove_cv_t<T>;
    std::shared_ptr<U> wsv_ptr(const_cast<U *>(&wsv), [](U *) {});
    WorkspaceVariableStruct wsvs;
    wsvs.initialized = true;
    wsvs.wsv = wsv_ptr;
    ws[i].push(wsvs);
    return {ws[i]};
  }

  /** Put a new WSV onto its stack.
   *
   * @see push
   *
   * @param[in] i WSV index.
   * @param[in] Pointer to the WSV to put on the stack.
   */
  template <typename T>
  WorkspaceBorrowGuard borrow_uninitialized(Index i, T &wsv) {
    using U = std::remove_cv_t<T>;
    std::shared_ptr<U> wsv_ptr(const_cast<U *>(&wsv), [](U *) {});
    WorkspaceVariableStruct wsvs;
    wsvs.initialized = false;
    wsvs.wsv = wsv_ptr;
    ws[i].push(wsvs);
    return {ws[i]};
  }

  /** Move a WSV onto its stack.
   *
   * @see push_uninitialized
   *
   * @param[in] i WSV index.
   * @param[in] Pointer to variable that should be put on the stack.
   */
  template <typename T>
  void push_move(Index i, std::shared_ptr<T> &&wsv_ptr) {
    WorkspaceVariableStruct wsvs;
    wsvs.initialized = true;
    wsvs.wsv = std::forward<std::shared_ptr<T>>(wsv_ptr);
    ws[i].push(wsvs);
  }

  /** Get the number of workspace variables. */
  [[nodiscard]] Index nelem() const { return ws.nelem(); }

  /** Add a new variable to existing workspace and to the static maps */
  Index add_wsv_inplace(const WsvRecord &wsv);

  /** Retrieve a pointer to the given WSV. */
  std::shared_ptr<void> operator[](Index i);

  /** Retrieve a value ptr if it exist (FIXME: C++20 allows const char* as template argument) */
  template <class T>
  T* get(const char *name) {
    if (const Index pos = WsvMap.at(name); is_initialized(pos)) {
      return static_cast<T*>(ws[pos].top().wsv.get());
    }
    return nullptr;
  }

  /** Swap with another workspace */
  void swap(Workspace &other);
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
