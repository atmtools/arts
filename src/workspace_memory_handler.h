/* Copyright (C) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

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

/** The WorkspaceMemoryHandler
 *
 * @file   workspace_memory_handler.h
 * @author Simon Pfreundschuh
 * @date   2020-09-30
 */
#ifndef __ARTS_WORKSPACE_MEMORY_HANDLER__
#define __ARTS_WORKSPACE_MEMORY_HANDLER__
#include "arts.h"
#include <memory>
#include <vector>

/** Handling of workspace memory.
 *
 * The WorkspaceMemoryHandler class is the interface between the workspace
 * and the implementations of the WS group classes. It dispatches calls
 * to create, destroy and duplicate workspace variables to the corresponding
 * classes.
 */
class WorkspaceMemoryHandler {
 public:
  WorkspaceMemoryHandler() {};

  /** Allocate workspace WSV of given group.
   * @param group_index: Index of the group to allocate.
   * @return Void pointer to newly allocated group instance.
   */
  std::shared_ptr<void> allocate(Index group_index) { return allocation_ptrs_[group_index](); }

  /** Duplicate workspace variable of given group.
     * @param group_index The index of the group of the WSV.
     * @param Pointer to the WSV.
      WSV group with the given Index.
  */
  std::shared_ptr<void> duplicate(Index group_index, const std::shared_ptr<void>& ptr) {
    return duplication_ptrs_[group_index](ptr);
  }

  void initialize();

 private:
  std::vector<std::shared_ptr<void> (*)()> allocation_ptrs_;
  std::vector<std::shared_ptr<void> (*)(const std::shared_ptr<void>&)> duplication_ptrs_;
};
#endif
