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
