#pragma once

#include <memory>
#include <unordered_map>

#include "auto_wsg.h"

enum class WorkspaceInitialization : bool { FromGlobalDefaults, Empty };

class Workspace {
  std::unordered_map<std::string, std::shared_ptr<Wsv>> wsv;

 public:
  Workspace(WorkspaceInitialization how_to_initialize =
                WorkspaceInitialization::FromGlobalDefaults);

  //! Returns a shared pointer to the workspace variable with the given name.
  [[nodiscard]] const std::shared_ptr<Wsv>& share(
      const std::string& name) const;

  //! Returns a copy of the workspace variable with the given name.
  [[nodiscard]] std::shared_ptr<Wsv> copy(const std::string& name) const;

  //! Sets the workspace variable with the given name to the given value.
  void set(const std::string& name, const std::shared_ptr<Wsv>& data);

  //! Copy the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, const T& data);

  //! Move the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, T&& data);

  //! Borrows the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, T* data);

  //! Overwrites a variable with another variable of the same type
  void overwrite(const std::string& name, const std::shared_ptr<Wsv>& data);

  //! Copy the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, const T& data);

  //! Move the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, T&& data);

  //! Borrows the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, T* data);

  //! Returns a type directly based on the name of the workspace variable.
  template <WorkspaceGroup T>
  [[nodiscard]] T& get(const std::string& name) const;

  //! Returns a type directly based on the name of the workspace variable, creating it in-place if it is not there
  template <WorkspaceGroup T>
  [[nodiscard]] const std::shared_ptr<T>& share_or(const std::string& name);

  //! Returns a type directly based on the name of the workspace variable, creating it in-place if it is not there
  template <WorkspaceGroup T>
  [[nodiscard]] T& get_or(const std::string& name);

  [[nodiscard]] bool contains(const std::string& name) const;

  friend std::ostream& operator<<(std::ostream& os, const Workspace& ws);

  [[nodiscard]] auto begin() { return wsv.begin(); }

  [[nodiscard]] auto end() { return wsv.end(); }

  [[nodiscard]] auto begin() const { return wsv.begin(); }

  [[nodiscard]] auto end() const { return wsv.end(); }

  void init(const std::string& name);

  [[nodiscard]] Workspace deepcopy() const;
};
