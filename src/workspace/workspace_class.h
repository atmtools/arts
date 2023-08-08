#pragma once

#include <any>
#include <exception>
#include <iomanip>
#include <memory>
#include <optional>
#include <unordered_map>

#include "auto_wsg.h"


enum class WorkspaceInitialization : bool {FromGlobalDefaults, Empty};

class Workspace {
  std::unordered_map<std::string, std::shared_ptr<Wsv>> wsv;

 public:
  Workspace(WorkspaceInitialization how_to_initialize = WorkspaceInitialization::FromGlobalDefaults);

  //! Returns a shared pointer to the workspace variable with the given name.
  [[nodiscard]] std::shared_ptr<Wsv> share(const std::string& name) const;

  //! Returns a copy of the workspace variable with the given name.
  [[nodiscard]] std::shared_ptr<Wsv> copy(const std::string& name) const;

  //! Sets the workspace variable with the given name to the given value.
  void set(const std::string& name, const std::shared_ptr<Wsv>& data);

  //! Copy the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, const T& data) {
    set(name, std::make_shared<Wsv>(data));
  }

  //! Move the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, T&& data) {
    set(name, std::make_shared<Wsv>(std::move(data)));
  }

  //! Borrows the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, T* data) {
    set(name, std::make_shared<Wsv>(data));
  }

  //! Overwrites a variable with another variable of the same type
  void overwrite(const std::string& name, const std::shared_ptr<Wsv>& data);

  //! Copy the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, const T& data) {
    overwrite(name, std::make_shared<Wsv>(data));
  }

  //! Move the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, T&& data) {
    overwrite(name, std::make_shared<Wsv>(std::move(data)));
  }

  //! Borrows the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, T* data) {
    overwrite(name, std::make_shared<Wsv>(data));
  }

  //! Returns a type directly based on the name of the workspace variable.
  template <WorkspaceGroup T>
  [[nodiscard]] T& get(const std::string& name) const try {
    return wsv.at(name) -> get<T>();
  } catch (std::out_of_range&) {
    throw std::runtime_error(var_string("Undefined workspace variable ", std::quoted(name)));
  } catch (std::exception& e) {
    throw std::runtime_error(var_string("Error getting workspace variable ", std::quoted(name), ":\n", e.what()));
  }

  //! Returns a type directly based on the name of the workspace variable, creating it in-place if it is not there
  template <WorkspaceGroup T>
  [[nodiscard]] T& get_or(const std::string& name) {
    if (auto ptr = wsv.find(name); ptr not_eq wsv.end()) {
      return ptr->second->get<T>();
    }

    std::shared_ptr<Wsv> out;
    if constexpr (std::is_same_v<T, Agenda>) {
      out = std::make_shared<Wsv>(T{name});
    } else {
      out = std::make_shared<Wsv>(T{});
    }

    set(name, out);
    return out->get<T>();
  }

  [[nodiscard]] bool contains(const std::string& name) const;

  friend std::ostream& operator<<(std::ostream& os, const Workspace& ws);

  [[nodiscard]] auto begin() const { return wsv.begin(); }

  [[nodiscard]] auto end() const { return wsv.end(); }

  void empty(const std::string& name);
};
