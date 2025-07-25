#pragma once

#include <optional>

#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

struct SetWsv {
  std::string name;
  std::optional<std::string> other{std::nullopt};
  std::optional<Wsv> wsv{std::nullopt};

  SetWsv(std::string n);
  SetWsv(const char* const n) : SetWsv(std::string{n}) {}
  SetWsv(std::string n, WorkspaceGroup auto wsv_value)
      : name(std::move(n)), wsv(std::move(wsv_value)) {}
  SetWsv(std::string n, const char* wsv_name)
      : name(std::move(n)), other(wsv_name) {}
};

class AgendaCreator {
  Agenda a;

 public:
  AgendaCreator(std::string name) : a(std::move(name)) {}
  AgendaCreator(AgendaCreator&&)      = delete;
  AgendaCreator(const AgendaCreator&) = delete;

  AgendaCreator& set(const std::string& name, WorkspaceGroup auto v) {
    a.add(Method{name, std::move(v)});
    return *this;
  }

  AgendaCreator& add(const std::string& name, std::vector<SetWsv>&& v);

  AgendaCreator& ignore(const std::string& name);

  template <std::convertible_to<SetWsv>... T>
  AgendaCreator& add(const std::string& name, T&&... v) {
    return add(name, std::vector<SetWsv>{std::forward<T>(v)...});
  }

  //! Finalize from rvalue.  fix = true will add missing variables to the agenda.
  Agenda finalize(bool fix) &&;
};
