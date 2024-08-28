#include "workspace_class.h"

#include <ranges>
#include <stdexcept>
#include <type_traits>

#include <auto_wsv.h>
#include "debug.h"

const auto& wsv_data = workspace_variables();

Workspace::Workspace(WorkspaceInitialization how_to_initialize) : wsv{} {
  if (WorkspaceInitialization::FromGlobalDefaults == how_to_initialize) {
    for (const auto& [name, record] : wsv_data) {
      if (record.default_value.has_value()) {
        wsv[name] = Wsv{record.default_value.value()};
      }
    }
  }
}

const Wsv& Workspace::share(const std::string& name) const try {
  return wsv.at(name);
} catch (std::out_of_range&) {
  throw std::runtime_error(
      var_string("Undefined workspace variable ", '"', name, '"'));
}

Wsv Workspace::copy(const std::string& name) const {
  return share(name).copy();
}

void Workspace::set(const std::string& name, const Wsv& data) try {
  auto ptr = wsv.find(name);

  if (ptr == wsv.end()) {
    if (auto wsv_ptr = wsv_data.find(name);
        wsv_ptr != wsv_data.end() and wsv_ptr->second.type != data.type_name())
      throw wsv_ptr->second.type;

    wsv[name] = data;
  } else {
    std::visit(
        [&data](auto& v) {
          // std::get may throw std::bad_variant_access
          *v = *std::get<std::remove_cvref_t<decltype(v)>>(data.value);
        },
        ptr->second.value);
  }
} catch (const std::bad_variant_access&) {
  throw std::runtime_error(var_string("Workspace variable ",
                                      '"',
                                      name,
                                      '"',
                                      " is of type ",
                                      '"',
                                      wsv.at(name).type_name(),
                                      '"',
                                      ".\nIt cannot be set to ",
                                      '"',
                                      data.type_name(),
                                      '"'));
} catch (const std::string& type) {
  throw std::runtime_error(var_string("Cannot set built-in workspace variable ",
                                      '"',
                                      name,
                                      '"',
                                      " of workspace group ",
                                      '"',
                                      type,
                                      '"',
                                      " to ",
                                      '"',
                                      data.type_name(),
                                      '"'));
}

void Workspace::overwrite(const std::string& name, const Wsv& data) try {
  if (auto wsv_ptr = wsv_data.find(name);
      wsv_ptr != wsv_data.end() and wsv_ptr->second.type != data.type_name())
    throw wsv_ptr->second.type;
  wsv[name] = data;
} catch (const std::string& type) {
  throw std::runtime_error(var_string("Cannot set built-in workspace variable ",
                                      '"',
                                      name,
                                      '"',
                                      " of workspace group ",
                                      '"',
                                      type,
                                      '"',
                                      " to ",
                                      '"',
                                      data.type_name(),
                                      '"'));
}

std::ostream& operator<<(std::ostream& os, const Workspace& ws) {
  return os << std::format("{:s,}", ws);
}

bool Workspace::contains(const std::string& name) const {
  return wsv.contains(name);
}

bool Workspace::wsv_and_contains(const std::string& name) const {
  ARTS_USER_ERROR_IF(
      std::ranges::none_of(wsv_data | std::views::keys,
                           [&name](const auto& n) { return n == name; }),
      "Invalid workspace variable: ",
      '"',
      name,
      '"')
  return contains(name);
}

void Workspace::init(const std::string& name) try {
  set(name, Wsv::from_named_type(workspace_variables().at(name).type));
} catch (std::out_of_range&) {
  throw std::runtime_error(
      var_string("Undefined workspace variable ", '"', name, '"'));
}

Workspace Workspace::deepcopy() const {
  Workspace ws = *this;
  for (auto& [_, value] : ws.wsv) {
    value = value.copy();
  }
  return ws;
}
