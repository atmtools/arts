#include "workspace_class.h"

#include <iomanip>
#include <ranges>
#include <stdexcept>
#include <type_traits>

#include "auto_wsv.h"
#include "debug.h"

const auto& wsv_data = workspace_variables();

Workspace::Workspace(WorkspaceInitialization how_to_initialize) : wsv{} {
  if (WorkspaceInitialization::FromGlobalDefaults == how_to_initialize) {
    for (const auto& [name, record] : wsv_data) {
      if (record.default_value.has_value()) {
        wsv[name] = std::make_shared<Wsv>(record.default_value.value());
      }
    }
  }
}

std::shared_ptr<Wsv> Workspace::share(const std::string& name) const try {
  return wsv.at(name);
} catch (std::out_of_range&) {
  throw std::runtime_error(
      var_string("Cannot find workspace variable ", std::quoted(name)));
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Cannot share workspace variable ", std::quoted(name), '\n', e.what()));
}

std::shared_ptr<Wsv> Workspace::copy(const std::string& name) const {
  return std::visit([](auto x){return std::make_shared<Wsv>(std::move(*x));}, share(name) -> value);
}

void Workspace::set(const std::string& name,
                    const std::shared_ptr<Wsv>& data) try {
  auto ptr = wsv.find(name);

  if (ptr == wsv.end()) {
    if (auto wsv_ptr = wsv_data.find(name);
        wsv_ptr != wsv_data.end() and wsv_ptr->second.type != data->type_name())
      throw wsv_ptr->second.type;

    wsv[name] = std::make_shared<Wsv>(data->copy());
  } else {
    std::visit(
        [&data](auto& v) {
          // std::get may throw std::bad_variant_access
          *v = *std::get<std::remove_cvref_t<decltype(v)>>(data->value);
        },
        ptr->second->value);
  }
} catch (const std::bad_variant_access& e) {
  throw std::runtime_error(var_string("Workspace variable ",
                                      std::quoted(name),
                                      " is of type ",
                                      std::quoted(wsv.at(name)->type_name()),
                                      ".\nIt cannot be set to ",
                                      std::quoted(data->type_name())));
} catch (const std::string& type) {
  throw std::runtime_error(var_string("Cannot set built-in workspace variable ",
                                      std::quoted(name),
                                      " of workspace group ",
                                      std::quoted(type),
                                      " to ",
                                      std::quoted(data->type_name())));
}

void Workspace::overwrite(const std::string& name,
                          const std::shared_ptr<Wsv>& data) try {
  if (auto wsv_ptr = wsv_data.find(name);
      wsv_ptr != wsv_data.end() and wsv_ptr->second.type != data->type_name())
    throw wsv_ptr->second.type;
  wsv[name] = data;
} catch (const std::string& type) {
  throw std::runtime_error(var_string("Cannot set built-in workspace variable ",
                                      std::quoted(name),
                                      " of workspace group ",
                                      std::quoted(type),
                                      " to ",
                                      std::quoted(data->type_name())));
}

std::ostream& operator<<(std::ostream& os, const Workspace& ws) {
  os << "Workspace containing:";

  if (ws.wsv.empty()) {
    os << " nothing";
    return os;
  }

  os << "\n  ";

  constexpr std::size_t linebreak = 77;
  auto str = var_string(std::quoted(ws.wsv.begin() -> first), " : ", ws.wsv.begin() -> second ->type_name());
  os << str;
  std::size_t count = str.size();
  for (auto& v: std::ranges::drop_view{ws.wsv, 1}) {
    str = var_string(std::quoted(v.first), " : ", v.second -> type_name());
    count += str.size();
    if (count >= linebreak) {
      os << ",\n  ";
      count = 0;
    } else {
      os << ", ";
    }
    os << str;
  }
  return os;
}

bool Workspace::contains(const std::string& name) const {
  return wsv.contains(name);
}
