#include "workspace_class.h"

#include <iomanip>
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

const std::shared_ptr<Wsv>& Workspace::share(const std::string& name) const try {
  return wsv.at(name);
} catch (std::out_of_range&) {
  throw std::runtime_error(
      var_string("Undefined workspace variable ", std::quoted(name)));
}

std::shared_ptr<Wsv> Workspace::copy(const std::string& name) const {
  return std::make_shared<Wsv>(share(name) -> copy());
}

void Workspace::set(const std::string& name,
                    const std::shared_ptr<Wsv>& data) try {
  auto ptr = wsv.find(name);

  if (ptr == wsv.end()) {
    if (auto wsv_ptr = wsv_data.find(name);
        wsv_ptr != wsv_data.end() and wsv_ptr->second.type != data->type_name())
      throw wsv_ptr->second.type;

    wsv[name] = data;
  } else {
    std::visit(
        [&data](auto& v) {
          // std::get may throw std::bad_variant_access
          *v = *std::get<std::remove_cvref_t<decltype(v)>>(data->value);
        },
        ptr->second->value);
  }
} catch (const std::bad_variant_access&) {
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

  const std::vector<std::string> names = [](const Workspace& w) {
    std::vector<std::string> n;
    n.reserve(w.wsv.size());
    for (const auto& [name, _] : w.wsv) n.push_back(name);
    std::sort(n.begin(), n.end());
    return n;
  }(ws);

  for (auto& name : names) {
    const Wsv& wsv = *ws.wsv.at(name);

    os << "\n  " << name << " : " << wsv.type_name() << " = ";

    std::string varvalue =
        std::visit([](auto& val) { return var_string(*val); }, wsv.value);
    constexpr std::size_t len = 50;
    const bool is_long = varvalue.size() > len;
    if (is_long) varvalue = varvalue.substr(0, len);
    std::replace(varvalue.begin(), varvalue.end(), '\n', ' ');

    os << varvalue;
    if (is_long) os << "...";
  }
  return os;
}

bool Workspace::contains(const std::string& name) const {
  return wsv.contains(name);
}

void Workspace::init(const std::string& name) try {
  set(name, std::make_shared<Wsv>(Wsv::from_named_type(workspace_variables().at(name).type)));
} catch (std::out_of_range&) {
  throw std::runtime_error(
      var_string("Undefined workspace variable ", std::quoted(name)));
}

Workspace Workspace::deepcopy() const {
  Workspace ws = *this;
  for (auto& [_, value] : ws.wsv) {
    value = std::make_shared<Wsv>(value->copy());
  }
  return ws;
}
