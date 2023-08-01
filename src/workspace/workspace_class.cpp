#include "workspace_class.h"

#include <iomanip>

#include "auto_wsv.h"

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
  return std::make_shared<Wsv>(*share(name));
}

std::shared_ptr<Wsv> Workspace::copy_type(const std::string& name) const {
  return std::make_shared<Wsv>(share(name)->copy_type());
}

void Workspace::set(const std::string& name, const std::shared_ptr<Wsv>& data) {
  auto ptr = wsv.find(name);

  if (ptr == wsv.end()) {
    if (auto wsv_ptr = wsv_data.find(name);
        wsv_ptr != wsv_data.end() and wsv_ptr->second.type != data->type_name())
      throw std::runtime_error(
          var_string("Cannot set built-in workspace variable ",
                     std::quoted(name),
                     " of workspace group ",
                     std::quoted(wsv_ptr->second.type),
                     " to ",
                     std::quoted(data->type_name())));
    wsv[name] = data;
  } else {
    if (not ptr->second->holds_same(*data))
      throw std::runtime_error(
          var_string("Cannot set existing workspace variable ",
                     std::quoted(name),
                     " of workspace group ",
                     std::quoted(ptr->second->type_name()),
                     " to ",
                     std::quoted(data->type_name())));
    ptr->second = data;
  }
}

std::ostream& operator<<(std::ostream& os, const Workspace&) {
  os << "Workspace: PLACEHOLDER\n";
  return os;
}
