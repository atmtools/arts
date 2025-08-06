#include "workspace_class.h"

#include <auto_wsv.h>

#include <ranges>
#include <stdexcept>
#include <type_traits>

#include "debug.h"

Workspace::Workspace(WorkspaceInitialization how_to_initialize) : wsv{} {
  if (WorkspaceInitialization::FromGlobalDefaults == how_to_initialize) {
    for (const auto& [name, record] : workspace_variables()) {
      if (record.default_value.has_value()) {
        wsv[name] = Wsv{record.default_value.value()};
      }
    }
  }
}

Workspace::Workspace(std::unordered_map<std::string, Wsv> wsv)
    : wsv(std::move(wsv)) {};

const Wsv& Workspace::share(const std::string& name) const try {
  return wsv.at(name);
} catch (std::out_of_range&) {
  throw std::runtime_error(
      std::format("Undefined workspace variable \"{}\"", name));
}

Wsv Workspace::copy(const std::string& name) const {
  return share(name).copied();
}

void Workspace::set(const std::string& name, const Wsv& data) {
  auto ptr = wsv.find(name);

  if (ptr == wsv.end()) {
    if (auto wsv_ptr = workspace_variables().find(name);
        wsv_ptr != workspace_variables().end() and
        wsv_ptr->second.type != data.type_name()) {
      throw std::runtime_error(std::format(
          R"(Cannot set built-in workspace variable "{}" of workspace group "{}" to type "{}")",
          name,
          wsv_ptr->second.type,
          data.type_name()));
    }

    wsv[name] = data;
  } else {
    if (not ptr->second.holds_same(data)) {
      throw std::runtime_error(std::format(
          R"(Workspace variable "{}" is of type "{}". It cannot be set to be a "{}")",
          name,
          ptr->second.type_name(),
          data.type_name()));
    }

    ptr->second = data;
  }
}

void Workspace::overwrite(const std::string& name, const Wsv& data) try {
  if (auto wsv_ptr = workspace_variables().find(name);
      wsv_ptr != workspace_variables().end() and
      wsv_ptr->second.type != data.type_name())
    throw wsv_ptr->second.type;
  wsv[name] = data;
} catch (const std::string& type) {
  throw std::runtime_error(std::format(
      R"(Cannot set built-in workspace variable "{}" of workspace group "{}" to type "{}")",
      name,
      type,
      data.type_name()));
}

std::ostream& operator<<(std::ostream& os, const Workspace& ws) {
  return os << std::format("{:s,}", ws);
}

bool Workspace::contains(const std::string& name) const {
  return wsv.contains(name);
}

bool Workspace::wsv_and_contains(const std::string& name) const {
  ARTS_USER_ERROR_IF(
      std::ranges::none_of(workspace_variables() | std::views::keys,
                           [&name](const auto& n) { return n == name; }),
      "Invalid workspace variable: \"{}\"",
      name)
  return contains(name);
}

void Workspace::init(const std::string& name) try {
  set(name, Wsv::from_named_type(workspace_variables().at(name).type));
} catch (std::out_of_range&) {
  throw std::runtime_error(
      std::format("Undefined workspace variable \"{}\"", name));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error setting '{}'\n{}", name, e.what()));
}

Workspace Workspace::deepcopy() const {
  Workspace ws = *this;
  for (auto& [_, value] : ws.wsv) {
    value = value.copied();
  }
  return ws;
}

void xml_io_stream<Workspace>::write(std::ostream& os,
                                     const Workspace& x,
                                     bofstream* pbofs,
                                     std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.wsv, pbofs, "WSVs"sv);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Workspace>::read(std::istream& is,
                                    Workspace& x,
                                    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.wsv, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
