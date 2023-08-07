#include <workspace.h>

#include <cstdio>

#include "python_interface/pydocs.h"

void groups(std::ostream&) {

}

void variable(std::ostream& os,
              const std::string& name,
              const WorkspaceVariableRecord& wsv) {
  os << "ws.def_property(\"" << name << "\",";
  os << R"--(
    py::cpp_function([](Workspace& w) -> WorkspaceVariable {
      return WorkspaceVariable{w, ")--"
     << name << R"--("};
    }, py::return_value_policy::reference_internal, py::keep_alive<0, 1>()),
    [](Workspace& w, )--"
     << wsv.type << R"--( val) -> void {
      // FIXME!!! FIX AGENDAS
      w.get_or<)--"
     << wsv.type << R"--(>(")--" << name << R"--(") = std::move(val);
    }, py::doc(R"-x-(:class:`~pyarts.arts.WorkspaceVariable` - holds type :class:`~pyarts.arts.)--"
     << wsv.type << "` " << unwrap_stars(wsv.desc) << "\n\n";

  if (wsv.type == "Agenda" or wsv.type == "ArrayOfAgenda") {
    os << get_agenda_io(name);
  }

  if (wsv.default_value) {
    os << "\nDefault value\n"
          "-------------\n\n``"
       << to_defval_str(*wsv.default_value) << "``\n\n";
  }

  os << variable_used_by(name) << '\n';

  os << "\n\nGeneric workspace methods that can generate or use " << name
     << "\n"
     << String(51 + name.size(), '-') << "\n\nSee :class:`~pyarts.arts."
     << wsv.type << "`";
  if (wsv.type not_eq "Any") os << " and/or :class:`~pyarts.arts.Any`\n";
}

void variables(std::vector<std::ofstream>& oss) {
  const auto wsvs = workspace_variables();
  {
    int i = 0;
    for (auto& os : oss)
      os << R"--(#include <python_interface.h

#include <pybind11/functional.h>

namespace Python {
void void py_auto_workspace_wsv_)--"
         << i++
         << "py::class_<Workspace, std::shared_ptr<Workspace>>& ws [[maybe_unused]]) {";
  }

  auto ptr = oss.begin();
  for (auto& [name, wsv] : wsvs) {
    auto& os = *ptr;
    
    if (++ptr == oss.end()) ptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}  // namespace Python\n";
}

void methods(std::ostream& os) {
  const auto wsms = internal_workspace_methods();

  for (const auto& [name, method] : wsms) {
    
  }
}

int main() {
  std::vector<std::ofstream> var_oss(1);
  var_oss[0] = std::ofstream("test.cpp");
  variables(var_oss);
}
