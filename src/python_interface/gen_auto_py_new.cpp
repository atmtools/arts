#include <workspace.h>

#include <cstdio>
#include <fstream>
#include <ostream>

#include "python_interface/pydocs.h"
#include "workspace_methods.h"

std::string variable(const std::string& name,
                     const WorkspaceVariableRecord& wsv) try {
  std::ostringstream os;

  os << "  ws.def_property(\"" << name << "\",";
  os << R"--(
    py::cpp_function([](Workspace& w) -> )--" << wsv.type << R"--(& {
      return w.get_or<)--"
     << wsv.type << R"--(>(")--" << name << R"--(");
    }, py::return_value_policy::reference_internal, py::keep_alive<0, 1>()),
    [](Workspace& w, )--"
     << wsv.type << R"--( val) -> void {
      auto& ws_val = w.get_or<)--"
     << wsv.type << R"--(>(")--" << name << R"--(");
      ws_val = std::move(val);
)--";

  if (wsv.type == "Agenda")
    os << "      ws_val.set_name(\"" << name
       << "\");\n      ws_val.finalize();\n";
  if (wsv.type == "ArrayOfAgenda")
    os << "      for (auto&ws_value: ws_val) {\n        ws_value.set_name(\""
       << name << "\");\n        ws_value.finalize();\n      }\n";

  os << "    }, py::doc(R\"-x-(:class:`~pyarts.arts."
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

  return fix_newlines(os.str()) + ")-x-\"));\n\n";
} catch (std::exception& e) {
  std::cerr << "Error in variable " << std::quoted(name) << ":\n"
            << e.what() << '\n';
  std::exit(1);
}

void variables(std::vector<std::ofstream>& oss) {
  const auto& wsvs = workspace_variables();

  for (std::size_t i = 0; i < oss.size(); i++) {
    oss[i]
        << R"--(#include <python_interface.h>

namespace Python {
void py_auto_workspace_wsv_)--"
        << i++
        << "(py::class_<Workspace, std::shared_ptr<Workspace>>& ws [[maybe_unused]]) {\n";
  }

  auto ptr = oss.begin();
  for (auto& [name, wsv] : wsvs) {
    *ptr << variable(name, wsv) << std::flush;
    if (ptr == oss.end()) ptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}  // namespace Python\n";
}

std::string method(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  os << "  ws.def(\"" << name << "\",";
  os << "nullptr";
  return os.str() + ");\n\n";
}

void methods(std::vector<std::ofstream>& oss) {
  const auto wsms = internal_workspace_methods();

  for (std::size_t i = 0; i < oss.size(); i++) {
    oss[i]
        << R"--(#include <python_interface.h>

namespace Python {
void py_auto_workspace_wsv_)--"
        << i++
        << "(py::class_<Workspace, std::shared_ptr<Workspace>>& ws [[maybe_unused]]) {\n";
  }

  auto ptr = oss.begin();
  for (auto& [name, wsv] : wsms) {
    *ptr << method(name, wsv) << std::flush;
    if (ptr == oss.end()) ptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}  // namespace Python\n";
}

int main() {
  std::vector<std::ofstream> var_oss(1);
  var_oss[0] = std::ofstream("test.cpp");
  variables(var_oss);
  
  std::vector<std::ofstream> met_oss(1);
  met_oss[0] = std::ofstream("test2.cpp");
  methods(met_oss);
}
