#include <workspace.h>

#include <iostream>

#include "pydocs.h"

std::vector<std::string> errors;
#define ERRORAPPEND                \
  catch (std::exception & e) {     \
    errors.emplace_back(e.what()); \
  }

std::ofstream& select_ofstream(std::vector<std::ofstream>& ofs, int i) {
  return ofs[i % ofs.size()];
}

std::string share_type(const std::string& name) {
  const auto& wsgs = internal_workspace_groups();
  if (wsgs.at(name).value_type) return "ValueHolder<" + name + ">";
  return "std::shared_ptr<" + name + ">";
}

std::string variable(const std::string& name,
                     const WorkspaceVariableRecord& wsv) try {
  const auto& wsgs = internal_workspace_groups();
  std::ostringstream os;

  os << "  ws.def_prop_rw(\"" << name << "\",";
  os << R"--(
    [](Workspace& w) -> )--"
     << share_type(wsv.type) << R"--( {
      return w.share_or<)--"
     << wsv.type << ">(\"" << name << R"--(");
    },
    [](Workspace& w, )--"
     << share_type(wsv.type) << R"--( val) -> void {
      w.set(")--"
     << name << R"--(", Wsv{std::move(val)--";
  if (wsgs.at(wsv.type).value_type) os << ".val";
  os << R"--()});
)--";

  if (wsv.type == "Agenda")
    os << "      auto& ws_val = w.get<" << wsv.type << ">(\"" << name
       << "\");\n      ws_val.set_name(\"" << name
       << "\");\n      ws_val.finalize();\n";

  os << "    }, R\"-x-(" << unwrap_stars(wsv.desc) << "\n\n";

  if (wsv.type == "Agenda") {
    os << get_agenda_io(name);
  }

  if (wsv.default_value) {
    os << "\n.. rubric:: Default value\n\n"
       << to_defval_str(*wsv.default_value, "``"sv) << "\n\n";
  }

  os << variable_used_by(name) << '\n';

  return fix_newlines(os.str()) + "\n.. :class:`~pyarts3.arts." + wsv.type +
         "`\n)-x-\");\n\n";
} catch (std::exception& e) {
  std::cerr << "Error in variable " << '"' << name << '"' << ":\n"
            << e.what() << '\n';
  std::exit(1);
}

void variables(const int nfiles) {
  const auto& wsvs = workspace_variables();

  std::vector<std::ofstream> ofs(nfiles);
  for (int i = 0; i < nfiles; i++) {
    ofs[i].open(std::format("py_auto_wsv_{}.cpp", i));
  }

  for (int i = 0; i < nfiles; i++) {
    select_ofstream(ofs, i)
        << R"--(#include <python_interface.h>

#include <nanobind/stl/shared_ptr.h>

namespace Python {
void py_auto_wsv_)--"
        << i << "(py::class_<Workspace>& ws [[maybe_unused]]) {\n";
  }

  int ifile = 0;
  for (auto& [name, wsv] : wsvs) {
    try {
      select_ofstream(ofs, ifile++) << variable(name, wsv) << std::flush;
    }
    ERRORAPPEND;
  }

  for (int i = 0; i < nfiles; i++) {
    select_ofstream(ofs, i) << "}\n}  // namespace Python\n";
  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <NUM VARIABLE FILES>\n";
    return 1;
  }

  const int num_variables = std::stoi(argv[1]);
  variables(num_variables);

  if (errors.size()) {
    std::cerr << "Errors (" << errors.size() << "):\n";
    for (auto& e : errors) {
      std::cerr << e << '\n';
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
