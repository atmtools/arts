#include <workspace.h>

#include "pydocs.h"
#include "workspace_groups.h"

void default_groups(const std::string& fname) {
  const auto& wsgs = internal_workspace_groups();

  {
    std::ofstream hh(fname + ".h");

    hh << "#pragma once\n\n#include <python_interface.h>\n#include <type_traits>\nnamespace Python {\nnamespace py = pybind11;\n\n";

    for (const auto& [name, wsg] : wsgs) {
      const std::string T = wsg.value_type ? "ValueHolder<" + name + ">" : name;
      const std::string dpyT =
          wsg.is_array
              ? var_string("decltype(artsarray<", T, ">(m, \"", name, "\"))")
              : "artsclass<" + T + ">";

      hh << "auto py_static" << name << "(py::module_& m) -> " << dpyT
         << "&;\n";
    }
    hh << "}  // namespace Python\n";
  }
  {
    std::ofstream cc(fname + ".cpp");
    cc << "#include \"" << fname
       << ".h\"\n#include <py_macros.h>\n\nnamespace Python{\n";
    for (const auto& [name, wsg] : wsgs) {
      const std::string T = wsg.value_type ? "ValueHolder<" + name + ">" : name;
      const std::string pyT = std::string{"arts"} +
                              (wsg.is_array ? "array" : "class") + "<" + T +
                              ">";
      const std::string dpyT =
          wsg.is_array
              ? var_string("decltype(artsarray<", T, ">(m, \"", name, "\"))")
              : "artsclass<" + T + ">";

      cc << "auto py_static" << name << "(py::module_& m) -> " << dpyT
         << "& {\n  static " << dpyT << " g = " << pyT << "(m, \"" << name
         << "\"";
      if (wsg.is_simple_contiguous) {
        cc << ", py::buffer_protocol()";
      }
      cc << ");\n";

      cc << '\n';
      cc << "  static bool done = false;\n  if (done) return g;\n  done = true;\n";

      cc << '\n';
      cc << "  g.def(py::init<>(), py::doc(\"Default constructor\"));\n";
      cc << "  g.def(py::init([](const " << T
         << "& x) { return std::make_shared<" << T << ">(x); }), "
         << "py::arg(\"val\"), "
         << " py::doc(\"Copy instance\"));\n";
      if (wsg.value_type) {
        cc << "  g.def(py::init([](const " << name
           << "& x) { return std::make_shared<" << T << ">(x); }), "
           << "py::arg(\"val\"), "
           << " py::doc(\"Copy instance\"));\n";
      }

      cc << '\n';
      cc << "  g.PythonInterfaceBasicRepresentation(" << T << ");\n";
      cc << "  g.PythonInterfaceCopyValue(" << T << ");\n";
      cc << "  g.PythonInterfaceFileIO2(" << T << ", " << name << ");\n";

      if (wsg.is_simple_contiguous and not wsg.is_array) {
        cc << R"-x-(
  g.def(py::pickle(
    [](const py::object& self) {
      return py::make_tuple(self.attr("value"));
    },
    [](const py::tuple& t) {
      ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
      return py::type::of<)-x-"
           << T << ">()(t[0]).cast<" << T << ">();\n    }\n  ));\n";
      }

      cc << '\n';
      cc << "  g.doc() = R\"-x-("
         << unwrap_stars(var_string(internal_workspace_groups().at(name).desc,
                                    '\n',
                                    Python::group_generics_inout(name),
                                    Python::group_workspace_types(name)))
         << ")-x-\";\n";

      cc << "  return g;\n}\n\n";
    }
    cc << "}  // namespace Python\n";
  }
}