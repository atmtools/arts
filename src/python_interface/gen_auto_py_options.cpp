#include <arts_options.h>
#include <nonstd.h>
#include <workspace.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <set>

#include "pydocs.h"

void enum_option(std::ostream& os, const EnumeratedOption& wso) {
  os << "  py::class_<" << wso.name << "> _g" << wso.name << "(m, \""
     << wso.name << "\");\n";

  os << "  xml_interface(_g" << wso.name << ");\n";
  os << "  str_interface(_g" << wso.name << ");\n";

  os << "  _g" << wso.name << ".def(py::init<>())\n";

  os << "      .def(py::init<" << wso.name << ">())\n";

  os << "      .def(\"__init__\", [](" << wso.name
     << "*y, const std::string& x) {new (y) " << wso.name << "{to<" << wso.name
     << ">(x)};}, \"String constructor\")\n";

  os << "      .def(\"__hash__\", [](const " << wso.name
     << "& x) {return std::hash<" << wso.name << ">{}(x);})\n";

  os << "      .def(\"__copy__\", [](" << wso.name << " t) -> " << wso.name
     << " {return t;})\n";

  os << "      .def(\"__deepcopy__\", [](" << wso.name << " t, py::dict&) -> "
     << wso.name << " { return t; })\n";

  os << "      .def(py::self == py::self)\n";
  os << "      .def(py::self != py::self)\n";
  os << "      .def(py::self <= py::self)\n";
  os << "      .def(py::self >= py::self)\n";
  os << "      .def(py::self < py::self)\n";
  os << "      .def(py::self > py::self)\n";

  os << "      .def(\"__getstate__\",\n        [](" << wso.name
     << "& t) {\n"
        "          return std::tuple<std::string>{String{toString(t)}};\n"
        "      })\n"
        "      .def(\"__setstate__\",\n        []("
     << wso.name
     << "* e, const std::tuple<std::string>& state) {\n"
        "           new (e) "
     << wso.name << "{to<" << wso.name
     << ">(std::get<0>(state))};\n"
        "      })\n";

  os << "      .def_static(\"get_options\", [](){return enumtyps::" << wso.name
     << "Types;}, \"Get a list of all options\")\n";
  os << "      .def_static(\"get_options_as_strings\", [](){return enumstrs::"
     << wso.name << "Names<>;}, \"Get a list of all options as strings\")\n";

  static std::array pykeywords{"None", "any", "all", "print"};
  constexpr std::string_view ignore_str =
      "-+={§±!@#$%^&*()-+=]}[{\\|'\";:?/.>,<`~}] ";
  static std::set<char> ignore(ignore_str.begin(), ignore_str.end());
  auto contains_invalid_chars = [](const std::string& str) {
    for (auto ch : str)
      if (std::ranges::binary_search(ignore, ch)) return true;
    return false;
  };
  for (auto& value : wso.values_and_desc) {
    // Skip last element in value which contains the description
    for (Size i = 0; i < value.size() - 1; i++) {
      auto& x = value[i];
      if (nonstd::isdigit(x.front()) or contains_invalid_chars(x) or
          std::ranges::any_of(value | std::views::take(i), Cmp::eq(x)))
        continue;

      os << "      .def_prop_ro_static(\"" << x;

      if (std::ranges::any_of(pykeywords, Cmp::eq(x))) {
        os << '_';
      }

      // .front is the actual value, .back the description
      os << "\", [](py::object&){return " << wso.name << "::" << value.front()
         << ";}, R\"-x-(" << value.back() << ")-x-\")\n";
    }
  }

  os << "      .doc() = R\"-x-(" << unwrap_stars(wso.docs()) << ")-x-\";\n";

  os << "  py::implicitly_convertible<std::string, " << wso.name << ">();\n\n";
}

void enum_options(const std::string& fname) {
  const auto& wsos = internal_options();

  auto cc = std::ofstream(fname + ".cpp");
  auto hh = std::ofstream(fname + ".h");

  cc << R"-x-(#include "python_interface.h"
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string.h>

#include <hpy_arts.h>

#include")-x-"
     << fname << ".h\"\n\n";

  hh << "#pragma once\n\n"
        "#include <enums.h>\n"
        "#include <nanobind/nanobind.h>\n\n";
  for (auto& wso : wsos) {
    hh << "NB_MAKE_OPAQUE(" << wso.name << ");\n";
  }

  cc << R"-x-(

namespace Python {
void py_auto_options(py::module_& m) try {
)-x-";

  for (auto& opt : wsos) {
    enum_option(cc, opt);
  }
  cc << R"-x-(} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize automatic options\n", e.what()));
}
}  // namespace Python
)-x-";
}

int main() { enum_options("py_auto_options"); }
