#include <workspace.h>
#include <arts_options.h>
#include <iostream>
#include <set>
#include <nonstd.h>
#include <algorithm>
#include <ranges>

#include "pydocs.h"

void enum_option(std::ostream& os, const EnumeratedOption& wso) {
  os << "  py::class_<" << wso.name << "> _g"<< wso.name <<"(m, \"" << wso.name << "\");\n  _g"<< wso.name;

  os << ".def(py::init<>())\n";

  os << "      .def(py::init<"<< wso.name <<">())\n";

  os << "      .def(\"__init__\", []("<< wso.name <<"*y, const std::string& x) {new (y) "<< wso.name <<"{to<" << wso.name
     << ">(x)};}, \"String constructor\")\n";

  os << "      .def(\"__hash__\", [](const " << wso.name
     << "& x) {return std::hash<" << wso.name << ">{}(x);})\n";

  os << "      .def(\"__copy__\", [](" << wso.name << " t) -> " << wso.name
     << " {return t;})\n";

  os << "      .def(\"__deepcopy__\", [](" << wso.name << " t, py::dict&) -> "
     << wso.name << " { return t; })\n";

  os << "      .def(\"__str__\", [](" << wso.name << " t) -> std::string";
  if (wso.name == "SpeciesEnum")
    os << " { return String{toString<1>(t)}; })\n";
  else
    os << " { return String{toString(t)}; })\n";

  os << "      .def(\"__repr__\", [](" << wso.name << " t) -> std::string";
  if (wso.name == "SpeciesEnum")
    os << " { return String{toString<1>(t)}; })\n";
  else
    os << " { return String{toString(t)}; })\n";

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
        "      .def(\"__setstate__\",\n        [](" << wso.name << "* e, const std::tuple<std::string>& state) {\n"
        "           new (e) "<< wso.name << "{to<"
     << wso.name
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
    for (auto& x : value | std::views::take(value.size() - 1)) {
      if (nonstd::isdigit(x.front()) || contains_invalid_chars(x)) continue;

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

  cc << R"-x-(#include "python_interface.h"
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>

#include <enums.h>

)-x-";

  for (auto& wso : wsos) {
    cc << "NB_MAKE_OPAQUE(" << wso.name << ");\n";
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

int main() {
  enum_options("py_auto_options");
}