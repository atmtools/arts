#include <arts_options.h>
#include <nonstd.h>
#include <workspace.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <set>

#include "pydocs.h"

void enum_option(std::ostream& os, const EnumeratedOption& wso) {
os << std::format(R"-x-(
void enum_{0}(py::module_& m) {{
   py::class_<{0}> _g{0}(m, "{0}");

   _g{0}.doc() = R"-ENUMDOC-({1})-ENUMDOC-";

   xml_interface(_g{0});
   
   _g{0}.def("__str__", [](const {0}& x){{return std::format("{{}}", x);}});
   _g{0}.def("__repr__", [](const {0}& x){{return std::format("\"{{}}\"", x);}});

   _g{0}.def(py::init<>());

   _g{0}.def(py::init<{0}>());

   _g{0}.def("__init__", []({0} *y, const std::string& x){{
      new (y) {0}{{to<{0}>(x)}};
   }});
   py::implicitly_convertible<std::string, {0}>();

   _g{0}.def("__hash__", [](const {0}& x){{return std::hash<{0}>{{}}(x);}}, "Allows hashing");

   _g{0}.def("__copy__", []({0} t) -> {0}{{return t;}});

   _g{0}.def("__deepcopy__", []({0} t, py::dict&) -> {0}{{return t;}});

   _g{0}.def(py::self == py::self, "`self == other`");
   _g{0}.def(py::self != py::self, "`self != other`");
   _g{0}.def(py::self <= py::self, "`self <= other`");
   _g{0}.def(py::self >= py::self, "`self >= other`");
   _g{0}.def(py::self < py::self, "`self < other`");
   _g{0}.def(py::self > py::self, "`self > other`");

   _g{0}.def("__getstate__", []({0}& t) {{
      return std::tuple<std::string>{{String{{toString(t)}}}};
   }});

   _g{0}.def("__setstate__", []({0}* e, const std::tuple<std::string>& state) {{
      new (e) {0}{{to<{0}>(std::get<0>(state))}};
   }});

   _g{0}.def_static("get_options", [](){{return enumtyps::{0}Types;}}, "Get a list of all options");

   _g{0}.def_static("get_options_as_strings", [](){{return enumstrs::{0}Names<>;}}, "Get a list of all options as strings");

)-x-", wso.name, unwrap_stars(wso.docs()));

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

      os << "  _g" << wso.name << ".def_prop_ro_static(\"" << x;

      if (std::ranges::any_of(pykeywords, Cmp::eq(x))) {
        os << '_';
      }

      os << std::format(R"-x-(", [](py::object&){{return {}::{};}}, R"-ENUMDOC-({})-ENUMDOC-");
)-x-", wso.name, value.front(), unwrap_stars(value.back()));
    }
  }

  os << "}\n";
}

void enum_options(const std::string& fname) {
  const auto& wsos = internal_options();

  auto cc = std::ofstream(fname + ".cpp");
  auto hh = std::ofstream(fname + ".h");

  cc << R"-x-(#include "python_interface.h"
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>

#include <hpy_arts.h>

#include")-x-"
     << fname << ".h\"\n\n";

  hh << "#pragma once\n\n"
        "#include <enums.h>\n"
        "#include <nanobind/nanobind.h>\n\n";
  for (auto& wso : wsos) {
    hh << "NB_MAKE_OPAQUE(" << wso.name << ");\n";
  }

  cc << "namespace Python {\n";

  for (auto& opt : wsos) {
    enum_option(cc, opt);
  }
  cc << "void py_auto_options(py::module_& m) try {\n";
  for (auto& opt : wsos) {
    cc << "  enum_" << opt.name << "(m);\n";
  }
  cc << R"-x-(} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize automatic options\n", e.what()));
}
}  // namespace Python
)-x-";
}

int main() { enum_options("py_auto_options"); }
