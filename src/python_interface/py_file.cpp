#include <enums-common-helper.h>
#include <enumsFileType.h>
#include <file.h>
#include <isotopologues.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <parameters.h>
#include <py_auto_options.h>

#include <optional>

#include "hpy_arts.h"
#include "python_interface.h"

extern Parameters parameters;

namespace Python {
namespace py = nanobind;
using namespace py::literals;

void enum_FileType(py::module_& m) {
  py::class_<FileType> _gFileType(m, "FileType");
  _gFileType.doc() = PythonWorkspaceGroupInfo<FileType>::desc();
  xml_interface(_gFileType);
  _gFileType.def("__str__", [](const FileType& x) { return std::format("{}", x); });
  _gFileType.def("__repr__", [](const FileType& x) { return std::format("\"{}\"", x); });
  _gFileType.def(py::init<>());
  _gFileType.def(py::init<FileType>());
  _gFileType.def("__init__", [](FileType* y, const std::string& x) { new (y) FileType{to<FileType>(x)}; });
  py::implicitly_convertible<std::string, FileType>();
  _gFileType.def("__hash__", [](const FileType& x) { return std::hash<FileType>{}(x); }, "Allows hashing");
  _gFileType.def("__copy__", [](FileType t) -> FileType { return t; });
  _gFileType.def("__deepcopy__", [](FileType t, py::dict&) -> FileType { return t; });
  _gFileType.def(py::self == py::self, "`self == other`");
  _gFileType.def(py::self != py::self, "`self != other`");
  _gFileType.def(py::self <= py::self, "`self <= other`");
  _gFileType.def(py::self >= py::self, "`self >= other`");
  _gFileType.def(py::self < py::self, "`self < other`");
  _gFileType.def(py::self > py::self, "`self > other`");
  _gFileType.def_static("get_options", [] { return enumtyps::FileTypeTypes; }, "Get a list of all options");
  _gFileType.def_static(
      "get_options_as_strings",
      [](Size i) {
        if (i == 0) return enumstrs::FileTypeNames<0>;
        if (i == 1) return enumstrs::FileTypeNames<1>;
        if (i == 2) return enumstrs::FileTypeNames<2>;
        if (i == 3) return enumstrs::FileTypeNames<3>;
        throw std::invalid_argument("Valid input [0, 4)");
      },
      "i"_a = Size{0},
      "Get a list of all options as strings");
  _gFileType.def_prop_ro_static(
      "ascii",
      [](py::object&) { return FileType::ascii; },
      R"-ENUMDOC-(Save as ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "ASCII",
      [](py::object&) { return FileType::ascii; },
      R"-ENUMDOC-(Save as ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "Ascii",
      [](py::object&) { return FileType::ascii; },
      R"-ENUMDOC-(Save as ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "text",
      [](py::object&) { return FileType::ascii; },
      R"-ENUMDOC-(Save as ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "zascii",
      [](py::object&) { return FileType::zascii; },
      R"-ENUMDOC-(Save as zipped ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "ZASCII",
      [](py::object&) { return FileType::zascii; },
      R"-ENUMDOC-(Save as zipped ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "Zip",
      [](py::object&) { return FileType::zascii; },
      R"-ENUMDOC-(Save as zipped ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "zip",
      [](py::object&) { return FileType::zascii; },
      R"-ENUMDOC-(Save as zipped ASCII
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "binary",
      [](py::object&) { return FileType::binary; },
      R"-ENUMDOC-(Save as binary data
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "BINARY",
      [](py::object&) { return FileType::binary; },
      R"-ENUMDOC-(Save as binary data
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "Binary",
      [](py::object&) { return FileType::binary; },
      R"-ENUMDOC-(Save as binary data
)-ENUMDOC-");
  _gFileType.def_prop_ro_static(
      "bin",
      [](py::object&) { return FileType::binary; },
      R"-ENUMDOC-(Save as binary data
)-ENUMDOC-");
}

void py_file(py::module_& m) try {
  enum_FileType(m);

  auto file =
      m.def_submodule("file", R"--(Contain methods to handle files available to ARTS via its path environment.)--");

  file.def(
      "find",
      [](const String& filename) -> std::optional<String> {
        ArrayOfString files;
        if (find_file2(files, filename, parameters.datapath)) return files[0];
        return std::nullopt;
      },
      R"(Find file in paths available to ARTS.

Parameters
----------
filename : str
    The name of the file to find minus ".xml".

Returns
-------
str or None
    The full path to the file if found, otherwise None.
)",
      "filename"_a);

  file.def(
      "find_xml",
      [](const String& filename) -> std::optional<String> {
        ArrayOfString files;
        if (find_file(files, filename, parameters.datapath, {".xml"})) return files[0];
        return std::nullopt;
      },
      R"(Find xml-file in paths available to ARTS.

Parameters
----------
filename : str
    The name of the file to find.

Returns
-------
str or None
    The full path to the file if found, otherwise None.
)",
      "filename"_a);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize files\n{}", e.what()));
}
}  // namespace Python
