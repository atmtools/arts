#include <file.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <parameters.h>
#include <isotopologues.h>
#include <py_auto_options.h>

#include <optional>

extern Parameters parameters;

namespace Python {
namespace py = nanobind;
using namespace py::literals;

void py_file(py::module_& m) try {
  auto file = m.def_submodule(
      "file",
      R"--(Contain methods to handle files available to ARTS via its path environment.)--");

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
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize constant\n", e.what()));
}
}  // namespace Python
