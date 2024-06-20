#pragma once

#include <nanobind/nanobind.h>
#include <xml_io.h>

namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename T, typename U = T>
void xml_interface(py::class_<T>& c) {
  c.def(
      "savexml",
      [](const U& x,
         const char* const file,
         const char* const type,
         bool clobber) {
        xml_write_to_file(file, x, to<FileType>(type), clobber ? 0 : 1);
      },
      "file"_a.none(false),
      "type"_a.none(false) = "ascii",
      "clobber"_a          = true,
      "Saves variable to file\n"
      "\n"
      "Parameters:\n"
      "    file (str): The path to which the file is written."
      " Note that several of the options might modify the"
      " name or write more files\n"
      "    type (str): Type of file to save (ascii. zascii,"
      " or binary)\n"
      "    clobber (bool): Overwrite existing files or add new"
      " file with modified name?\n"
      "\n"
      "On Error:\n"
      "    Throws RuntimeError for any failure to save");

  c.def(
      "readxml",
      [](U& x, const char* const file) { xml_read_from_file(file, x); },
      "file"_a.none(false),
      "Read variable from file\n"
      "\n"
      "Parameters:\n"
      "    file (str): A file that can be read\n"
      "\n"
      "On Error:\n"
      "    Throws RuntimeError for any failure to read");

  c.def_static(
      "fromxml",
      [](const char* const file) -> T {
        U x;
        xml_read_from_file(file, x);
        return x;
      },
      "file"_a.none(false),
      "Create variable from file\n"
      "\n"
      "Parameters:\n"
      "    file (str): A file that can be read\n"
      "\n"
      "On Error:\n"
      "    Throws RuntimeError for any failure to read");
}
}  // namespace Python
