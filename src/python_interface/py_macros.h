#ifndef py_macros_h
#define py_macros_h

#include <python_interface.h>
#include <xml_io.h>

#include <functional>

namespace Python {
namespace py = nanobind;
}  // namespace Python

constexpr Index negative_clamp(const Index i, const Index n) noexcept {
  return (i < 0) ? i + n : i;
}

#define PythonInterfaceFileIO2(PythonType, RealType)                    \
  def(                                                                  \
      "savexml",                                                        \
      [](const PythonType& x,                                           \
         const char* const file,                                        \
         const char* const type,                                        \
         bool clobber) {                                                \
        xml_write_to_file(file,                                         \
                          static_cast<const RealType&>(x),              \
                          to<FileType>(type),                           \
                          clobber ? 0 : 1);                             \
      },                                                                \
      py::arg("file").none(false),                                      \
      py::arg("type").none(false) = "ascii",                            \
      py::arg("clobber")          = true,                               \
      py::doc("Saves :class:`" #RealType "` to file\n"                  \
              "\n"                                                      \
              "Parameters:\n"                                           \
              "    file (str): The path to which the file is written."  \
              " Note that several of the options might modify the"      \
              " name or write more files\n"                             \
              "    type (str): Type of file to save (ascii. zascii,"    \
              " or binary)\n"                                           \
              "    clobber (bool): Overwrite existing files or add new" \
              " file with modified name?\n"                             \
              "\n"                                                      \
              "On Error:\n"                                             \
              "    Throws RuntimeError for any failure to save"))       \
      .def(                                                             \
          "readxml",                                                    \
          [](PythonType& x, const char* const file) {                   \
            xml_read_from_file(file, x);                                \
          },                                                            \
          py::arg("file").none(false),                                  \
          py::doc("Read :class:`" #RealType "` from file\n"             \
                  "\n"                                                  \
                  "Parameters:\n"                                       \
                  "    file (str): A file that can be read\n"           \
                  "\n"                                                  \
                  "On Error:\n"                                         \
                  "    Throws RuntimeError for any failure to read"))   \
      .def_static(                                                      \
          "fromxml",                                                    \
          [](const char* const file) -> PythonType {                    \
            RealType x;                                                 \
            xml_read_from_file(file, x);                                \
            return x;                                                   \
          },                                                            \
          py::arg("file").none(false),                                  \
          py::doc("Create :class:`" #RealType "` from file\n"           \
                  "\n"                                                  \
                  "Parameters:\n"                                       \
                  "    file (str): A file that can be read\n"           \
                  "\n"                                                  \
                  "On Error:\n"                                         \
                  "    Throws RuntimeError for any failure to read"))

#define PythonInterfaceFileIO(Type) PythonInterfaceFileIO2(Type, Type)

#define PythonInterfaceIndexItemAccess(Type)                                 \
  def("__len__", [](const Type& x) { return x.size(); })                     \
      .def(                                                                  \
          "__getitem__",                                                     \
          [](Type& x, Index i)                                               \
              -> std::shared_ptr<std::remove_cvref_t<decltype(x[i])>> {      \
            i = negative_clamp(i, x.size());                                 \
            if (x.size() <= static_cast<Size>(i) or i < 0)                   \
              throw std::out_of_range(var_string("Bad index access: ",       \
                                                 i,                          \
                                                 " in object of range [0, ", \
                                                 x.size(),                   \
                                                 ")"));                      \
            return std::shared_ptr<std::remove_cvref_t<decltype(x[i])>>(     \
                &x[i], [](void*) {});                                        \
          },                                                                 \
          py::rv_policy::reference_internal,                                 \
          py::keep_alive<0, 1>())                                            \
      .def("__setitem__", [](Type& x, Index i, decltype(x[i]) y) {           \
        i = negative_clamp(i, x.size());                                     \
        if (x.size() <= static_cast<Size>(i) or i < 0)                       \
          throw std::out_of_range(var_string("Bad index access: ",           \
                                             i,                              \
                                             " in object of range [0, ",     \
                                             x.size(),                       \
                                             ")"));                          \
        x[i] = std::move(y);                                                 \
      })

#define PythonInterfaceBasicRepresentation(Type) \
  def("__str__", [](const Type& x) {             \
    return var_string(x);                        \
  }).def("__repr__", [](const Type& x) { return var_string(x); })

#define PythonInterfaceCopyValue(Type)  \
  def("__copy__", [](Type& t) -> Type { \
    return t;                           \
  }).def("__deepcopy__", [](Type& t, py::dict&) -> Type { return t; })

#endif
