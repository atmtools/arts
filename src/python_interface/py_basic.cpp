#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl_bind.h>
#include <python_interface.h>

#include <algorithm>
#include <iomanip>
#include <memory>

#include "auto_wsg.h"
#include "mystring.h"
#include "py_macros.h"
#include "python_interface_value_type.h"
#include "xml_io.h"
#include "xml_io_general_types.h"

namespace Python {
namespace Extra {
template <typename T>
void binary_ops(auto& class_) {
  static constexpr auto copycast = [](const ValueHolder<T>& x) {
    if constexpr (std::same_as<Numeric, T>) {
      return py::float_(*x.val);
    } else if constexpr (std::same_as<String, T>) {
      return py::str(*x.val);
    } else if constexpr (std::same_as<Index, T>) {
      return py::int_(*x.val);
    } else
      return py::object(*x.val);
  };

  constexpr std::array binops{
      "__add__",      "__radd__",      "__sub__",     "__rsub__",
      "__mul__",      "__rmul__",      "__div__",     "__rdiv__",
      "__matmul__",   "__rmatmul__",   "__truediv__", "__rtruediv__",
      "__floordiv__", "__rfloordiv__", "__divmod__",  "__rdivmod__",
      "__mod__",      "__rmod__",      "__pow__",     "__rpow__",
      "__lshift__",   "__rlshift__",   "__rshift__",  "__rrshift__",
      "__and__",      "__rand__",      "__xor__",     "__rxor__",
      "__or__",       "__ror__",       "__lt__",      "__le__",
      "__eq__",       "__ne__",        "__gt__",      "__ge__",
  };

  for (auto& op : binops) {
    class_.def(
        op,
        [op](const ValueHolder<T>& a, const ValueHolder<T>& b) {
          return copycast(a).attr(op)(copycast(b));
        },
        py::is_operator());

    class_.def(
        op,
        [op](const ValueHolder<T>& a, const py::object& b) {
          return copycast(a).attr(op)(b);
        },
        py::is_operator());
  }
}
}  // namespace Extra

template <typename T>
auto value_holder_artsclass(py::module_& m, const char* name, bool use_buffer) {
  auto out = [use_buffer, m, name]() {
    if (use_buffer) {
      return artsclass<ValueHolder<T>>(m, name, py::buffer_protocol())
          .def(py::init([]() { return std::make_shared<ValueHolder<T>>(); }),
               "Create default");
    }
    return artsclass<ValueHolder<T>>(m, name).def(
        py::init([]() { return std::make_shared<ValueHolder<T>>(); }),
        "Create default");
  }();

  out.def(py::init([](const ValueHolder<T>& a) {
            return std::make_shared<ValueHolder<T>>(*a.val);
          }),
          py::arg("value").none(false),
          "Copy constructor");

  if constexpr (std::convertible_to<Numeric, T>) {
    out.def(py::init([](const Numeric& a) -> ValueHolder<T> {
              return std::make_shared<T>(a);
            }),
            py::arg("value").none(false),
            "Create from value");
    py::implicitly_convertible<Numeric, ValueHolder<T>>();
  }

  if constexpr (std::convertible_to<String, T>) {
    out.def(py::init([](const String& a) -> ValueHolder<T> {
              return std::make_shared<T>(a);
            }),
            py::arg("value").none(false),
            "Create from value");
    py::implicitly_convertible<String, ValueHolder<T>>();
  }

  if constexpr (std::convertible_to<Index, T>) {
    out.def(py::init([](const Index& a) -> ValueHolder<T> {
              return std::make_shared<T>(a);
            }),
            py::arg("value").none(false),
            "Create from value");
    py::implicitly_convertible<Index, ValueHolder<T>>();
  }

  out.def("__copy__",
          [](const ValueHolder<T>& a) -> ValueHolder<T> { return *a.val; });

  out.def("__deepcopy__",
          [](const ValueHolder<T>& a, const py::dict&) -> ValueHolder<T> {
            return *a.val;
          });

  out.def("__str__",
          [](const ValueHolder<T>& a) { return var_string(*a.val); });

  if constexpr (std::same_as<String, T>) {
    out.def("__repr__",
            [](const ValueHolder<T>& a) { return py::repr(py::str(*a.val)); });
  } else {
    out.def("__repr__",
            [](const ValueHolder<T>& a) { return var_string(*a.val); });
  }

  if constexpr (std::same_as<T, String>) {
    out.def("__hash__",
            [](const ValueHolder<T>& a) { return py::hash(py::str(*a.val)); });
  } else if constexpr (std::same_as<T, Index>) {
    out.def("__hash__",
            [](const ValueHolder<T>& a) { return py::hash(py::int_(*a.val)); });
  } else if constexpr (std::same_as<T, Numeric>) {
    out.def("__hash__", [](const ValueHolder<T>& a) {
      return py::hash(py::float_(*a.val));
    });
  }

  if constexpr (std::same_as<String, T>) {
    out.def("__int__", [](const ValueHolder<T>& a) {
      return py::int_(std::stoi(*a.val));
    });
    out.def("__float__", [](const ValueHolder<T>& a) {
      return py::float_(std::stod(*a.val));
    });
    out.def("__bool__",
            [](const ValueHolder<T>& a) { return a.val->size() > 0; });
  } else {
    out.def("__int__", [](const ValueHolder<T>& a) {
      return py::int_(static_cast<Index>(*a.val));
    });
    out.def("__float__", [](const ValueHolder<T>& a) {
      return py::float_(static_cast<Numeric>(*a.val));
    });
    out.def("__complex__", [](const ValueHolder<T>& a) {
      return static_cast<Complex>(static_cast<Numeric>(*a.val));
    });
    out.def("__bool__", [](const ValueHolder<T>& a) {
      return static_cast<bool>(*a.val) and (*a.val == *a.val);
    });
  }

  if constexpr (std::same_as<String, T>) {
    out.def("__len__", [](const ValueHolder<T>& a) { return a.val->size(); });
    out.def("__getitem__", [](const ValueHolder<T>& a, const py::object& i) {
      return py::str(*a.val).attr("__getitem__")(i);
    });
  }

  out.def(
      "readxml",
      [](ValueHolder<T>& a, const char* const filename) {
        xml_read_from_file(filename, *a.val);
      },
      py::arg("file").none(false),
      py::doc(var_string("Read :class:`",
                         WorkspaceGroupInfo<T>::name,
                         "` from file\n"
                         "\n"
                         "Parameters:\n"
                         "    file (str): A file that can be read\n"
                         "\n"
                         "On Error:\n"
                         "    Throws RuntimeError for any failure to read")
                  .c_str()));

  out.def_static(
      "fromxml",
      [](const char* const filename) {
        ValueHolder<T> a{};
        xml_read_from_file(filename, *a.val);
        return a;
      },
      py::arg("file").none(false),
      py::doc(var_string("Create :class:`",
                         WorkspaceGroupInfo<T>::name,
                         "` from file\n"
                         "\n"
                         "Parameters:\n"
                         "    file (str): A file that can be read\n"
                         "\n"
                         "On Error:\n"
                         "    Throws RuntimeError for any failure to read")
                  .c_str()));

  out.def(
      "savexml",
      [](const ValueHolder<T>& a,
         const char* const filename,
         const char* const file_format,
         bool clobber) {
        xml_write_to_file(
            filename, *a.val, to<FileType>(file_format), clobber ? 0 : 1);
      },
      py::arg("file").none(false),
      py::arg("type").none(false) = "ascii",
      py::arg("clobber") = true,
      py::doc(
          var_string("Saves :class:`",
                     WorkspaceGroupInfo<T>::name,
                     "` to file\n"
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
                     "    Throws RuntimeError for any failure to save")
              .c_str()));

  out.def(py::pickle(
      [](const ValueHolder<T>& self) { return py::make_tuple(*self.val); },
      [](const py::tuple& t) -> ValueHolder<T> {
        ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
        return std::make_shared<T>(t[0].cast<T>());
      }));

  Extra::binary_ops<T>(out);

  out.doc() = var_string(WorkspaceGroupInfo<T>::desc).c_str();

  return out;
}

void py_basic(py::module_& m) try {
  value_holder_artsclass<String>(m, "String", false);

  value_holder_artsclass<Numeric>(m, "Numeric", true)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](ValueHolder<Numeric>& x) -> py::buffer_info {
        return {x.val.get(),
                sizeof(Numeric),
                py::format_descriptor<Numeric>::format(),
                0,
                {},
                {}};
      })
      .def_property(
          "value",
          py::cpp_function(
              [](ValueHolder<Numeric>& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ValueHolder<Numeric>& x, ValueHolder<Numeric>& y) {
            *x.val = *y.val;
          },
          "Operate on type as if :class:`numpy.ndarray` type");

  value_holder_artsclass<Index>(m, "Index", true)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](ValueHolder<Index>& x) -> py::buffer_info {
        return {x.val.get(),
                sizeof(Index),
                py::format_descriptor<Index>::format(),
                0,
                {},
                {}};
      })
      .def_property(
          "value",
          py::cpp_function(
              [](ValueHolder<Index>& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ValueHolder<Index>& x, ValueHolder<Index>& y) { *x.val = *y.val; },
          "Operate on type as if :class:`numpy.ndarray` type");

  artsarray<ArrayOfString, ArrayOptions::index>(m, "ArrayOfString")
      .def(py::init(
          [](const std::vector<std::variant<ValueHolder<String>, String>>& x) {
            auto out = std::make_shared<ArrayOfString>(x.size());
            std::transform(x.begin(), x.end(), out->begin(), [](auto& s) {
              return std::visit([](auto& v) -> const String& { return v; }, s);
            });
            return out;
          }))
      .PythonInterfaceFileIO(ArrayOfString)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfString);
  py::implicitly_convertible<
      std::vector<std::variant<ValueHolder<String>, String>>,
      ArrayOfString>();

  artsarray<ArrayOfIndex, ArrayOptions::index>(
      m, "ArrayOfIndex", py::buffer_protocol())
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_property(
          "value",
          py::cpp_function(
              [](ArrayOfIndex& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ArrayOfIndex& x, ArrayOfIndex& y) { x = y; },
          "Operate on type as if :class:`numpy.ndarray` type")
      .PythonInterfaceFileIO(ArrayOfIndex)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfIndex);

  artsarray<ArrayOfNumeric>(m, "ArrayOfNumeric", py::buffer_protocol())
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_property(
          "value",
          py::cpp_function(
              [](ArrayOfNumeric& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ArrayOfNumeric& x, ArrayOfNumeric& y) { x = y; },
          "Operate on type as if :class:`numpy.ndarray` type")
      .doc() = "List of :class:`~Numeric` values";

  artsarray<ArrayOfArrayOfString>(m, "ArrayOfArrayOfString")
      .PythonInterfaceFileIO(ArrayOfArrayOfString)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfString);

  artsarray<ArrayOfArrayOfIndex>(m, "ArrayOfArrayOfIndex")
      .PythonInterfaceFileIO(ArrayOfArrayOfIndex)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfIndex);

  artsclass<Any>(m, "Any")
      .def(py::init([]() { return std::make_shared<Any>(); }), "Create empty")
      .def(py::init([](const py::args&, const py::kwargs&) {
             return std::make_shared<Any>();
           }),
           "Create empty")
      .PythonInterfaceBasicRepresentation(Any)
      .PythonInterfaceFileIO(Any)
      .def(py::pickle([](const py::object&) { return py::make_tuple(); },
                      [](const py::tuple& t) {
                        ARTS_USER_ERROR_IF(t.size() != 0, "Invalid state!")
                        return std::make_shared<Any>();
                      }))
      .PythonInterfaceWorkspaceDocumentation(Any);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize basic\n", e.what()));
}
}  // namespace Python
