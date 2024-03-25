#include <py_auto_wsg_init.h>
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
auto& fix_value_holder_artsclass(artsclass<ValueHolder<T>>& out) {
  if constexpr (std::convertible_to<Numeric, T>) {
    py::implicitly_convertible<Numeric, ValueHolder<T>>();
  }

  if constexpr (std::convertible_to<String, T>) {
    py::implicitly_convertible<String, ValueHolder<T>>();
  }

  if constexpr (std::convertible_to<Index, T>) {
    py::implicitly_convertible<Index, ValueHolder<T>>();
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
  fix_value_holder_artsclass(py_staticString(m));

  fix_value_holder_artsclass(py_staticNumeric(m))
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

  fix_value_holder_artsclass(py_staticIndex(m))
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

  py_staticArrayOfString(m)
      .def("index",
           [](const ArrayOfString& v, const String& f) {
             auto ptr = std::ranges::find(v, f);
             if (ptr == v.end())
               throw py::value_error(var_string(f, " is not in list"));
             return std::distance(v.begin(), ptr);
           })
      .def(py::init(
          [](const std::vector<std::variant<ValueHolder<String>, String>>& x) {
            auto out = std::make_shared<ArrayOfString>(x.size());
            std::transform(x.begin(), x.end(), out->begin(), [](auto& s) {
              return std::visit([](auto& v) -> const String& { return v; }, s);
            });
            return out;
          }));
  py::implicitly_convertible<
      std::vector<std::variant<ValueHolder<String>, String>>,
      ArrayOfString>();

  py_staticArrayOfIndex(m)
      .def("index",
           [](const ArrayOfIndex& v, const Index& f) {
             auto ptr = std::ranges::find(v, f);
             if (ptr == v.end())
               throw py::value_error(var_string(f, " is not in list"));
             return std::distance(v.begin(), ptr);
           })
      .def_buffer([](ArrayOfIndex& x) -> py::buffer_info {
        return py::buffer_info(x.data(),
                               sizeof(Index),
                               py::format_descriptor<Index>::format(),
                               1,
                               {static_cast<ssize_t>(x.size())},
                               {sizeof(Index)});
      })
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
          "Operate on type as if :class:`numpy.ndarray` type");

  artsarray<ArrayOfNumeric>(m, "ArrayOfNumeric", py::buffer_protocol())
      .def("index",
           [](const ArrayOfNumeric& v, const Numeric& f) {
             auto ptr = std::ranges::find(v, f);
             if (ptr == v.end())
               throw py::value_error(var_string(f, " is not in list"));
             return std::distance(v.begin(), ptr);
           })
      .def_buffer([](ArrayOfNumeric& x) -> py::buffer_info {
        return py::buffer_info(x.data(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {static_cast<ssize_t>(x.size())},
                               {sizeof(Numeric)});
      })
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

  py_staticArrayOfArrayOfString(m);

  py_staticArrayOfArrayOfIndex(m);

  py_staticAny(m)
      .def(py::init([](const py::args&, const py::kwargs&) {
             return std::make_shared<Any>();
           }),
           "Create empty")
      .def(py::pickle([](const py::object&) { return py::make_tuple(); },
                      [](const py::tuple& t) {
                        ARTS_USER_ERROR_IF(t.size() != 0, "Invalid state!")
                        return std::make_shared<Any>();
                      }));
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize basic\n", e.what()));
}
}  // namespace Python
