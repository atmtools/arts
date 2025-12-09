#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <py_auto_wsgdocs.h>
#include <xml_io.h>

#include <concepts>
#include <format>

#include "python_interface_value_type.h"

namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename T, typename U = T>
void xml_interface(py::class_<T>& c) {
  c.def(
      "savexml",
      [](const T& x,
         const char* const file,
         const char* const type,
         bool clobber) {
        return xml_write_to_file(file,
                                 static_cast<const U&>(x),
                                 to<FileType>(type),
                                 clobber ? 0 : 1);
      },
      "file"_a.none(false),
      "type"_a.none(false) = "ascii",
      "clobber"_a          = true,
      R"(Saves variable to file.

Parameters
----------
file : str
    The path to which the file is written. Note that several of the options might modify the name or write more files.
type : str, optional
    Type of file to save.  See :class:`FileType` for options.  Defaults is "ascii".
clobber : bool, optional
    Overwrite existing files or add new file with modified name?  Defaults is True.

Raises
------
  RuntimeError
      For any failure to write.

Return
------
file : str
    The file saved.  May differ from input.
)");

  c.def(
      "readxml",
      [](T& x, const char* const file) {
        return xml_read_from_file(file, static_cast<U&>(x));
      },
      "file"_a.none(false),
      R"(Read variable from file.

Parameters
----------
file : str
    A file that can be read.

Raises
------
  RuntimeError
      For any failure to read.

Return
------
file : str
    The file path found (may differ from input due to environment variables).
)");

  if constexpr (arts_xml_extendable<U>) {
    c.def(
        "extendxml",
        [](T& x, const char* const file) {
          return xml_extend_from_file(file, static_cast<U&>(x));
        },
        "file"_a.none(false),
        R"(Extend variable from file.

The content of the file is added to the existing variable.

Parameters
----------
file : str
    A file that can be read.

Raises
------
  RuntimeError
      For any failure to read.

Return
------
file : str
    The file path found (may differ from input due to environment variables).
)");
  }

  if constexpr (arts_xml_appendable<U>) {
    c.def(
        "appendxml",
        [](T& x, const char* const file) {
          return xml_append_from_file(file, static_cast<U&>(x));
        },
        "file"_a.none(false),
        R"(Append variable from file.

The content of the file is added to the existing variable.

Parameters
----------
file : str
    A file that can be read.

Raises
------
  RuntimeError
      For any failure to read.

Return
------
file : str
    The file path found (may differ from input due to environment variables).
)");
  }

  c.def_static(
      "fromxml",
      [](const char* const file) -> T {
        U x;
        xml_read_from_file(file, x);
        return x;
      },
      "file"_a.none(false),
      R"(Create variable from file.

Parameters
----------
file : str
    A file that can be read

Raises
------
  RuntimeError
      For any failure to read.

Return
------
artstype : T
    The variable created from the file.
)");
}

static constexpr std::array binops{
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

template <typename T>
constexpr auto copycast(const ValueHolder<T>& x) {
  if constexpr (std::same_as<Numeric, T>) {
    return py::float_(*x.val);
  } else if constexpr (std::same_as<String, T>) {
    return py::str(String{x}.c_str());
  } else if constexpr (std::same_as<Index, T>) {
    return py::int_(*x.val);
  } else
    return py::object(*x.val);
};

template <typename T>
void value_holder_interface(py::class_<ValueHolder<T>>& c) {
  for (auto& op : binops) {
    c.def(
        op,
        [op](const ValueHolder<T>& a, const ValueHolder<T>& b) {
          return copycast(a).attr(op)(copycast(b));
        },
        py::is_operator());

    c.def(
        op,
        [op](const ValueHolder<T>& a, const py::object& b) {
          return copycast(a).attr(op)(b);
        },
        py::is_operator());
  }

  if constexpr (std::same_as<String, T>) {
    c.def(
        "__int__",
        [](const ValueHolder<T>& a) { return py::int_(std::stoi(*a.val)); },
        "Allows conversion to int");
    c.def(
        "__float__",
        [](const ValueHolder<T>& a) { return py::float_(std::stod(*a.val)); },
        "Allows conversion to float");
    c.def(
        "__bool__",
        [](const ValueHolder<T>& a) { return a.val->size() > 0; },
        "Allows conversion to bool");
  } else {
    c.def(
        "__int__",
        [](const ValueHolder<T>& a) {
          return py::int_(static_cast<Index>(*a.val));
        },
        "Allows conversion to int");
    c.def(
        "__float__",
        [](const ValueHolder<T>& a) {
          return py::float_(static_cast<Numeric>(*a.val));
        },
        "Allows conversion to float");
    c.def("__complex__", [](const ValueHolder<T>& a) {
      return static_cast<Complex>(static_cast<Numeric>(*a.val));
    });
    c.def(
        "__bool__",
        [](const ValueHolder<T>& a) {
          return static_cast<bool>(*a.val) and (*a.val == *a.val);
        },
        "Allows conversion to bool");
  }

  c.def(
      "__hash__",
      [](const ValueHolder<T>& a) { return std::hash<T>{}(*a.val); },
      "Allows hashing");

  if constexpr (std::same_as<String, T>) {
    c.def(
        "__len__",
        [](const ValueHolder<T>& a) { return a.val->size(); },
        "Length of the string");
    c.def(
        "__getitem__",
        [](const ValueHolder<T>& a, const py::object& i) {
          return py::str(String{a}.c_str()).attr("__getitem__")(i);
        },
        "Get item from string");
  }
}

template <typename T>
void str_interface(py::class_<T>& c) {
  c.def("__format__", [](const T& x, std::string fmt) {
    if constexpr (std::formattable<T, char>) {
      fmt = std::format("{}{}{}", "{:"sv, fmt, "}"sv);
      return std::vformat(fmt, std::make_format_args(x));
    } else {
      if (not fmt.empty()) throw std::format_error("Not formattable");
      return std::format("{}", x);
    }
  });

  c.def("__str__", [](const T& x) {
    if constexpr (std::formattable<T, char>) {
      return std::format("{:B,}", x);
    } else {
      return std::format("{}", x);
    }
  });

  c.def("__repr__", [](const T& x) {
    if constexpr (std::formattable<T, char>) {
      return std::format("{:sB,}", x);
    } else {
      return std::format("{}", x);
    }
  });
}

template <typename T, class... E>
void boolean_compare(py::class_<T, E...>& c [[maybe_unused]]) {
  if constexpr (std::three_way_comparable<T>) {
    c.def(py::self == py::self);
    c.def(py::self != py::self);
    c.def(py::self <= py::self);
    c.def(py::self >= py::self);
    c.def(py::self < py::self);
    c.def(py::self > py::self);
  }
}

template <WorkspaceGroup T, class... E>
void generic_interface(py::class_<ValueHolder<T>, E...>& c) {
  using U = ValueHolder<T>;

  c.def(py::init<>());
  c.def(py::init_implicit<T>());
  c.def(py::init<U>());

  c.def("__copy__", [](const U& t) -> U { return t; });
  c.def("__deepcopy__", [](const U& t, py::dict&) -> U { return t; });

  str_interface(c);
  xml_interface<U, T>(c);

  c.doc() = std::string{PythonWorkspaceGroupInfo<T>::desc()};
}

template <typename T, class... E>
void generic_interface(py::class_<T, E...>& c) {
  if constexpr (std::is_default_constructible_v<T>) c.def(py::init<>());

  if constexpr (std::is_copy_constructible_v<T>) {
    c.def(py::init<T>());
    c.def("__copy__", [](const T& t) -> T { return t; });
    c.def("__deepcopy__", [](const T& t, py::dict&) -> T { return t; });
  }

  if constexpr (arts_formattable_or_value_type<T>) str_interface(c);

  if constexpr (arts_xml_ioable<T>) xml_interface(c);

  boolean_compare(c);

  if constexpr (requires { PythonWorkspaceGroupInfo<T>::desc(); }) {
    c.doc() = std::string{PythonWorkspaceGroupInfo<T>::desc()};
  }
}
}  // namespace Python
