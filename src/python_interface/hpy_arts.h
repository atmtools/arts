#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <py_auto_wsgdocs.h>
#include <xml_io.h>

#include <concepts>
#include <format>
#include <iomanip>

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
        xml_write_to_file(file,
                          static_cast<const U&>(x),
                          to<FileType>(type),
                          clobber ? 0 : 1);
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
      [](T& x, const char* const file) {
        xml_read_from_file(file, static_cast<U&>(x));
      },
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
    c.def("__int__",
          [](const ValueHolder<T>& a) { return py::int_(std::stoi(*a.val)); });
    c.def("__float__", [](const ValueHolder<T>& a) {
      return py::float_(std::stod(*a.val));
    });
    c.def("__bool__",
          [](const ValueHolder<T>& a) { return a.val->size() > 0; });
  } else {
    c.def("__int__", [](const ValueHolder<T>& a) {
      return py::int_(static_cast<Index>(*a.val));
    });
    c.def("__float__", [](const ValueHolder<T>& a) {
      return py::float_(static_cast<Numeric>(*a.val));
    });
    c.def("__complex__", [](const ValueHolder<T>& a) {
      return static_cast<Complex>(static_cast<Numeric>(*a.val));
    });
    c.def("__bool__", [](const ValueHolder<T>& a) {
      return static_cast<bool>(*a.val) and (*a.val == *a.val);
    });
  }

  c.def("__hash__",
        [](const ValueHolder<T>& a) { return std::hash<T>{}(*a.val); });

  c.def("__getstate__",
        [](const ValueHolder<T>& self) { return std::make_tuple(*self.val); });
  c.def("__setstate__", [](ValueHolder<T>* self, const std::tuple<T>& state) {
    new (self) ValueHolder<T>{std::get<0>(state)};
  });

  if constexpr (std::same_as<String, T>) {
    c.def("__len__", [](const ValueHolder<T>& a) { return a.val->size(); });
    c.def("__getitem__", [](const ValueHolder<T>& a, const py::object& i) {
      return py::str(String{a}.c_str()).attr("__getitem__")(i);
    });
  }
}

template <typename T>
void str_interface(py::class_<T>& c) {
  c.def("__format__", [](const T& x, std::string fmt) {
    if constexpr (std::formattable<T, char>) {
      fmt = std::format("{}{}{}", "{:"sv, fmt, "}"sv);
      return std::vformat(fmt.c_str(), std::make_format_args(x));
    } else {
      if (not fmt.empty()) {
        throw std::format_error(
            var_string("Cannot support options: ", '"', fmt, '"'));
      }

      return var_string(x);
    }
  });

  c.def("__str__", [](const T& x) {
    if constexpr (std::formattable<T, char>) {
      return std::format("{:qNB,}", x);
    } else {
      return var_string(x);
    }
  });

  c.def("__repr__", [](const T& x) {
    if constexpr (std::formattable<T, char>) {
      return std::format("{:sqNB,}", x);
    } else {
      return var_string(x);
    }
  });
}

template <WorkspaceGroup T>
void workspace_group_interface(py::class_<T>& c) {
  c.def(py::init<>());
  c.def(py::init<T>());

  c.def("__copy__", [](const T& t) -> T { return t; });
  c.def("__deepcopy__", [](const T& t, py::dict&) -> T { return t; });

  str_interface(c);
  xml_interface<T>(c);

  c.doc() = PythonWorkspaceGroupInfo<T>::desc;
}

template <WorkspaceGroup T>
void workspace_group_interface(py::class_<ValueHolder<T>>& c) {
  using U = ValueHolder<T>;

  c.def(py::init<>());
  c.def(py::init_implicit<T>());
  c.def(py::init<U>());

  c.def("__copy__", [](const U& t) -> U { return t; });
  c.def("__deepcopy__", [](const U& t, py::dict&) -> U { return t; });

  str_interface(c);
  xml_interface<U, T>(c);

  c.doc() = PythonWorkspaceGroupInfo<T>::desc;
}
}  // namespace Python
