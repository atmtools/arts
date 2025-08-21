#pragma once

#include <lagrange_interp.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <workspace.h>

#include <iterator>

#include "python_interface_value_type.h"

NB_MAKE_OPAQUE(
    std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>);
NB_MAKE_OPAQUE(std::vector<lbl::line_shape::species_model>);
NB_MAKE_OPAQUE(std::vector<lbl::line>)
NB_MAKE_OPAQUE(Array<lagrange_interp::lag_t<-1, lagrange_interp::identity>>);
NB_MAKE_OPAQUE(Array<lagrange_interp::lag_t<-1, lagrange_interp::loncross>>);
NB_MAKE_OPAQUE(Array<Array<SpeciesTag>>);

namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename T, typename... Ts>
void vector_interface(py::class_<Array<T>, Ts...> &c) {
  using Vec = Array<T>;

  c.def("__getstate__", [](const Vec &v) { return std::tuple<Vec>{v}; });

  c.def("__setstate__",
        [](Vec &v, const std::tuple<Vec> &x) { new (&v) Vec{std::get<0>(x)}; });
}

// THIS IS A DIRECT COPY OF THE FUNCTION FROM bind_vector.h
// SOME MODIFICATIONS ARE DONE, SPECIFICALLY TO GET THE INIT
// FUNCTION TO WORK WITH THE ARRAY OF STRING USING OUR VALUEHOLDER
// IF THIS FAILS, WE MUST PROBABLY REVERT TO THE ORIGINAL FUNCTION
// OR REWRITE THIS ENTIRERLY.
// IT ALSO WRAPS vector_interface FOR SIMPLICITY
template <typename Value>
void value_holder_vector_interface(py::class_<Array<Value>> &cl) {
  using namespace nanobind;
  using ValueRef = Value &;

  cl.def(init<>(), "Default constructor")

      .def(
          "__len__",
          [](const Array<Value> &v) { return v.size(); },
          "Return the number of items in the list.")

      .def(
          "__bool__",
          [](const Array<Value> &v) { return !v.empty(); },
          "Check whether the vector is nonempty")

      .def(
          "__iter__",
          [](Array<Value> &v) {
            return make_iterator(
                type<Array<Value>>(), "Iterator", v.begin(), v.end());
          },
          keep_alive<0, 1>(),
          "Allows `iter(self)`")

      .def(
          "__getitem__",
          [](Array<Value> &v, Py_ssize_t i) -> ValueRef {
            return v[detail::wrap(i, v.size())];
          },
          rv_policy::automatic_reference,
          "i"_a,
          "Return the i-th item.")

      .def(
          "clear",
          [](Array<Value> &v) { v.clear(); },
          "Remove all items from list.");

  if constexpr (detail::is_copy_constructible_v<Value>) {
    cl.def(init<const Array<Value> &>(), "Copy constructor");

    cl.def(
        "__init__",
        [](Array<Value> *v, typed<iterable, ValueHolder<Value>> seq) {
          new (v) Array<Value>();
          v->reserve(len_hint(seq));
          for (handle h : seq) v->push_back(cast<ValueHolder<Value>>(h));
        },
        "Construct from an iterable object");

    implicitly_convertible<iterable, Array<Value>>();

    cl.def(
          "append",
          [](Array<Value> &v, const Value &value) { v.push_back(value); },
          "value"_a,
          "Append `arg` to the end of the list.")

        .def(
            "insert",
            [](Array<Value> &v, Py_ssize_t i, const Value &x) {
              if (i < 0) i += (Py_ssize_t)v.size();
              if (i < 0 || (size_t)i > v.size()) throw index_error();
              v.insert(v.begin() + i, x);
            },
            "Insert object `arg1` before index `arg0`.")

        .def(
            "pop",
            [](Array<Value> &v, Py_ssize_t i) {
              size_t index = detail::wrap(i, v.size());
              Value result = std::move(v[index]);
              v.erase(v.begin() + index);
              return result;
            },
            arg("index") = -1,
            "Remove and return item at `index` (default last).")

        .def(
            "extend",
            [](Array<Value> &v, const Array<Value> &src) {
              v.insert(v.end(), src.begin(), src.end());
            },
            "Extend `self` by appending elements from `arg`.")

        .def(
            "__setitem__",
            [](Array<Value> &v, Py_ssize_t i, const Value &value) {
              v[detail::wrap(i, v.size())] = value;
            },
            "i"_a,
            "value"_a,
            "Set the i-th item to `value`.")

        .def(
            "__delitem__",
            [](Array<Value> &v, Py_ssize_t i) {
              v.erase(v.begin() + detail::wrap(i, v.size()));
            },
            "i"_a,
            "Remove the i-th item from the list.")

        .def(
            "__getitem__",
            [](const Array<Value> &v, const slice &slice) -> Array<Value> * {
              auto [start, stop, step, length] = slice.compute(v.size());
              auto *seq                        = new Array<Value>();
              seq->reserve(length);

              for (size_t i = 0; i < length; ++i) {
                seq->push_back(v[start]);
                start += step;
              }

              return seq;
            },
            "slice"_a,
            "Return a slice of the list.")

        .def(
            "__setitem__",
            [](Array<Value> &v, const slice &slice, const Array<Value> &value) {
              auto [start, stop, step, length] = slice.compute(v.size());

              if (length != value.size())
                throw index_error(
                    "The left and right hand side of the slice "
                    "assignment have mismatched sizes!");

              for (size_t i = 0; i < length; ++i) {
                v[start]  = value[i];
                start    += step;
              }
            },
            "slice"_a,
            "value"_a,
            "Set a slice of the list.")

        .def(
            "__delitem__",
            [](Array<Value> &v, const slice &slice) {
              auto [start, stop, step, length] = slice.compute(v.size());
              if (length == 0) return;

              stop = start + (length - 1) * step;
              if (start > stop) {
                std::swap(start, stop);
                step = -step;
              }

              if (step == 1) {
                v.erase(v.begin() + start, v.begin() + stop + 1);
              } else {
                for (size_t i = 0; i < length; ++i) {
                  v.erase(v.begin() + stop);
                  stop -= step;
                }
              }
            },
            "slice"_a,
            "Remove a slice of the list.");
  }

  if constexpr (detail::is_equality_comparable_v<Value>) {
    cl.def(self == self, "`self == other`");
    cl.def(self != self, "`self != other`");

    cl.def(
          "__contains__",
          [](const Array<Value> &v, const Value &x) {
            return std::find(v.begin(), v.end(), x) != v.end();
          },
          "Return true if `arg` is in the list.")

        .def(
            "__contains__",  // fallback for incompatible types
            [](const Array<Value> &, handle) { return false; },
            "Return false.")

        .def(
            "count",
            [](const Array<Value> &v, const Value &x) {
              return std::count(v.begin(), v.end(), x);
            },
            "Return number of occurrences of `arg`.")

        .def(
            "index",
            [](const Array<Value> &v,
               const Value &value,
               Size start,
               Size end) {
              auto s = v.begin() + std::min(start, v.size());
              auto e = v.begin() + std::min(end, v.size());
              auto p = std::find(s, e, value);
              if (p == e)
                throw std::invalid_argument(
                    std::format("{} is not in list", value));
              return Index{std::distance(v.begin(), p)};
            },
            "value"_a,
            "start"_a = Size{0},
            "end"_a   = std::numeric_limits<Size>::max(),
            "Return first occurence of value between `start` and `end` - defaults to full range")

        .def(
            "remove",
            [](Array<Value> &v, const Value &x) {
              auto p = std::find(v.begin(), v.end(), x);
              if (p != v.end())
                v.erase(p);
              else
                throw value_error();
            },
            "Remove first occurrence of `arg`.");
  }

  vector_interface(cl);
}
}  // namespace Python
